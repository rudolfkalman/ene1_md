"""Raster finite-volume copper sheet model. Python 3 + Pillow; solve.cpp does IC-PCG.
No original design files are modified. All lengths in geometry JSON are mm.
"""
from PIL import Image,ImageDraw,ImageChops,ImageFont
from pathlib import Path
from array import array
import json,math,sys,struct,subprocess,collections
ROOT=Path(__file__).resolve().parent
X0,Y0,X1,Y1=97.5,55.5,170,148
RHO=1.724e-8
THICK=[35e-6,17e-6,17e-6,35e-6]
Z=[0,.126,1.383,1.509] # conductor midplanes, mm
class Raster:
 def __init__(self,geom,step):
  self.g=geom;self.step=step;self.w=round((X1-X0)/step);self.h=round((Y1-Y0)/step);self.n=self.w*self.h
 def draw(self,polys):
  # 4x sampling with area-majority classification minimizes edge inflation.
  s=4;w,h=self.w*s,self.h*s
  im=Image.new('L',(w,h));d=ImageDraw.Draw(im)
  def coords(r): return [((x-X0)/self.step*s-.5,(y-Y0)/self.step*s-.5) for x,y in r]
  for p in polys:
   if not p['holes']: d.polygon(coords(p['outer']),fill=255)
   else:
    tmp=Image.new('L',(w,h));td=ImageDraw.Draw(tmp);td.polygon(coords(p['outer']),fill=255)
    for hole in p['holes']: td.polygon(coords(hole),fill=0)
    im=ImageChops.lighter(im,tmp);d=ImageDraw.Draw(im)
  return im.resize((self.w,self.h),Image.Resampling.BOX).point(lambda x:255 if x>=128 else 0)
 def holes(self):
  s=4; im=Image.new('L',(self.w*s,self.h*s));d=ImageDraw.Draw(im)
  for h in self.g['holes']:
   x,y=h['pos'];dx,dy=h['size'];rot=math.radians(h['rotation']);pts=[]
   # Actual board holes are round or oval. Polygon capsule for oval holes.
   if abs(dx-dy)<1e-8:
    for i in range(40):a=2*math.pi*i/40;pts.append((x+dx/2*math.cos(a),y+dy/2*math.sin(a)))
   else:
    # Ellipse approximation only for NPTH mechanical slots; not used in power vias.
    for i in range(40):
     a=2*math.pi*i/40;u=dx/2*math.cos(a);v=dy/2*math.sin(a);pts.append((x+u*math.cos(rot)+v*math.sin(rot),y-u*math.sin(rot)+v*math.cos(rot)))
   d.polygon([((a-X0)/self.step*s-.5,(b-Y0)/self.step*s-.5) for a,b in pts],fill=255)
  return im.resize((self.w,self.h),Image.Resampling.BOX).point(lambda x:255 if x>=128 else 0)
 def masks(self,net):
  holes=self.holes()
  return [ImageChops.subtract(self.draw(p),holes) for p in self.g['copper'][net]]
 def index(self,x,y):return int((y-Y0)/self.step)*self.w+int((x-X0)/self.step)
 def pos(self,i):return (X0+(i%self.w+.5)*self.step,Y0+(i//self.w+.5)*self.step)

def build(g,net,step,anchor,terms,dest):
 r=Raster(g,step);ms=r.masks(net);raw=[m.tobytes() for m in ms];N=r.n
 # Collect copper cell graph; merge cells underneath each device terminal to an ideal contact.
 ids=array('i',[-1])*(4*N);parents=[];positions=[]
 for l in range(4):
  for i,b in enumerate(raw[l]):
   if b:ids[l*N+i]=len(parents);parents.append(len(parents));positions.append(l*N+i)
 def root(a):
  while parents[a]!=a:parents[a]=parents[parents[a]];a=parents[a]
  return a
 def merge(a,b):parents[root(b)]=root(a)
 terminal_nodes=[]
 for ref,pin in [anchor]+terms:
  nodes=[]
  for p in g['pads']:
   if p['ref']==ref and p['pin']==pin:
    for l in p['layers']:
     # SMD power pads are all F.Cu. No cross-layer short is introduced for PTH here.
     if l!=0:continue
     pm=r.draw(p['polys'][str(l)]).tobytes()
     nodes.extend(ids[l*N+i] for i,b in enumerate(pm) if b and ids[l*N+i]>=0)
  if not nodes:raise RuntimeError(('terminal has no raster copper',ref,pin))
  for q in nodes[1:]:merge(nodes[0],q)
  terminal_nodes.append(nodes[0])
 edges=[]
 # Horizontal sheet conductance between adjacent equal-sized cells: t/rho.
 for l in range(4):
  gg=THICK[l]/RHO;off=l*N
  for i,b in enumerate(raw[l]):
   if not b:continue
   a=ids[off+i]
   if i%r.w+1<r.w and raw[l][i+1]:edges.append((a,ids[off+i+1],gg))
   if i+r.w<N and raw[l][i+r.w]:edges.append((a,ids[off+i+r.w],gg))
 # Vias/PTH barrels. Annulus is an ideal local contact; actual barrel R is retained.
 barrels=[v for v in g['vias'] if v['net']==net]
 barrels += [{'pos':p['pos'],'drill':min(p['drill']),'diameter':min(p['size']),'layers':p['layers']} for p in g['pads'] if p['net']==net and p['plated'] and min(p['drill'])>0]
 via_meta=[]
 for v in barrels:
  x,y=v['pos'];rad=v['diameter']/2;dr=v['drill'];ls=[]
  for l in v['layers']:
   ns=[]
   for iy in range(max(0,int((y-rad-Y0)/step)),min(r.h,int((y+rad-Y0)/step)+1)):
    for ix in range(max(0,int((x-rad-X0)/step)),min(r.w,int((x+rad-X0)/step)+1)):
     i=iy*r.w+ix;px,py=r.pos(i)
     if (px-x)**2+(py-y)**2<=(rad+step*.1)**2 and ids[l*N+i]>=0:ns.append(ids[l*N+i])
   if ns:
    for q in ns[1:]:merge(ns[0],q)
    ls.append((l,ns[0]))
  for (l,a),(ll,b) in zip(ls,ls[1:]):
   area=math.pi*(dr*1e-3+20e-6)*20e-6
   gg=area/(RHO*((Z[ll]-Z[l])*1e-3));edges.append((a,b,gg));via_meta.append((a,b,gg,x,y,l,ll))
 # Compact merged nodes and discard electrically floating islands.
 remap={};pre=[]
 for a in range(len(parents)):
  rt=root(a)
  if rt not in remap:remap[rt]=len(remap)
  pre.append(remap[rt])
 ed={}
 for a,b,gg in edges:
  a,b=pre[a],pre[b]
  if a==b:continue
  if a>b:a,b=b,a
  ed[a,b]=ed.get((a,b),0)+gg
 neighbors=[[] for _ in remap]
 for a,b in ed:neighbors[a].append(b);neighbors[b].append(a)
 terminal_nodes=[pre[q] for q in terminal_nodes];seen={terminal_nodes[0]};stack=list(seen)
 while stack:
  for b in neighbors[stack.pop()]:
   if b not in seen:seen.add(b);stack.append(b)
 if not all(t in seen for t in terminal_nodes):raise RuntimeError('Disconnected terminal at raster resolution')
 new={a:i for i,a in enumerate(sorted(seen))};node_map=array('i',[-1])*(4*N)
 for i,a in enumerate(ids):
  if a>=0:node_map[i]=new.get(pre[a],-1)
 final=[(new[a],new[b],gg) for (a,b),gg in ed.items() if a in seen]
 terminals=[new[t] for t in terminal_nodes]
 dest=Path(dest);dest.mkdir(parents=True,exist_ok=True)
 with open(dest/'graph.bin','wb') as f:
  f.write(struct.pack('4i',len(new),len(final),len(terms),terminals[0]))
  for a,b,gg in final:f.write(struct.pack('iid',a,b,gg))
  f.write(struct.pack(f'{len(terms)}i',*terminals[1:]))
 with open(dest/'map.bin','wb') as f:node_map.tofile(f)
 meta={'net':net,'step_mm':step,'w':r.w,'h':r.h,'nodes':len(new),'anchor':anchor,'terms':terms,'terminal_nodes':terminals,'x0':X0,'y0':Y0,'via_count':len(barrels),'floating_cells':sum(a>=0 and new.get(pre[a],-1)<0 for a in ids)}
 (dest/'meta.json').write_text(json.dumps(meta,indent=2))
 for l,m in enumerate(ms):m.save(dest/f'copper{l}.png')
 print(dest,meta,flush=True)
 result=subprocess.run([str(ROOT/'solve'),str(dest/'graph.bin'),str(dest/'volt.bin')],capture_output=True,text=True)
 (dest/'solver.log').write_text(result.stderr);print(result.stderr,flush=True);result.check_returncode()

if __name__=='__main__':
 geom=json.load(open(sys.argv[1]));step=float(sys.argv[2]);tag=sys.argv[3];which=sys.argv[4]
 specs={'batt':('+BATT',('R46','1'),[('Q1','1'),('Q3','1'),('Q5','1')]),'gnd':('GND',('J6','1'),[('R1','2'),('R3','2'),('R4','2')]),'input':('Net-(J1-Pin_1)',('J1','1'),[('R46','2')]),'spa':('/GATE DRIVE/SPA',('R1','1'),[('Q2','3')]),'spb':('/GATE DRIVE/SPB',('R3','1'),[('Q4','3')]),'spc':('/GATE DRIVE/SPC',('R4','1'),[('Q6','3')])}
 build(geom,*specs[which][:1],step,*specs[which][1:],ROOT/'work'/f'{tag}_{which}_{step}')
