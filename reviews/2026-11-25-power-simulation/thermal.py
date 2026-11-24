"""Coarse four-sheet thermal network, not CFD or junction-temperature qualification.
Includes actual copper coverage and vias, FR4 interlayer coupling, two exposed faces.
"""
from model import *
import argparse

def read_case(path,c):
 path=Path(path);m=json.loads((path/'meta.json').read_text());mp=array('i');mp.frombytes((path/'map.bin').read_bytes())
 v=array('d');v.frombytes((path/'volt.bin').read_bytes());n=m['nodes'];v=v[c*n:(c+1)*n]
 count=array('i',[0])*n;sx=array('d',[0])*n;sy=array('d',[0])*n;sl=array('d',[0])*n;wh=m['w']*m['h']
 for ri,ni in enumerate(mp):
  if ni<0:continue
  l,pi=divmod(ri,wh);iy,ix=divmod(pi,m['w']);count[ni]+=1;sx[ni]+=m['x0']+(ix+.5)*m['step_mm'];sy[ni]+=m['y0']+(iy+.5)*m['step_mm'];sl[ni]+=l
 for i in range(n):
  if count[i]:sx[i]/=count[i];sy[i]/=count[i];sl[i]/=count[i]
 powers=array('d',[0])*n
 data=(path/'graph.bin').read_bytes();_,ne,_,_=struct.unpack_from('4i',data)
 for a,b,g in struct.iter_unpack('iid',data[16:16+16*ne]):
  p=g*(v[a]-v[b])**2;powers[a]+=p*.5;powers[b]+=p*.5
 return [(sx[i],sy[i],round(sl[i]),p) for i,p in enumerate(powers) if p>1e-15]

def run(g,h,transient=False):
 step=.5;r=Raster(g,step);N=r.n;area=(step*1e-3)**2
 board=r.draw(g['outline']).tobytes();valid=[i for i,x in enumerate(board) if x];compact={i:j for j,i in enumerate(valid)};nn=len(valid);n=4*nn+1;ambient=n-1
 copper=[]
 for l in range(4):
  polys=[p for ps in g['copper'].values() for p in ps[l]]
  im=ImageChops.subtract(r.draw(polys),r.holes());copper.append(im.tobytes())
 # FR4 values are assumed, not identified production laminate data.
 kfr=.4;kcu=385;fr=[.05e-3,.67e-3,.67e-3,.05e-3];cvfr=1.8e6;cvcu=3.45e6
 conduct=[];capacity=array('d',[0])*n
 for l in range(4):
  conduct.append([kfr*fr[l]+kcu*THICK[l]*(copper[l][i]/255) for i in range(N)])
  for i in valid:capacity[l*nn+compact[i]]=area*(cvfr*fr[l]+cvcu*THICK[l]*(copper[l][i]/255))
 edges=[]
 for l in range(4):
  for i in valid:
   a=l*nn+compact[i]
   for j in [i+1 if i%r.w+1<r.w else -1,i+r.w]:
    if j not in compact:continue
    aa,bb=conduct[l][i],conduct[l][j];gg=2*aa*bb/(aa+bb);edges.append((a,l*nn+compact[j],gg))
   if l<3:edges.append((a,(l+1)*nn+compact[i],kfr*area/([.1e-3,1.24e-3,.1e-3][l])))
   if l in [0,3]:edges.append((a,ambient,h*area))
 barrels=g['vias']+[{'pos':p['pos'],'drill':min(p['drill']),'layers':p['layers']} for p in g['pads'] if p['plated'] and min(p['drill'])>0]
 for v in barrels:
  i=r.index(*v['pos'])
  if i not in compact:continue
  va=math.pi*(v['drill']*1e-3+20e-6)*20e-6
  for l,ll in zip(v['layers'],v['layers'][1:]):edges.append((l*nn+compact[i],ll*nn+compact[i],kcu*va/((Z[ll]-Z[l])*1e-3)))
 def deposit(points,scale=1):
  q=array('d',[0])*n
  for x,y,l,p in points:
   i=r.index(x,y)
   if i in compact:q[l*nn+compact[i]]+=p*scale
  return q
 def spot(q,ref,power):
  pads=[p for p in g['pads'] if p['ref']==ref and 0 in p['layers'] and (p['pin']=='1' or not ref.startswith('Q'))]
  ps=[poly for p in pads for poly in p['polys']['0']];mask=r.draw(ps).tobytes();ns=[compact[i] for i in valid if mask[i]]
  for i in ns:q[i]+=power/len(ns)
 # Copper at 80 C: resistance x1.236 versus the electrical model's 20 C.
 hot=1+.00393*60
 batt_points=read_case(ROOT/'work/new_batt_0.1',2)
 q15=deposit(batt_points,15**2*hot)
 qs=[];labels=[]
 for cur in [5,15,30]:qs.append(array('d',(q*(cur/15)**2 for q in q15)));labels.append(f'BATT_only_{cur}A_RMS_hotCu')
 # Deliberately specified sensitivity loads, not inferred actual switching loss.
 # 8 A DC input; 15 A RMS phases; average half-time low-side shunt conduction.
 # Copper: 15 A through far BATT/GND paths is a conservative layout heat source.
 qbase=array('d',q15)
 for key,c,current in [('gnd',2,15),('input',0,8),('spa',0,15/math.sqrt(2)),('spb',0,15/math.sqrt(2)),('spc',0,15/math.sqrt(2))]:
  sp=.2 if key=='gnd' else .1;points=read_case(ROOT/f'work/new_{key}_{sp}',c);qq=deposit(points,current**2*hot)
  for i in range(n):qbase[i]+=qq[i]
 for ref in ['R1','R3','R4']:spot(qbase,ref,.5*15**2*.002)
 spot(qbase,'R46',8**2*.001)
 # Electronics loss assumption, 1 W total (not a measured figure).
 for ref,p in [('U7',.35),('U4',.25),('U1',.2),('U2',.2)]:spot(qbase,ref,p)
 for fet_total in [2.73375,5.73375,8.73375]:
  q=array('d',qbase)
  for ref in ['Q1','Q2','Q3','Q4','Q5','Q6']:spot(q,ref,fet_total/6)
  qs.append(q);labels.append(f'phase15Arms_FETtotal{fet_total:.2f}W')
 if transient:
  # Constant 30 A BATT-copper-only heat, initially at ambient, 60 steps x 1 s.
  qs=[qs[2]]*60;labels=[f't={i+1}s' for i in range(60)]
  edges.extend((i,ambient,capacity[i]) for i in range(n-1))
 d=ROOT/'work'/f'thermal_h{h}_{"60s" if transient else "steady"}'
 d.mkdir(parents=True,exist_ok=True)
 with open(d/'graph.bin','wb') as f:
  f.write(struct.pack('4i',n,len(edges),len(qs),ambient))
  for a,b,gg in edges:f.write(struct.pack('iid',a,b,gg))
  f.write(struct.pack(f'{len(qs)}i',*([0]*len(qs))))
 with open(d/'rhs.bin','wb') as f:
  for q in qs:q.tofile(f)
 with open(d/'capacity.bin','wb') as f:capacity.tofile(f)
 args=[str(ROOT/'solve'),str(d/'graph.bin'),str(d/'temp.bin'),str(d/'rhs.bin')]
 if transient:args.append(str(d/'capacity.bin'))
 rr=subprocess.run(args,capture_output=True,text=True);(d/'solver.log').write_text(rr.stderr);rr.check_returncode()
 t=array('d');t.frombytes((d/'temp.bin').read_bytes());res=[]
 battmask=r.masks('+BATT')[3].tobytes();battnodes=[3*nn+compact[i] for i in valid if battmask[i]]
 for c,label in enumerate(labels):
  temp=t[c*n:(c+1)*n];source=qs[c];bal=sum(h*area*temp[l*nn+compact[i]] for i in valid for l in [0,3]);peak=max(temp)
  res.append({'case':label,'power_W':sum(source),'peak_rise_K':peak,'BATT_peak_rise_K':max(temp[i] for i in battnodes),'ambient_heat_flow_W':bal})
 meta={'step_mm':step,'h_each_face_W_m2K':h,'k_FR4_W_mK':kfr,'Cu_k_W_mK':kcu,'FR4_volumetric_heat_capacity_J_m3K':cvfr,'board_area_m2':nn*area,'nodes':n,'valid_cells':valid,'nn':nn,'w':r.w,'h':r.h,'results':res}
 (d/'meta.json').write_text(json.dumps(meta));print(json.dumps(res if not transient else res[-1:],indent=2),flush=True)
 return meta

if __name__=='__main__':run(json.load(open(sys.argv[1])),float(sys.argv[2]),len(sys.argv)>3)
