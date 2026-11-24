from model import *
from thermal import read_case

def electrical(d):
 m=json.loads((d/'meta.json').read_text());mp=array('i');mp.frombytes((d/'map.bin').read_bytes());v=array('d');v.frombytes((d/'volt.bin').read_bytes())
 return m,mp,v

def cuts(d,case,ys):
 m,mp,v=electrical(d);w,h,s,n=m['w'],m['h'],m['step_mm'],m['nodes'];N=w*h;result=[]
 for y in ys:
  iy=round((y-Y0)/s-.5);lr=[]
  for l in range(4):
   runs=[];active=[]
   for ix in range(w):
    i=iy*w+ix;a=mp[l*N+i];b=mp[l*N+i+w]
    if a>=0 and b>=0:active.append((ix,(v[case*n+a]-v[case*n+b])*THICK[l]/RHO))
    elif active:runs.append(active);active=[]
   if active:runs.append(active)
   lr.append([{'x_start_mm':X0+a[0][0]*s,'x_end_mm':X0+(a[-1][0]+1)*s,'fraction_signed':sum(z[1] for z in a)} for a in runs])
  result.append({'y_mm':y,'layer_intervals':lr,'layer_totals':[sum(a['fraction_signed'] for a in lr[l]) for l in range(4)]})
 return result

def font(size):return ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf',size)
def color(v):
 # Perceptually ordered dark purple -> orange -> pale yellow, shared scale.
 stops=[(30,24,70),(86,43,110),(170,60,90),(231,120,55),(252,222,120)]
 a=min(1,max(0,v))*(len(stops)-1);j=min(len(stops)-2,int(a));t=a-j
 return tuple(round(stops[j][k]*(1-t)+stops[j+1][k]*t) for k in range(3))

def current_image(d,case,layer,view,scale_max):
 m,mp,v=electrical(d);w,h,s,n=m['w'],m['h'],m['step_mm'],m['nodes'];N=w*h;im=Image.new('RGB',(w,h),(224,230,235));pix=im.load()
 for i in range(N):
  a=mp[layer*N+i]
  if a<0:continue
  x=i%w;y=i//w;ixs=[];iys=[]
  for off,target in [(-1,ixs),(1,ixs),(-w,iys),(w,iys)]:
   j=i+off
   if j<0 or j>=N or (abs(off)==1 and j//w!=y):continue
   b=mp[layer*N+j]
   if b>=0:target.append(abs(v[case*n+a]-v[case*n+b])*THICK[layer]/RHO*15/s)
  jx=sum(ixs)/max(1,len(ixs));jy=sum(iys)/max(1,len(iys));mag=math.hypot(jx,jy);pix[x,y]=color(mag/scale_max)
 x0,y0,x1,y1=view
 return im.crop((round((x0-X0)/s),round((y0-Y0)/s),round((x1-X0)/s),round((y1-Y0)/s))).resize((350,630))

def plot_currents():
 view=(130,77,152,131);im=Image.new('RGB',(790,760),'white');d=ImageDraw.Draw(im)
 d.text((20,10),'BATT copper, 15 A to phase W (B.Cu)',font=font(23),fill='#17252f')
 for k,tag in enumerate(['old','new']):
  panel=current_image(ROOT/'work'/f'{tag}_batt_0.1',2,3,view,3)
  im.paste(panel,(25+k*390,75));d.text((25+k*390,48),'Nov 24' if k==0 else 'Nov 25',font=font(18),fill='black')
  yy=75+(116-view[1])/(view[3]-view[1])*630
  d.line((25+k*390,yy,375+k*390,yy),fill='#00bad0',width=2)
  if k==1:d.text((430,yy+5),'56%          44%',font=font(17),fill='white',stroke_width=2,stroke_fill='black')
 for x in range(250):d.line((25+x,721,25+x,735),fill=color(x/249))
 d.text((290,716),'0 to 3 A/mm sheet current; >3 clipped',font=font(15),fill='black')
 d.text((25,740),'Cross-section: y = 116 mm. DC model; not a switching/EMI simulation.',font=font(13),fill='black')
 im.save(ROOT/'batt_current_comparison.png')

def plot_thermal():
 p=ROOT/'work/thermal_h10.0_steady';m=json.loads((p/'meta.json').read_text());t=array('d');t.frombytes((p/'temp.bin').read_bytes());n=m['nodes'];nn=m['nn'];w=m['w'];h=m['h'];valid=m['valid_cells']
 cases=[(1,20,'BATT only: 15 A RMS'),(4,120,'Specified total heat: 8.32 W')]
 im=Image.new('RGB',(780,650),'white');d=ImageDraw.Draw(im);d.text((15,10),'Approximate board temperature rise (F.Cu)',font=font(22),fill='black')
 for k,(c,top,title) in enumerate(cases):
  panel=Image.new('RGB',(w,h),(225,230,235));pix=panel.load()
  for j,i in enumerate(valid):pix[i%w,i//w]=color(t[c*n+j]/top)
  panel=panel.resize((350,450));im.paste(panel,(20+390*k,85));d.text((20+390*k,52),title,font=font(17),fill='black')
  for x in range(250):d.line((20+390*k+x,550,20+390*k+x,567),fill=color(x/249))
  d.text((20+390*k,575),f'0 to {top} K rise; peak {m["results"][c]["peak_rise_K"]:.1f} K',font=font(16),fill='black')
 d.text((20,610),'h = 10 W/m2/K per face; no cable/chassis cooling. Sensitivity model.',font=font(15),fill='black')
 d.text((20,632),'Right panel uses assumed device losses, not measured operating losses.',font=font(14),fill='black')
 im.save(ROOT/'thermal_sensitivity.png')

if __name__=='__main__':
 results={'electrical':{},'thermal':{}}
 for p in sorted((ROOT/'work').iterdir()):
  if (p/'volt.bin').exists():
   m,mp,v=electrical(p);n=m['nodes']
   results['electrical'][p.name]={'terminals':m['terms'],'R_ohm':[v[i*n+t] for i,t in enumerate(m['terminal_nodes'][1:])],'nodes':n,'step_mm':m['step_mm']}
  if (p/'temp.bin').exists():
   m=json.loads((p/'meta.json').read_text());results['thermal'][p.name]={k:v for k,v in m.items() if k not in ['valid_cells']}
 results['BATT_W_cuts']=cuts(ROOT/'work/new_batt_0.1',2,[100,108,116,118,125])
 results['GND_W_cuts']=cuts(ROOT/'work/new_gnd_0.2',2,[90,105,120,130])
 (ROOT/'results.json').write_text(json.dumps(results,indent=2))
 plot_currents();plot_thermal()
 print('W return layer fractions:',[(a['y_mm'],a['layer_totals']) for a in results['GND_W_cuts']])
