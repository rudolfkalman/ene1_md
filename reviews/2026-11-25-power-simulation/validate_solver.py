"""Independent analytical checks of the numerical solver, not PCB certification."""
from model import *

def solve(name,n,edges,anchor,terminal,rhs=None):
 d=ROOT/'work'/name;d.mkdir(exist_ok=True,parents=True)
 with open(d/'graph.bin','wb') as f:
  f.write(struct.pack('4i',n,len(edges),1,anchor))
  for e in edges:f.write(struct.pack('iid',*e))
  f.write(struct.pack('i',terminal))
 args=[str(ROOT/'solve'),str(d/'graph.bin'),str(d/'out.bin')]
 if rhs is not None:
  with open(d/'rhs.bin','wb') as f:array('d',rhs).tofile(f)
  args.append(str(d/'rhs.bin'))
 rr=subprocess.run(args,capture_output=True,text=True);rr.check_returncode()
 v=array('d');v.frombytes((d/'out.bin').read_bytes());return v,rr.stderr

checks=[]
v,log=solve('check_parallel',2,[(0,1,1000),(0,1,500)],0,1)
expected=1/1500;checks.append({'check':'1 mOhm || 2 mOhm','actual_ohm':v[1],'expected_ohm':expected,'relative_error':abs(v[1]/expected-1)})
# 40 x 10 mm sheet, 35 um Cu, 0.1 mm cells; ideal end contacts.
nx=401;ny=100;ids={};nextid=2
for x in range(nx):
 for y in range(ny):
  if x==0:ids[x,y]=0
  elif x==nx-1:ids[x,y]=1
  else:ids[x,y]=nextid;nextid+=1
edges=[];gg=35e-6/RHO
for (x,y),a in ids.items():
 for xy in [(x+1,y),(x,y+1)]:
  if xy in ids and ids[xy]!=a:edges.append((a,ids[xy],gg))
v,log=solve('check_strip',nextid,edges,0,1);expected=RHO*.04/(.01*35e-6)
checks.append({'check':'40 x 10 mm, 35 um sheet','actual_ohm':v[1],'expected_ohm':expected,'relative_error':abs(v[1]/expected-1)})
# Uniform 100 x 100 mm thermal plate, 1 W evenly deposited, h=10 both sides.
nx=20;N=nx*nx;edges=[];area=.01/N;rhs=[1/N]*N+[0]
for i in range(N):
 edges.append((i,N,2*10*area))
 for j in [i+1 if i%nx<nx-1 else -1,i+nx]:
  if 0<=j<N:edges.append((i,j,.02))
v,log=solve('check_thermal',N+1,edges,N,0,rhs);expected=1/(2*10*.01)
checks.append({'check':'uniform thermal plate 1 W','actual_max_K':max(v),'expected_K':expected,'relative_error':abs(max(v)/expected-1)})
assert all(c['relative_error']<1e-6 for c in checks)
(ROOT/'solver_validation.json').write_text(json.dumps(checks,indent=2));print(json.dumps(checks,indent=2))
