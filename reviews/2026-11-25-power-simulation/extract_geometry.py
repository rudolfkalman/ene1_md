"""Run using KiCad's Python. Read only; export actual filled copper and pad terminals."""
import pcbnew as k
import json,sys

b=k.LoadBoard(sys.argv[1]); layers=[k.F_Cu,k.In1_Cu,k.In2_Cu,k.B_Cu]
def xy(v): return [v.x/1e6,v.y/1e6]
def poly(p):
    def ring(c): return [xy(c.CPoint(i)) for i in range(c.PointCount())]
    return [{'outer':ring(p.Outline(i)), 'holes':[ring(p.Hole(i,j)) for j in range(p.HoleCount(i))]} for i in range(p.OutlineCount())]
def shape(obj,layer):
    p=k.SHAPE_POLY_SET(); obj.TransformShapeToPolygon(p,layer,0,3000,k.ERROR_INSIDE)
    return poly(p)
out={'layers':[b.GetLayerName(l) for l in layers],'copper':{},'pads':[],'vias':[],'holes':[],'zones':[]}
def add(net,layer,p):
    if p: out['copper'].setdefault(net,[[],[],[],[]])[layers.index(layer)].extend(p)
for z in b.Zones():
    for l in layers:
        if z.GetLayerSet().Contains(l):
            p=poly(z.GetFilledPolysList(l)); add(z.GetNetname(),l,p)
            out['zones'].append({'net':z.GetNetname(),'layer':layers.index(l),'polys':p})
for f in b.GetFootprints():
    for p in f.Pads():
        ls=[l for l in layers if p.GetLayerSet().Contains(l)]
        if not ls: continue
        polys={layers.index(l):shape(p,l) for l in ls}
        for l in ls: add(p.GetNetname(),l,polys[layers.index(l)])
        out['pads'].append({'ref':f.GetReference(),'pin':p.GetNumber(),'net':p.GetNetname(),'pos':xy(p.GetPosition()),'size':xy(p.GetSize()),'rotation':p.GetOrientationDegrees(),'layers':[layers.index(l) for l in ls],'polys':polys,'drill':xy(p.GetDrillSize()),'plated':p.GetAttribute()==k.PAD_ATTRIB_PTH})
        if p.GetDrillSize().x:
            out['holes'].append({'pos':xy(p.GetPosition()),'size':xy(p.GetDrillSize()),'rotation':p.GetOrientationDegrees()})
for t in b.GetTracks():
    ls=[l for l in layers if t.GetLayerSet().Contains(l)]
    for l in ls: add(t.GetNetname(),l,shape(t,l))
    if t.GetClass()=='PCB_VIA':
        out['vias'].append({'pos':xy(t.GetPosition()),'net':t.GetNetname(),'layers':[layers.index(l) for l in ls],'drill':t.GetDrillValue()/1e6,'diameter':t.GetWidth(k.F_Cu)/1e6})
        out['holes'].append({'pos':xy(t.GetPosition()),'size':[t.GetDrillValue()/1e6]*2,'rotation':0})
p=k.SHAPE_POLY_SET(); b.GetBoardPolygonOutlines(p,False,None,False,True); out['outline']=poly(p)
with open(sys.argv[2],'w') as f: json.dump(out,f,separators=(',',':'))
print(sys.argv[2],len(out['pads']),'pads',len(out['vias']),'vias',len(out['copper']),'nets')
