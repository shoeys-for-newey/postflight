import numpy as np, json, sys
# Betaflight 2025.12 quadratic-Bezier throttle curve replica (fc/rc.c initRcProcessing).
# usage: python curve.py <thr_mid> <thr_expo> <thr_hover> <measured_hover_out_pct> [alt1_mid,expo,hover ...]
def qb(x,p0x,p1x,p2x,p0y,p1y,p2y):
    a=p0x-2*p1x+p2x; b=2*(p1x-p0x); c=p0x-x
    if abs(a)<1e-9: t=-c/b
    else:
        d=max(b*b-4*a*c,0); t=(-b+np.sqrt(d))/(2*a)
        if t<0 or t>1: t=(-b-np.sqrt(d))/(2*a)
    t=min(max(t,0),1)
    return (1-t)**2*p0y+2*(1-t)*t*p1y+t*t*p2y
def table(mid,expo,hover):
    mid/=100; expo/=100; hover/=100
    cp1x=mid*0.5; cp1y=hover*0.5*(1+expo); cp2x=(1+mid)*0.5; cp2y=1+((hover-1)*0.5*(1+expo))
    L=[]
    for i in range(12):
        x=i/11
        y=qb(x,0,cp1x,mid,0,cp1y,hover) if x<=mid else qb(x,mid,cp2x,1,hover,cp2y,1)
        L.append(int(round(1000+1000*y)))
    return L
def lookup(L,tmp):  # tmp 0..1000
    sc=tmp*11; idx=int(sc//1000); rem=sc%1000
    if idx>=11: return L[11]
    return L[idx]+(L[idx+1]-L[idx])*rem/1000
def out_pct(mid,expo,hover,stick_pct,mincheck=1050):
    L=table(mid,expo,hover)
    rc=1000+10*stick_pct
    tmp=max(0,(min(rc,2000)-mincheck))*1000/(2000-mincheck)
    return (lookup(L,tmp)-1000)/10
def stick_for_output(mid,expo,hover,target):
    xs=np.linspace(0,100,10001); ys=np.array([out_pct(mid,expo,hover,x) for x in xs])
    return float(xs[np.argmin(np.abs(ys-target))])
def slope(mid,expo,hover,s):
    return (out_pct(mid,expo,hover,s+1)-out_pct(mid,expo,hover,s-1))/2
if __name__=='__main__':
    a=sys.argv[1:]
    cur=tuple(int(v) for v in a[:3]) if len(a)>=3 else (30,40,27)
    tgt=float(a[3]) if len(a)>=4 else 34.5
    props=[tuple(int(v) for v in s.split(',')) for s in a[4:]] or [(30,40,34),(30,40,35),(35,40,35),(37,40,35),(35,50,35)]
    print("current mid/expo/hover=%s table:"%(cur,),table(*cur))
    print("current curve: stick->out", [(s,round(out_pct(*cur,s),1)) for s in range(0,101,10)])
    sc=stick_for_output(*cur,tgt); print("stick position giving %.1f%% output under CURRENT curve: %.1f%% (raw rc %d)"%(tgt,sc,1000+10*sc))
    print("current sensitivity at hover stick (%.1f%%): %.2f out%%/stick%%; at stick 30%%: %.2f; at 50%%: %.2f"%(sc,slope(*cur,sc),slope(*cur,30),slope(*cur,50)))
    J={'cur':[(s,out_pct(*cur,s)) for s in range(0,101)]}
    for k,prop in enumerate(props):
        s=stick_for_output(*prop,tgt)
        print("proposed mid/expo/hover=%s: hover stick %.1f%%, sensitivity at hover %.2f, at 50%% stick %.2f, out@25/50/75 stick = %.1f/%.1f/%.1f, table %s"%(prop,s,slope(*prop,s),slope(*prop,50),out_pct(*prop,25),out_pct(*prop,50),out_pct(*prop,75),table(*prop)))
        J['prop' if k==0 else 'prop%d'%(k+1)]=[(s,out_pct(*prop,s)) for s in range(0,101)]
    json.dump(J,open('curve.json','w'))
