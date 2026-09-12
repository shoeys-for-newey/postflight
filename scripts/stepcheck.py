import numpy as np, json
from common import *
from scipy import signal
# Independent cross-check of stepfix.py: time-domain least-squares FIR (setpoint -> gyro) at ~1 kHz, cumsum -> step.
# Also re-derives hover throttle restricted to near-level attitude, using body-Z and earth-vertical accel.
AX=['roll','pitch','yaw']
OUT={}
for li in LOGS:
    OUT[li]={'step':{},'hover':{}}
    c,t,fs,H=load(li); thr=c['rcCommand[3]']; idx=np.where(thr>1150)[0]; fl=slice(idx[0],idx[-1])
    dec=max(1,int(round(fs/1000))); fs1=fs/dec
    L=int(0.3*fs1)  # 300 ms FIR
    print(f"== LOG {li}: LS-FIR step check (fs {fs1:.0f} Hz, {L} taps)")
    for i,a in enumerate(AX):
        sp=signal.decimate(c['setpoint[%d]'%i][fl],dec,zero_phase=True); gy=signal.decimate(c['gyroADC[%d]'%i][fl],dec,zero_phase=True)
        w=int(fs1); keep=np.zeros(len(sp),bool)
        for s in range(0,len(sp)-w,w//2):
            if np.abs(sp[s:s+w]).max()>40: keep[s:s+w]=True
        RtR=np.zeros((L,L)); Rty=np.zeros(L)
        idxk=np.where(keep)[0]; idxk=idxk[idxk>=L]
        if len(idxk) < 2*L:   # axis never excited (calm hover): solve would be singular
            print(f"  {a:5s} n={len(idxk)}: insufficient excitation - skipped"); continue
        for s in range(0,len(idxk),20000):
            ii=idxk[s:s+20000]
            X=np.stack([sp[ii-k] for k in range(L)],1)
            RtR+=X.T@X; Rty+=X.T@gy[ii]
        lam=1e-3*np.trace(RtR)/L
        h=np.linalg.solve(RtR+lam*np.eye(L),Rty)
        st=np.cumsum(h); tt=np.arange(L)/fs1*1000
        lim=int(0.2*fs1); pk=st[:lim].max(); tpk=tt[st[:lim].argmax()]
        def tc(v):
            j=np.where(st>=v)[0]; return f"{tt[j[0]]:.1f}" if len(j) else "n/a"
        ss=st[int(0.2*fs1):].mean()
        OUT[li]['step'][a]={'n':int(len(idxk)),'t50':tc(0.5),'t90':tc(0.9),'t100':tc(1.0),'peak':round(float(pk),3),'t_peak':round(float(tpk),1),'overshoot_pct':round(float(100*(pk-1)),1),'ss':round(float(ss),3),'fs':fs1,'curve':[round(float(v),4) for v in st]}
        print(f"  {a:5s} n={len(idxk)}: t50 {tc(0.5)} t90 {tc(0.9)} t100 {tc(1.0)} ms  peak {pk:.2f} @{tpk:.0f}ms  overshoot {100*(pk-1):+.0f}%  ss(200-300ms) {ss:.2f}  | curve: "+' '.join(f"{ms}:{st[int(ms/1000*fs1)]:.2f}" for ms in (10,20,30,50,80,120,200,280)))
    r22,aupA,aupB=quat_R22_and_aup(c); tilt=np.degrees(np.arccos(np.clip(r22,-1,1)))
    thp=(thr-1000)/10; az=c['accSmooth[2]']/2048.0
    g=np.stack([c['gyroADC[%d]'%k] for k in range(3)],1); gmax=np.abs(g).max(1)
    thr_l=lp(thr,fs,4); az_l=lp(az,fs,4); g_l=lp(gmax,fs,4); aup_l=lp(aupA,fs,4)
    m=(g_l<120)&(thp>8)&(thp<75); m[:idx[0]]=False; m[idx[-1]:]=False
    for name,mm in (('all',m),('tilt<5deg',m&(tilt<5)),('tilt<10deg',m&(tilt<10)),('tilt 10-25',m&(tilt>=10)&(tilt<25))):
        x=(thr_l[mm]-1000)/10; y=az_l[mm]; yu=aup_l[mm]
        bins=np.arange(10,60,2.5); ib=np.digitize(x,bins)
        med=[(bins[j-1]+1.25,np.median(y[ib==j]),np.median(yu[ib==j]),(ib==j).sum()) for j in range(1,len(bins)) if (ib==j).sum()>400]
        def cross(col):
            for j in range(len(med)-1):
                b0,b1=med[j],med[j+1]
                if (b0[col]-1)*(b1[col]-1)<=0 and b1[col]!=b0[col]: return b0[0]+(1-b0[col])*(b1[0]-b0[0])/(b1[col]-b0[col])
            return float('nan')
        OUT[li]['hover'][name]={'n':int(mm.sum()),'cross_body':float(cross(1)),'cross_earth':float(cross(2)),'bins':[(float(b),float(v),float(u),int(n)) for b,v,u,n in med]}
        print(f"  hover ({name:10s} n={mm.sum():6d}): accZ(body)=1G at {cross(1):.1f}%  earth-vertical=1G at {cross(2):.1f}%  | bins "+' '.join(f"{b:.1f}:{v:.2f}/{u:.2f}" for b,v,u,n in med))
json.dump(OUT,open('stepcheck_results.json','w'))
