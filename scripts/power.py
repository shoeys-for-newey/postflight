import numpy as np, json, sys
from common import *
# Power/thrust comparison metrics for one run folder, pooled over its logs.  usage: python ../../power.py <motor KV> <cells>
KV=float(sys.argv[1]); CELLS=int(sys.argv[2]) if len(sys.argv)>2 else 4
acc=dict(thr=[],az=[],aup=[],motp=[],rpm=[],vb=[],tilt=[],calm=[],amps=[])
for li in LOGS:
    c,t,fs,H=load(li); thr=c['rcCommand[3]']; idx=np.where(thr>1150)[0]; fl=slice(idx[0],idx[-1])
    r22,aupA,_=quat_R22_and_aup(c); tilt=np.degrees(np.arccos(np.clip(r22,-1,1)))
    g=np.stack([c[f'gyroADC[{i}]'] for i in range(3)],1); gmax=np.abs(g).max(1)
    mot=np.stack([c[f'motor[{i}]'] for i in range(4)],1); motp=((mot-48)/1999*100).mean(1)
    rpm=np.stack([c[f'eRPM[{i}]'] for i in range(4)],1).mean(1)*100/pole_pairs(H)
    az=lp(c['accSmooth[2]']/2048.0,fs,4); aup=lp(aupA,fs,4); thl=lp(thr,fs,4); gl=lp(gmax,fs,4); ml=lp(motp,fs,4)
    for k,v in (('thr',(thl-1000)/10),('az',az),('aup',aup),('motp',ml),('rpm',lp(rpm,fs,4)),('vb',c['vbatLatest']/100),('tilt',tilt),('calm',gl<120),('amps',c['amperageLatest']/100)):
        acc[k].append(v[fl])
A={k:np.concatenate(v) for k,v in acc.items()}
thr,az,aup,motp,rpm,vb,tilt,calm,amps=[A[k] for k in ('thr','az','aup','motp','rpm','vb','tilt','calm','amps')]
print(f"logs {LOGS}, {len(thr)/1000:.0f} k samples in flight, vbat {vb.min():.2f}-{vb.max():.2f} V, KV {KV:.0f}, no-load RPM at median vbat {KV*np.median(vb):.0f}")
hov=calm&(np.abs(az-1)<0.05)&(tilt<10)&(thr>8)
print(f"HOVER (|az-1|<0.05, tilt<10, calm; n={hov.sum()}): thr {np.median(thr[hov]):.1f}%  motor {np.median(motp[hov]):.1f}%  RPM {np.median(rpm[hov]):.0f}  vbat {np.median(vb[hov]):.2f}  amps {np.median(amps[hov]):.1f}")
m=calm&(thr>8)&(thr<75)
Af=np.polyfit(thr[m],az[m],1); Am=np.polyfit(motp[m],az[m],1)
print(f"accel-Z slope: {Af[0]:.4f} G per throttle-%   {Am[0]:.4f} G per motor-%   (zero-G intercepts thr {(-Af[1]/Af[0]):.1f}%, motor {(-Am[1]/Am[0]):.1f}%)")
# RPM vs motor% (calm): binned medians + fit
print(" RPM by motor-% bin:", ' '.join(f"{b}:{np.median(rpm[calm&(motp>=b-2.5)&(motp<b+2.5)]):.0f}" for b in range(20,100,10) if (calm&(motp>=b-2.5)&(motp<b+2.5)).sum()>200))
Ar=np.polyfit(motp[calm&(motp>15)&(motp<70)],rpm[calm&(motp>15)&(motp<70)],1)
print(f" RPM = {Ar[0]:.0f} * motor% + {Ar[1]:.0f}  (fit 15-70 %) -> extrapolated 100 %: {Ar[0]*100+Ar[1]:.0f} = {100*(Ar[0]*100+Ar[1])/(KV*np.median(vb)):.0f}% of no-load")
for lo in (80,90,95):
    mm=motp>lo
    if mm.sum()>50: print(f" motor>{lo}% (n={mm.sum()}, {mm.sum()/4000:.1f}s@4k): RPM median {np.median(rpm[mm]):.0f} max {rpm[mm].max():.0f}  vbat {np.median(vb[mm]):.2f} (min {vb[mm].min():.2f})  earth-up accel median {np.median(aup[mm]):.2f} G p90 {np.percentile(aup[mm],90):.2f}  body-Z p90 {np.percentile(az[mm],90):.2f}  amps median {np.median(amps[mm]):.1f}")
print(f" RPM 99.9pct {np.percentile(rpm,99.9):.0f}, max {rpm.max():.0f}; per-motor max {[float(v) for v in [rpm.max()]]}")
hr=np.median(rpm[hov]); top=np.percentile(rpm[motp>90],50) if (motp>90).sum()>50 else np.percentile(rpm,99.5)
print(f" thrust ratio (RPM_top/RPM_hover)^2 = ({top:.0f}/{hr:.0f})^2 = {(top/hr)**2:.2f}  (static, same prop => thrust ~ RPM^2)")
# sag: vbat drop from calm-hover to motor>90 within the same pack state
print(f" sag: vbat at hover {np.median(vb[hov]):.2f} vs motor>90 {np.median(vb[motp>90]) if (motp>90).any() else float('nan'):.2f}")
