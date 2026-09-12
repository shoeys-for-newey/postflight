import numpy as np, json
from common import *
from curve import out_pct
from scipy import signal
# Cruiser-oriented checks: throttle usage / raw-stick estimate, current, sag, attitude, I-term behaviour,
# low-frequency band RMS, per-motor RPM efficiency, and the pooled step-response curves from stepfix.py.
AX=['roll','pitch','yaw']
for li in LOGS:
    c,t,fs,H=load(li)
    thr=c['rcCommand[3]']; idx=np.where(thr>1150)[0]; fl=slice(idx[0],idx[-1]); thp=(thr[fl]-1000)/10; tf=t[fl]
    mid,expo,hov=int(H['thr_mid']),int(H['thr_expo']),int(H['thr_hover'])
    print(f"== LOG {li} ({tf[-1]-tf[0]:.1f}s in flight): curve mid/expo/hover = {mid}/{expo}/{hov}")
    edges=np.arange(0,105,5); hist,_=np.histogram(thp,edges); pct=100*hist/hist.sum()
    print(" post-curve throttle usage %: "+' '.join(f"{edges[i]:.0f}-{edges[i+1]:.0f}:{pct[i]:.1f}" for i in range(len(hist)) if pct[i]>=0.5))
    print(f" post-curve thr median {np.median(thp):.1f} IQR {np.percentile(thp,25):.1f}-{np.percentile(thp,75):.1f} p95 {np.percentile(thp,95):.1f} max {thp.max():.1f}")
    xs=np.linspace(0,100,2001); ys=np.array([out_pct(mid,expo,hov,x) for x in xs])
    stick=np.interp(thp,ys,xs)
    print(f" est. RAW stick under that curve: median {np.median(stick):.1f}% IQR {np.percentile(stick,25):.1f}-{np.percentile(stick,75):.1f} p95 {np.percentile(stick,95):.1f} max {stick.max():.1f}")
    amps=c['amperageLatest'][fl]/100; vb=c['vbatLatest'][fl]/100
    print(f" amperage: min {amps.min():.1f} median {np.median(amps):.1f} p95 {np.percentile(amps,95):.1f} max {amps.max():.1f} A; mAh over {tf[-1]-tf[0]:.0f}s = {np.trapezoid(amps,tf)/3.6:.0f}; power median {np.median(amps*vb):.0f} W")
    tb=np.arange(0,65,5)
    print("   amps by thr bin: "+' '.join(f"{tb[i]:.0f}%:{np.median(amps[(thp>=tb[i])&(thp<tb[i+1])]):.1f}A" for i in range(len(tb)-1) if ((thp>=tb[i])&(thp<tb[i+1])).sum()>fs))
    n1=int(fs); print(f" vbat: first-second {vb[:n1].mean():.2f} last-second {vb[-n1:].mean():.2f} min {vb.min():.2f} V ({vb.min()/6:.2f} V/cell)")
    x=c['imuQuaternion[0]'][fl]/32767; y=c['imuQuaternion[1]'][fl]/32767; z=c['imuQuaternion[2]'][fl]/32767; w=np.sqrt(np.clip(1-x*x-y*y-z*z,0,1))
    roll=np.degrees(np.arctan2(2*(w*x+y*z),1-2*(x*x+y*y))); pitch=np.degrees(np.arcsin(np.clip(2*(w*y-z*x),-1,1)))
    print(f" attitude (BF sign): pitch median {np.median(pitch):+.1f} deg (IQR {np.percentile(pitch,25):+.1f}..{np.percentile(pitch,75):+.1f}, extremes {pitch.min():+.1f}/{pitch.max():+.1f}); roll median {np.median(roll):+.1f}; |roll|>30: {100*np.mean(np.abs(roll)>30):.1f}%  |pitch|>30: {100*np.mean(np.abs(pitch)>30):.1f}%")
    ph=np.abs(pitch)
    for lo_,hi_ in ((0,5),(5,10),(10,20),(20,30),(30,90)):
        m=(ph>=lo_)&(ph<hi_)
        if m.sum()>fs: print(f"   |pitch| {lo_:2d}-{hi_:2d} deg: {100*m.mean():4.1f}% of flight, thr median {np.median(thp[m]):.1f}%, amps median {np.median(amps[m]):.1f}")
    I=[c['axisI[%d]'%i][fl] for i in range(3)]; sp=[c['setpoint[%d]'%i][fl] for i in range(3)]; gy=[c['gyroADC[%d]'%i][fl] for i in range(3)]
    for i,a in enumerate(AX):
        r=np.corrcoef(I[i],sp[i])[0,1]; r2=np.corrcoef(I[i],thp)[0,1]; r3=np.corrcoef(I[i],pitch)[0,1]
        big=np.abs(I[i])>100
        print(f" I-term {a:5s}: mean {I[i].mean():+6.1f} rms {np.sqrt(np.mean(I[i]**2)):5.1f} corr(sp) {r:+.2f} corr(thr) {r2:+.2f} corr(pitch) {r3:+.2f}; |I|>100: {100*big.mean():.1f}%"
              + (f" (sp median {np.median(sp[i][big]):+.0f} deg/s, err median {np.median((sp[i]-gy[i])[big]):+.1f})" if big.any() else ""))
    for lo_,hi_ in ((3,12),(12,30),(30,70),(70,110)):
        sos=signal.butter(4,[lo_/(fs/2),hi_/(fs/2)],'band',output='sos')   # sos: ba form overflows for 3-12 Hz at 2 kHz
        s=' '.join(f"{ax}:{np.sqrt(np.mean(signal.sosfiltfilt(sos,gy[i])**2)):.1f}/{np.sqrt(np.mean(signal.sosfiltfilt(sos,sp[i])**2)):.1f}" for i,ax in enumerate(['R','P','Y']))
        print(f" band RMS {lo_:3d}-{hi_:3d}Hz gyro/setpoint (deg/s): {s}")
    erpm=np.stack([c[f'eRPM[{i}]'][fl] for i in range(4)],1)*100/pole_pairs(H); mot=np.stack([c[f'motor[{i}]'][fl] for i in range(4)],1); motp=(mot-48)/1999*100
    calm=np.abs(np.stack(gy,1)).max(1)<60
    mr=np.median(erpm[calm],0); mp=np.median(motp[calm],0)
    print(" per-motor (calm) RR/FR/RL/FL: RPM "+' '.join(f"{v:.0f}" for v in mr)+" | motor% "+' '.join(f"{v:.1f}" for v in mp)+" | RPM per motor%: "+' '.join(f"{a/b:.0f}" for a,b in zip(mr,mp)))
    # RPM vs motor% slope per motor (linear fit over calm samples) -> motor/prop health
    for k in range(4):
        A=np.polyfit(motp[calm,k],erpm[calm,k],1); print(f"   motor{k+1} RPM = {A[0]:.0f}*motor% + {A[1]:.0f}")
S=json.load(open('step_results.json'))
pts=[0,5,10,15,20,25,30,40,50,60,80,100,120,150,200,250]
for a in AX:
    for k in ('lo','hi'):
        if k in S['pooled'][a]:
            cur=S['pooled'][a][k]['curve']; n=len(cur); fs_=n/0.25
            print(f" step {a:5s} {k}: "+' '.join(f"{ms}:{cur[int(ms/1000*fs_)]:.2f}" for ms in pts if int(ms/1000*fs_)<n))
