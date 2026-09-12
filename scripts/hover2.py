import numpy as np, json
from common import *
out={}
allx=[]; ally=[]; allv=[]
for li in LOGS:
    c, t, fs, H = load(li)
    thr = c['rcCommand[3]']; thp=(thr-1000)/10
    g = np.stack([c['gyroADC[0]'], c['gyroADC[1]'], c['gyroADC[2]']],1); gmax=np.abs(g).max(1)
    az = c['accSmooth[2]']/2048.0
    vb = c['vbatLatest']/100
    mot = np.stack([c[f'motor[{i}]'] for i in range(4)],1); motp=(mot-48)/1999*100
    erpm = np.stack([c[f'eRPM[{i}]'] for i in range(4)],1)
    inflight = np.zeros(len(t),bool)
    idx=np.where(thr>1150)[0]
    if len(idx)==0:
        print(f"LOG {li}: throttle never above 1150 (bench/idle log) - skipped"); continue
    inflight[idx[0]:idx[-1]]=True
    thr_l=lp(thr,fs,4); az_l=lp(az,fs,4); g_l=lp(gmax,fs,4)
    # rotation-induced accel at IMU offset -> exclude high rates; also exclude rapid throttle changes (thrust lag)
    dthr=np.abs(np.gradient(thr_l))*fs  # %/s in raw units
    m = inflight & (g_l<120) & (dthr<800) & (thp>8) & (thp<75)
    x=(thr_l[m]-1000)/10; y=az_l[m]
    allx.append(x); ally.append(y); allv.append(vb[m])
    bins=np.arange(10,72,2.5); idx=np.digitize(x,bins)
    med=[(bins[i-1]+1.25, np.median(y[idx==i]), (idx==i).sum()) for i in range(1,len(bins)) if (idx==i).sum()>400]
    bx=np.array([b for b,v,n in med]); by=np.array([v for b,v,n in med])
    # zero crossing of by-1 via interpolation
    cross=None
    for i in range(len(bx)-1):
        if (by[i]-1)*(by[i+1]-1)<=0 and by[i+1]!=by[i]:
            cross=bx[i]+(1-by[i])*(bx[i+1]-bx[i])/(by[i+1]-by[i]); break
    A=np.polyfit(x,y,1)
    print(f"LOG {li}: n={m.sum()} ({m.sum()/fs:.1f}s)  accZ=1G at thr% (binned) = {(cross if cross is not None else float('nan')):.1f}  | linear fit zero = {(1-A[1])/A[0]:.1f}  slope {A[0]:.4f} G/%")
    print("   bins: "+' '.join(f"{b:.1f}:{v:.2f}" for b,v,n in med))
    # hover-ish samples: accZ within 0.95-1.05 and low rates -> throttle distribution
    hm = m & (np.abs(az_l-1)<0.05)
    hp=(thr_l[hm]-1000)/10
    print(f"   near-1G samples: n={hm.sum()} thr% median {np.median(hp):.1f} IQR {np.percentile(hp,25):.1f}-{np.percentile(hp,75):.1f}; motor% median {np.median(motp[hm].mean(1)):.1f}; RPM median {np.median(erpm[hm].mean(1))*100/pole_pairs(H):.0f}; vbat median {np.median(vb[hm]):.2f}")
    # split by vbat halves
    for lo,hi in ((vb[hm].min(), np.median(vb[hm])),(np.median(vb[hm]), vb[hm].max())):
        mm=hm&(vb>=lo)&(vb<=hi)
        if mm.sum()>500: print(f"     vbat {lo:.2f}-{hi:.2f}: thr% median {np.median((thr_l[mm]-1000)/10):.1f} (n={mm.sum()})")
    out[li]={'cross':(float(cross) if cross is not None else None),'fit':float((1-A[1])/A[0]),'bins':[(float(b),float(v),int(n)) for b,v,n in med]}
X=np.concatenate(allx); Y=np.concatenate(ally)
bins=np.arange(10,72,2.5); idx=np.digitize(X,bins)
med=[(bins[i-1]+1.25, np.median(Y[idx==i]), (idx==i).sum()) for i in range(1,len(bins)) if (idx==i).sum()>400]
bx=np.array([b for b,v,n in med]); by=np.array([v for b,v,n in med]); cross=None
for i in range(len(bx)-1):
    if (by[i]-1)*(by[i+1]-1)<=0: cross=bx[i]+(1-by[i])*(bx[i+1]-bx[i])/(by[i+1]-by[i]); break
print(f"\nALL LOGS combined: accZ=1G at post-curve throttle {(cross if cross is not None else float('nan')):.1f}%")
print("   bins: "+' '.join(f"{b:.1f}:{v:.2f}" for b,v,n in med))
out['all']={'cross':(float(cross) if cross is not None else None),'bins':[(float(b),float(v),int(n)) for b,v,n in med]}
json.dump(out,open('hover_results.json','w'),indent=1)
