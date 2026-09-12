import numpy as np, sys
from common import *
# Load line: at matched RPM bins, how much duty*Vbat each motor needs, as a fraction of its no-load speed (KV*V*duty).
KV=float(sys.argv[1])
R=[];M=[];V=[];A=[]
for li in LOGS:
    c,t,fs,H=load(li); thr=c['rcCommand[3]']; idx=np.where(thr>1150)[0]; fl=slice(idx[0],idx[-1])
    g=np.stack([c[f'gyroADC[{i}]'] for i in range(3)],1); calm=lp(np.abs(g).max(1),fs,4)<120
    for k in range(4):
        rpm=lp(c[f'eRPM[{k}]'],fs,4)*100/pole_pairs(H); d=lp((c[f'motor[{k}]']-48)/1999*100,fs,4)
        R.append(rpm[fl][calm[fl]]); M.append(d[fl][calm[fl]]); V.append((c['vbatLatest']/100)[fl][calm[fl]]); A.append((c['amperageLatest']/100)[fl][calm[fl]])
R=np.concatenate(R);M=np.concatenate(M);V=np.concatenate(V);A=np.concatenate(A)
print(f"{'RPM bin':>8} {'n':>7} {'duty%':>6} {'vbat':>5} {'V*duty':>6} {'no-load':>7} {'fraction':>8} {'amps(total)':>11}")
for b in range(10000,50000,5000):
    m=(R>=b-1500)&(R<b+1500)
    if m.sum()<300: continue
    d=np.median(M[m]); v=np.median(V[m]); nl=KV*v*d/100
    print(f"{b:>8} {m.sum():>7} {d:>6.1f} {v:>5.2f} {v*d/100:>6.2f} {nl:>7.0f} {b/nl:>8.2f} {np.median(A[m]):>11.1f}")
