import numpy as np, json
from common import *
from scipy import signal
CD={}
LONGEST=longest_log()
for li in LOGS:
    c,t,fs,H=load(li); thr=c['rcCommand[3]']; idx=np.where(thr>1150)[0]; fl=slice(idx[0],idx[-1])
    tt=t[fl]; thp=(thr[fl]-1000)/10; gr=c['gyroADC[0]'][fl]; gu=c['gyroUnfilt[0]'][fl]; D=c['axisD[0]'][fl]
    b,a=signal.butter(4,[70/(fs/2),100/(fs/2)],'band'); x=signal.filtfilt(b,a,gr)
    W=int(0.25*fs); n=len(x)//W
    env=np.sqrt((x[:n*W].reshape(n,W)**2).mean(1)); tw=tt[:n*W].reshape(n,W).mean(1); th=thp[:n*W].reshape(n,W).mean(1)
    hot=env>15
    print("LOG %d: 70-100Hz roll band RMS: median %.1f, 90pct %.1f, max %.1f deg/s; windows >15deg/s: %.1f%% of flight"%(li,np.median(env),np.percentile(env,90),env.max(),100*hot.mean()))
    for lo,hi in ((0,20),(20,30),(30,40),(40,55),(55,100)):
        m=(th>=lo)&(th<hi)
        if m.sum()>4: print("   thr %2d-%3d%%: band RMS median %.1f, >15: %.0f%% (n=%d)"%(lo,hi,np.median(env[m]),100*hot[m].mean(),m.sum()))
    top=np.argsort(env)[::-1][:6]
    print("   top windows: "+', '.join("t=%.1fs %.0fdeg/s thr%.0f%%"%(tw[k],env[k],th[k]) for k in top))
    if li==LONGEST:
        k=top[0]; s=int(np.argmin(np.abs(tt-tw[k])))-int(0.15*fs); e=s+int(0.3*fs)
        CD['osc_seg']={'t':[round(float(v),4) for v in (tt[s:e]-tt[s])],'gyro':[round(float(v),1) for v in gr[s:e]],'unfilt':[round(float(v),1) for v in gu[s:e]],'dterm':[round(float(v),1) for v in D[s:e]],'thr':round(float(thp[s:e].mean()),1)}
    if li==LONGEST:
        sp=c['setpoint[0]'][fl]; k=int(np.argmax(np.abs(sp))); s=max(0,k-int(0.35*fs)); e=s+int(0.9*fs)
        CD['flip_seg']={'t':[round(float(v),4) for v in (tt[s:e:2]-tt[s])],'sp':[round(float(v),1) for v in sp[s:e:2]],'gyro':[round(float(v),1) for v in gr[s:e:2]],'thr':[round(float(v),1) for v in thp[s:e:2]]}
    h=np.load('heat%d.npz'%li); f=h['f']; k86=int(np.argmin(np.abs(f-86)))
    if li==LONGEST:
        print("   D-term roll power at %.0fHz by throttle bin (dB): "%f[k86]+' '.join("%.0f%%:%.0f"%(h['tb'][i],10*np.log10(h['d0'][i,k86]+1e-9)) for i in range(len(h['cnt'])) if h['cnt'][i]>3))
# chart data from the longest log
h=np.load('heat%d.npz'%LONGEST); f=h['f']; keep=f<=1000
def db(M): return [[round(float(v),1) for v in row] for row in 10*np.log10(M[:,keep]+1e-9)]
CD['heat']={'f':[round(float(v),1) for v in f[keep]],'tb':[float(v) for v in h['tb'][:-1]],'cnt':[int(v) for v in h['cnt']],'gu0':db(h['gu0']),'gy0':db(h['gy0']),'d0':db(h['d0'])}
P=json.load(open('pid_results.json'))[str(LONGEST)]['noise']
for a in ('roll','pitch','yaw'):
    ff=np.array(P[a]['f']); m=ff<=1000
    CD['psd_'+a]={'f':[round(float(v),1) for v in ff[m][::2]],'Pu':[round(float(v),1) for v in 10*np.log10(np.array(P[a]['Pu'])[m][::2]+1e-9)],'Pf':[round(float(v),1) for v in 10*np.log10(np.array(P[a]['Pf'])[m][::2]+1e-9)]}
    if 'Pd' in P[a]: CD['psd_'+a]['Pd']=[round(float(v),1) for v in 10*np.log10(np.array(P[a]['Pd'])[m][::2]+1e-9)]
mf=np.array(P['motor0']['f']); m=mf<=1000
CD['psd_motor']={'f':[round(float(v),1) for v in mf[m][::2]],'P':[round(float(v),1) for v in 10*np.log10(np.array(P['motor0']['P'])[m][::2]+1e-9)]}
CD['rpm_by_thr']=json.load(open('pid_results.json'))[str(LONGEST)]['rpm_by_thr']
CD['hover']=json.load(open('hover_results.json')); CD['curve']=json.load(open('curve.json'))
json.dump(CD,open('chart_data.json','w'))
import os; print("chart_data.json bytes:",os.path.getsize('chart_data.json'))
