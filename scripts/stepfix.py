import numpy as np, json
from common import *
from scipy import signal
from scipy.ndimage import gaussian_filter1d
AX=['roll','pitch','yaw']
def deconv(sp,gy,fs,seg_s=2.0,resp_s=0.5,cut=40,min_in=40):
    seg=int(seg_s*fs); resp=int(resp_s*fs); hop=seg//4
    win=signal.windows.tukey(seg,0.25); n=2*seg; f=np.fft.rfftfreq(n,1/fs)
    mask=np.clip(gaussian_filter1d((f<cut).astype(float),8),1e-9,1); reg=0.1/mask
    out={'lo':[],'hi':[]}; w={'lo':[],'hi':[]}
    for s in range(0,len(sp)-seg,hop):
        raw=sp[s:s+seg]; mx=np.abs(raw).max()
        if mx<min_in: continue
        X=np.fft.rfft(raw*win,n); Y=np.fft.rfft(gy[s:s+seg]*win,n)
        st=np.cumsum(np.fft.irfft(Y*np.conj(X)/(np.abs(X)**2+reg),n)[:resp])
        if not np.isfinite(st).all(): continue
        k='hi' if mx>500 else 'lo'; out[k].append(st); w[k].append(mx)
    res={}
    for k in out:
        if out[k]:
            L=np.array(out[k]); W=np.array(w[k]); res[k]=((L*W[:,None]).sum(0)/W.sum(),len(L))
    return res
def metrics(st,fs):
    tt=np.arange(len(st))/fs*1000; lim=int(0.2*fs)
    pk=st[:lim].max(); tpk=tt[st[:lim].argmax()]
    def tc(v):
        i=np.where(st>=v)[0]; return round(float(tt[i[0]]),1) if len(i) else None
    ss=float(np.mean(st[int(0.3*fs):int(0.5*fs)]))
    bad=np.where(np.abs(st[:int(0.3*fs)]-1)>0.05)[0]; settle=round(float(tt[bad[-1]]),1) if len(bad) else 0.0
    return {'t50':tc(0.5),'t90':tc(0.9),'t100':tc(1.0),'peak':round(float(pk),3),'t_peak':round(float(tpk),1),'overshoot_pct':round(float((pk-1)*100),1),'settle5_ms':settle,'ss':round(ss,3)}
OUT={}; pool={a:{'lo':[],'hi':[]} for a in AX}
for li in LOGS:
    c,t,fs,H=load(li); thr=c['rcCommand[3]']; idx=np.where(thr>1150)[0]; fl=slice(idx[0],idx[-1])
    OUT[li]={}
    print("== LOG %d =="%li)
    for i,a in enumerate(AX):
        res=deconv(c['setpoint[%d]'%i][fl],c['gyroADC[%d]'%i][fl],fs)
        OUT[li][a]={}
        for k,(st,n) in res.items():
            m=metrics(st,fs); OUT[li][a][k]={'n':n,'m':m,'curve':[round(float(v),4) for v in st[:int(0.25*fs)]]}
            pool[a][k].append((st,n))
            print("  %-5s %-2s n=%3d  t50 %5s  t90 %5s  t100 %5s  peak %.2f @%3.0fms  overshoot %+5.1f%%  settle5 %5.1fms  ss %.2f"%(a,k,n,m['t50'],m['t90'],m['t100'],m['peak'],m['t_peak'],m['overshoot_pct'],m['settle5_ms'],m['ss']))
    # signed I-term offsets + motor balance when calm & level
    g=np.stack([c['gyroADC[%d]'%i][fl] for i in range(3)],1); gmax=np.abs(g).max(1)
    r22,_,_=quat_R22_and_aup(c); tilt=np.degrees(np.arccos(np.clip(r22[fl],-1,1)))
    thp=(thr[fl]-1000)/10
    mot=np.stack([c['motor[%d]'%i][fl] for i in range(4)],1); motp=(mot-48)/1999*100
    I=np.stack([c['axisI[%d]'%i][fl] for i in range(3)],1)
    calm=gmax<60; lvl=calm&(tilt<6)&(thp>25)&(thp<45)
    OUT[li]['offsets']={}
    for nm,mk in (('calm',calm),('calm+level+25-45%thr',lvl)):
        if mk.sum()<200: continue
        mm=motp[mk].mean(0); fr=(mm[1]+mm[3])/2-(mm[0]+mm[2])/2; lr=(mm[2]+mm[3])/2-(mm[0]+mm[1])/2; dg=(mm[0]+mm[3])/2-(mm[1]+mm[2])/2
        Im=I[mk].mean(0)
        OUT[li]['offsets'][nm]={'n':int(mk.sum()),'motors':[round(float(v),1) for v in mm],'front_minus_rear':round(float(fr),1),'left_minus_right':round(float(lr),1),'diag14_minus_23':round(float(dg),1),'I_mean':[round(float(v),1) for v in Im]}
        print("  %-22s n=%6d motors RR/FR/RL/FL %s  front-rear %+.1f  left-right %+.1f  diag(1+4)-(2+3) %+.1f | I mean R/P/Y %s"%(nm,mk.sum(),np.round(mm,1),fr,lr,dg,np.round(Im,1)))
print("== POOLED (all logs, n-weighted) ==")
OUT['pooled']={}
for a in AX:
    OUT['pooled'][a]={}
    for k in ('lo','hi'):
        if pool[a][k]:
            L=np.array([s for s,n in pool[a][k]]); W=np.array([n for s,n in pool[a][k]],float); st=(L*W[:,None]).sum(0)/W.sum(); m=metrics(st,fs)
            OUT['pooled'][a][k]={'n':int(W.sum()),'m':m,'curve':[round(float(v),4) for v in st[:int(0.25*fs)]]}
            print("  %-5s %-2s n=%3d  t50 %5s  t90 %5s  t100 %5s  peak %.2f @%3.0fms  overshoot %+5.1f%%  settle5 %5.1fms  ss %.2f"%(a,k,int(W.sum()),m['t50'],m['t90'],m['t100'],m['peak'],m['t_peak'],m['overshoot_pct'],m['settle5_ms'],m['ss']))
json.dump(OUT,open('step_results.json','w'))
