import numpy as np, json
from common import *
from scipy import signal
AX=['roll','pitch','yaw']
R={}
def welch(x,fs,nper=2048):
    f,P=signal.welch(x,fs,nperseg=nper,noverlap=nper//2,window='hann'); return f,P
def step_resp(sp,gy,fs,seg_s=2.0,resp_s=0.5,min_in=40,hi_split=500):
    seg=int(seg_s*fs); resp=int(resp_s*fs); hop=seg//4
    win=signal.windows.tukey(seg,0.25)
    lo=[];hi=[];wl=[];wh=[]
    for s in range(0,len(sp)-seg,hop):
        x=sp[s:s+seg]; y=gy[s:s+seg]; mx=np.abs(x).max()
        if mx<min_in: continue
        n=2*seg; X=np.fft.rfft(x*win,n); Y=np.fft.rfft(y*win,n)
        sn=0.001*np.max(np.abs(X)**2)
        h=np.fft.irfft(Y*np.conj(X)/(np.abs(X)**2+sn),n)[:resp]
        st=np.cumsum(h)
        if not np.isfinite(st).all(): continue
        (hi if mx>hi_split else lo).append(st); (wh if mx>hi_split else wl).append(mx)
    def avg(L,W):
        if not L: return None
        L=np.array(L); W=np.array(W); return (L*W[:,None]).sum(0)/W.sum()
    return avg(lo,wl), avg(hi,wh), len(lo), len(hi)
def metrics(st,fs):
    if st is None: return None
    tt=np.arange(len(st))/fs*1000
    ss=np.mean(st[int(0.3*fs):int(0.5*fs)])
    pk=st[:int(0.15*fs)].max(); tpk=tt[st[:int(0.15*fs)].argmax()]
    def tcross(v):
        i=np.where(st>=v)[0]; return float(tt[i[0]]) if len(i) else None
    return {'t50':tcross(0.5),'t90':tcross(0.9),'t100':tcross(1.0),'peak':float(pk),'t_peak':float(tpk),'overshoot_pct':float((pk-1)*100),'ss':float(ss)}
def peaks_report(f,P,fmin,fmax,ratio=6):
    m=(f>=fmin)&(f<=fmax); ff=f[m]; PP=P[m]
    base=signal.medfilt(np.log10(PP+1e-20),kernel_size=31)
    idx,_=signal.find_peaks(np.log10(PP+1e-20)-base,height=np.log10(ratio))
    return [(float(ff[i]),float(10*np.log10(PP[i]))) for i in idx]
for li in LOGS:
    c,t,fs,H=load(li)
    thr=c['rcCommand[3]']; thp=(thr-1000)/10
    idx=np.where(thr>1150)[0]; fl=slice(idx[0],idx[-1])
    r={'dur':float(t[fl][-1]-t[fl][0])}
    sp=[c['setpoint[%d]'%i][fl] for i in range(3)]; gy=[c['gyroADC[%d]'%i][fl] for i in range(3)]; gu=[c['gyroUnfilt[%d]'%i][fl] for i in range(3)]
    P=[c['axisP[%d]'%i][fl] for i in range(3)]; I=[c['axisI[%d]'%i][fl] for i in range(3)]; D=[c['axisD[%d]'%i][fl] for i in range(2)]; F=[c['axisF[%d]'%i][fl] for i in range(3)]
    mot=np.stack([c['motor[%d]'%i][fl] for i in range(4)],1); motp=(mot-48)/1999*100
    erpm=np.stack([c['eRPM[%d]'%i][fl] for i in range(4)],1)*100/pole_pairs(H)
    thpf=thp[fl]
    print("\n================ LOG %d (%.1fs in flight) ================"%(li,r['dur']))
    r['step']={}
    for i in range(3):
        lo,hi,nl,nh=step_resp(sp[i],gy[i],fs)
        ml=metrics(lo,fs); mh=metrics(hi,fs)
        r['step'][AX[i]]={'low':ml,'high':mh,'n_low':nl,'n_high':nh,'curve_low':None if lo is None else lo[:int(0.3*fs)].tolist(),'curve_high':None if hi is None else hi[:int(0.3*fs)].tolist()}
        for nm,mm,n in (('low<500',ml,nl),('high>500',mh,nh)):
            if mm: print(" STEP %-5s %-9s n=%3d: t50=%sms t90=%sms t100=%sms peak=%.2f@%.0fms overshoot=%+.0f%% ss=%.2f"%(AX[i],nm,n,mm['t50'],mm['t90'],mm['t100'],mm['peak'],mm['t_peak'],mm['overshoot_pct'],mm['ss']))
    r['track']={}
    for i in range(3):
        e=sp[i]-gy[i]; act=np.abs(sp[i])>50; qu=np.abs(sp[i])<20
        r['track'][AX[i]]={'rms_active':float(np.sqrt(np.mean(e[act]**2))),'rms_quiet':float(np.sqrt(np.mean(e[qu]**2))),'sp_rms_active':float(np.sqrt(np.mean(sp[i][act]**2))),'sp_max':float(np.abs(sp[i]).max())}
        d=r['track'][AX[i]]
        print(" TRACK %-5s: err RMS active %.1f (sp RMS %.0f, max %.0f deg/s), quiet-stick err RMS %.1f deg/s"%(AX[i],d['rms_active'],d['sp_rms_active'],d['sp_max'],d['rms_quiet']))
    r['terms']={}
    for i in range(3):
        d={'P_rms':float(np.sqrt(np.mean(P[i]**2))),'I_rms':float(np.sqrt(np.mean(I[i]**2))),'I_max':float(np.abs(I[i]).max()),'I_gt150_pct':float(100*np.mean(np.abs(I[i])>150)),'F_rms':float(np.sqrt(np.mean(F[i]**2)))}
        s=" TERMS %-5s: P rms %.0f I rms %.0f (max %.0f, >150 %.1f%%) F rms %.0f"%(AX[i],d['P_rms'],d['I_rms'],d['I_max'],d['I_gt150_pct'],d['F_rms'])
        if i<2:
            d['D_rms']=float(np.sqrt(np.mean(D[i]**2)))
            b,a=signal.butter(2,80/(fs/2),'high'); Dh=signal.filtfilt(b,a,D[i]); d['D_hf_rms']=float(np.sqrt(np.mean(Dh**2)))
            s+=" D rms %.0f D>80Hz rms %.1f"%(d['D_rms'],d['D_hf_rms'])
        r['terms'][AX[i]]=d; print(s)
    r['noise']={}
    for i in range(3):
        f,Pu=welch(gu[i],fs); _,Pf=welch(gy[i],fs)
        def band(Pp,a,b):
            m=(f>=a)&(f<b); return float(10*np.log10(np.mean(Pp[m])+1e-20))
        pk_u=peaks_report(f,Pu,60,1000); pk_f=peaks_report(f,Pf,60,1000)
        r['noise'][AX[i]]={'f':f.tolist(),'Pu':Pu.tolist(),'Pf':Pf.tolist(),'peaks_unfilt':pk_u[:8],'peaks_filt':pk_f[:8],
            'unf_100_300':band(Pu,100,300),'unf_300_600':band(Pu,300,600),'unf_600_1000':band(Pu,600,1000),
            'flt_100_300':band(Pf,100,300),'flt_300_600':band(Pf,300,600),'flt_600_1000':band(Pf,600,1000)}
        print(" NOISE %-5s: unfilt dB 100-300:%.0f 300-600:%.0f 600-1k:%.0f | filt 100-300:%.0f 300-600:%.0f 600-1k:%.0f"%(AX[i],band(Pu,100,300),band(Pu,300,600),band(Pu,600,1000),band(Pf,100,300),band(Pf,300,600),band(Pf,600,1000)))
        print("    unfilt peaks: "+', '.join("%.0fHz(%.0fdB)"%(a,b) for a,b in pk_u[:8]))
        print("    filt   peaks: "+', '.join("%.0fHz(%.0fdB)"%(a,b) for a,b in pk_f[:8]))
    r['osc']={}
    for i in range(3):
        f,Pf=welch(gy[i],fs,4096); _,Ps=welch(sp[i],fs,4096)
        m=(f>=8)&(f<=200)
        ratio=Pf[m]/(Ps[m]+1e-9)
        base=signal.medfilt(np.log10(Pf[m]+1e-20),51)
        pk,_=signal.find_peaks(np.log10(Pf[m]+1e-20)-base,height=np.log10(4))
        cand=[(float(f[m][k]),float(10*np.log10(Pf[m][k])),float(ratio[k])) for k in pk if ratio[k]>2]
        r['osc'][AX[i]]=cand[:6]
        print(" OSC %-5s (8-200Hz gyro peaks not in setpoint): "%AX[i]+(', '.join("%.0fHz(%.0fdB,x%.0f)"%(a,b,cc) for a,b,cc in cand[:6]) or 'none'))
    for i in range(2):
        f,Pd=welch(D[i],fs); r['noise'][AX[i]]['Pd']=Pd.tolist()
        pk=peaks_report(f,Pd,40,1000)
        print(" DTERM %-5s peaks: "%AX[i]+', '.join("%.0fHz(%.0fdB)"%(a,b) for a,b in pk[:8]))
    f,Pm=welch(motp[:,0],fs); r['noise']['motor0']={'f':f.tolist(),'P':Pm.tolist()}
    b,a=signal.butter(2,100/(fs/2),'high'); mh=signal.filtfilt(b,a,motp,axis=0)
    r['motor_hf_rms']=[float(np.sqrt(np.mean(mh[:,k]**2))) for k in range(4)]
    print(" MOTOR >100Hz noise RMS (%% of output): "+' '.join("%.2f"%v for v in r['motor_hf_rms'])+"  motor0 peaks: "+', '.join("%.0fHz(%.0fdB)"%(a,b) for a,b in peaks_report(f,Pm,40,1000)[:6]))
    r['delay_ms']={}
    for i in range(3):
        b,a=signal.butter(2,[10/(fs/2),120/(fs/2)],'band'); xu=signal.filtfilt(b,a,gu[i]); xf=signal.filtfilt(b,a,gy[i])
        L=int(0.01*fs); N=len(xu)-L-1
        cc=[np.dot(xu[:N],xf[k:k+N]) for k in range(L)]
        k=int(np.argmax(cc)); r['delay_ms'][AX[i]]=k/fs*1000
    print(" FILTER DELAY (10-120Hz, unfilt->filt): roll %.2fms pitch %.2fms yaw %.2fms"%(r['delay_ms']['roll'],r['delay_ms']['pitch'],r['delay_ms']['yaw']))
    sat_hi=np.mean(motp.max(1)>=98)*100; sat_lo=np.mean(motp.min(1)<=0.5)*100
    r['motor']={'sat_hi_pct':float(sat_hi),'sat_lo_pct':float(sat_lo),'mean':[float(v) for v in motp.mean(0)],'rpm_min_inflight':float(np.percentile(erpm.min(1),1)),'rpm_median':float(np.median(erpm)),'rpm_max':float(erpm.max())}
    hb=thpf>60
    print(" MOTORS: any motor >=98%%: %.2f%% of time (at thr>60%%: %.1f%%); any motor at floor: %.2f%%; per-motor mean %% %s; RPM 1pct-min %.0f median %.0f max %.0f"%(sat_hi,(np.mean(motp[hb].max(1)>=98)*100 if hb.any() else 0),sat_lo,np.round(motp.mean(0),1),r['motor']['rpm_min_inflight'],r['motor']['rpm_median'],r['motor']['rpm_max']))
    gmax=np.abs(np.stack(gy,1)).max(1); q=(gmax<60)
    r['motor']['calm_mean']=[float(v) for v in motp[q].mean(0)]
    print("   motor mean%% when calm (n=%d): %s  (BF order: 1=RR,2=FR,3=RL,4=FL)"%(q.sum(),np.round(motp[q].mean(0),1)))
    nper=256; hop=256; tb=np.arange(0,102.5,2.5); nfb=nper//2+1
    heat={k:np.zeros((len(tb)-1,nfb)) for k in ('gu0','gy0','d0','gu1','gy1','m0')}; cnt=np.zeros(len(tb)-1)
    win=np.hanning(nper)
    srcs={'gu0':gu[0],'gy0':gy[0],'d0':D[0],'gu1':gu[1],'gy1':gy[1],'m0':motp[:,0]}
    for s in range(0,len(thpf)-nper,hop):
        bb=np.digitize(thpf[s:s+nper].mean(),tb)-1
        if bb<0 or bb>=len(cnt): continue
        cnt[bb]+=1
        for k,v in srcs.items():
            heat[k][bb]+=np.abs(np.fft.rfft(v[s:s+nper]*win))**2
    fh=np.fft.rfftfreq(nper,1/fs)
    np.savez_compressed('heat%d.npz'%li,**{k:(v/np.maximum(cnt,1)[:,None]) for k,v in heat.items()},f=fh,tb=tb,cnt=cnt)
    rb=[float(np.median(erpm.mean(1)[(thpf>=tb[i])&(thpf<tb[i+1])])) if cnt[i]>3 else None for i in range(len(tb)-1)]
    r['rpm_by_thr']=rb
    print(" RPM by throttle bin: "+' '.join("%.0f%%:%.0fHz"%(tb[i],rb[i]/60) for i in range(len(rb)) if rb[i]))
    R[li]=r
json.dump(R,open('pid_results.json','w'))
print("\nsaved pid_results.json")
