import numpy as np, json, sys, os
from common import *
from curve import out_pct, table
# Assemble page_data.json for the review page from the pipeline outputs + a few extra extracts.
# usage: python pagedata.py CUR_MID,EXPO,HOVER PROP_MID,EXPO,HOVER HOVER_TARGET [turn_t0 turn_dur]
cur=tuple(int(v) for v in sys.argv[1].split(',')); prop=tuple(int(v) for v in sys.argv[2].split(',')); target=float(sys.argv[3])
t0_turn=float(sys.argv[4]) if len(sys.argv)>4 else None; dur_turn=float(sys.argv[5]) if len(sys.argv)>5 else 8.0
li=longest_log()
c,t,fs,H=load(li); thr=c['rcCommand[3]']; idx=np.where(thr>1150)[0]; fl=slice(idx[0],idx[-1]); tf=t[fl]; thp=(thr[fl]-1000)/10
D={}
SC=json.load(open('stepcheck_results.json'))[str(li)]
HK=os.environ.get('HOVER_KEY','tilt<5deg')   # which tilt-restricted hover estimate the page headlines
# a nan crossing (level set too small to bin) would put NaN in the JSON and break JSON.parse
import math
def _finite(k):
    try: return math.isfinite(SC['hover'][k]['cross_body'])
    except Exception: return False
for k in [HK,'tilt<5deg','tilt<10deg']+[k for k in SC['hover'] if k!='all']+['all']:
    if k in SC['hover'] and _finite(k):
        if k!=HK: print(f"hover key '{HK}' has no finite crossing in log {li}; using '{k}'")
        HK=k; break
D['hoverLevel']=[[b,v,n] for b,v,u,n in SC['hover'][HK]['bins']]
D['hoverAll']=[[b,v,n] for b,v,u,n in SC['hover']['all']['bins']]
D['hoverMeas']=round(SC['hover'][HK]['cross_body'],1); D['hoverAllCross']=round(SC['hover']['all']['cross_body'],1)
D['curve']={'cur':[(s,round(out_pct(*cur,s),2)) for s in range(0,101)],'prop':[(s,round(out_pct(*prop,s),2)) for s in range(0,101)],'curLabel':'%d / %d / %d'%cur,'propLabel':'%d / %d / %d'%prop,'curTable':table(*cur),'propTable':table(*prop)}
D['curve']['curHover']=cur[2]; D['curve']['propHover']=prop[2]
ST=json.load(open('step_results.json'))['pooled']
D['step']={a:{'curve':ST[a]['lo']['curve'],'m':ST[a]['lo']['m'],'n':ST[a]['lo']['n']} for a in ('roll','pitch','yaw')}
D['stepFs']=fs
rows=[]
for a in ('roll','pitch','yaw'):
    m=ST[a]['lo']['m']; rows.append([a+' · deconvolution',m['t50'],m['t90'],'%.2f'%m['peak'],m['t_peak'],'%+.0f %%'%m['overshoot_pct'],ST[a]['lo']['n']])
for a in ('roll','pitch','yaw'):
    m=SC['step'][a]; rows.append([a+' · LS-FIR check',m['t50'],m['t90'],'%.2f'%m['peak'],m['t_peak'],'%+.0f %%'%m['overshoot_pct'],'%d s'%(m['n']/1000)])
D['stepRows']=rows
D['stepCheck']={a:{'curve':SC['step'][a]['curve'],'fs':SC['step'][a]['fs']} for a in ('roll','pitch','yaw')}
CD=json.load(open('chart_data.json'))
for k in ('flip_seg','psd_roll','heat','rpm_by_thr'): D[k]=CD[k]
SR=json.load(open('step_results.json'))[str(li)]['offsets']
D['motors']=[SR['calm']['motors'],SR['calm+level+25-45%thr']['motors']]; D['motorN']=[SR['calm']['n'],SR['calm+level+25-45%thr']['n']]
D['motorLabels']=['calm (|gyro| < 60 °/s)','calm + level + 25–45 % throttle']
edges=np.arange(0,65,5); hist,_=np.histogram(thp,edges); pct=100*hist/hist.sum()
D['thrUse']={'labels':['%d–%d'%(edges[i],edges[i+1]) for i in range(len(hist))],'pct':[round(float(v),1) for v in pct]}
amps=c['amperageLatest'][fl]/100
tb=np.arange(0,62.5,2.5); ab=[]
for i in range(len(tb)-1):
    m=(thp>=tb[i])&(thp<tb[i+1])
    if m.sum()>fs*2: ab.append([float(tb[i]+1.25),round(float(np.median(amps[m])),2),int(m.sum())])
D['amps']=ab
# stick estimate under the current curve
xs=np.linspace(0,100,2001); ys=np.array([out_pct(*cur,x) for x in xs]); stick=np.interp(thp,ys,xs)
D['stick']={'median':round(float(np.median(stick)),1),'q1':round(float(np.percentile(stick,25)),1),'q3':round(float(np.percentile(stick,75)),1),'p95':round(float(np.percentile(stick,95)),1),'max':round(float(stick.max()),1)}
# banked-turn segment (yaw setpoint / gyro / I-term / roll angle) at 100 Hz
if t0_turn is None:
    Iy=c['axisI[2]'][fl]; n=int(fs); N=len(Iy)//n; ser=np.abs(Iy[:N*n].reshape(N,n).mean(1)); k=int(np.argmax(ser)); t0_turn=max(0.0,k-dur_turn/2)
x=c['imuQuaternion[0]'][fl]/32767; y=c['imuQuaternion[1]'][fl]/32767; z=c['imuQuaternion[2]'][fl]/32767; w=np.sqrt(np.clip(1-x*x-y*y-z*z,0,1))
roll=np.degrees(np.arctan2(2*(w*x+y*z),1-2*(x*x+y*y)))
s=int(np.argmin(np.abs(tf-tf[0]-t0_turn))); e=min(len(tf),s+int(dur_turn*fs)); step=max(1,int(fs/100))
seg=slice(s,e,step)
D['turn']={'t':[round(float(v),3) for v in (tf[seg]-tf[s])],'sp':[round(float(v),1) for v in lp(c['setpoint[2]'][fl],fs,25)[seg]],'gyro':[round(float(v),1) for v in lp(c['gyroADC[2]'][fl],fs,25)[seg]],'I':[round(float(v),1) for v in c['axisI[2]'][fl][seg]],'roll':[round(float(v),1) for v in roll[seg]],'t0':round(float(tf[s]-tf[0]),1)}
D['hoverKey']=HK
D['meta']={'fs':round(float(fs),1),'craft':str(H.get('Craft name','')),'fw':str(H.get('Firmware revision','')),'dur':round(float(tf[-1]-tf[0]),1),'log':li}
# per-motor RPM + punch-out + saturation extras (need bidir-DShot eRPM; template hides what's missing)
try:
    mot=np.stack([c[f'motor[{i}]'][fl] for i in range(4)],1); motp=(mot-48)/1999*100
    gy=np.stack([c[f'gyroADC[{i}]'][fl] for i in range(3)],1); calm=np.abs(gy).max(1)<60
    if 'eRPM[0]' in c and calm.sum()>fs:
        erpm=np.stack([c[f'eRPM[{i}]'][fl] for i in range(4)],1)*100/pole_pairs(H)
        mr=np.median(erpm[calm],0); mp=np.median(motp[calm],0)
        D['motorRpm']=[round(float(v)) for v in mr]; D['rpmPerPct']=[round(float(a/b)) for a,b in zip(mr,mp)]
    hi=thp>60; sat=motp.max(1)>=98
    if hi.sum()>fs:
        D['sat']={'any_pct':round(float(100*sat.mean()),2),'above60_pct':round(float(100*sat[hi].mean()),1),'which':[round(float(100*np.mean(motp[hi,i]>=98)),1) for i in range(4)]}
    pk=np.where(thp>=99)[0]   # first full-throttle punch, window -1.5 .. +2.0 s at ~100 Hz
    if len(pk):
        s=max(0,int(pk[0])-int(1.5*fs)); e=min(len(tf),s+int(3.5*fs)); st_=max(1,int(round(fs/100))); seg=slice(s,e,st_)
        D['punch']={'t':[round(float(v),3) for v in (tf[seg]-tf[s])],'thr':[round(float(v),1) for v in thp[seg]],'m':[[round(float(v),1) for v in motp[seg,i]] for i in range(4)],'t0':round(float(tf[s]-tf[0]),1)}
except Exception as ex:
    print('extras skipped:',type(ex).__name__,ex)
json.dump(D,open('page_data.json','w'))
import os; print('page_data.json bytes',os.path.getsize('page_data.json'),'| hover level %.1f all %.1f | turn t0 %.1fs | psd range %.0f..%.0f dB'%(D['hoverMeas'],D['hoverAllCross'],D['turn']['t0'],min(min(D['psd_roll']['Pu']),min(D['psd_roll']['Pf']),min(D['psd_roll']['Pd'])),max(max(D['psd_roll']['Pu']),max(D['psd_roll']['Pf']),max(D['psd_roll']['Pd']))))
print('stick',D['stick'],'| curve prop hover stick:',min(D['curve']['prop'],key=lambda p:abs(p[1]-target)),'| cur:',min(D['curve']['cur'],key=lambda p:abs(p[1]-target)))
print('flip_seg sp range',min(D['flip_seg']['sp']),max(D['flip_seg']['sp']),'| turn I range',min(D['turn']['I']),max(D['turn']['I']),'| amps',D['amps'][:3],'...')
