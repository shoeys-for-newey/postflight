import numpy as np, time, json, sys
from orangebox import Parser
from orangebox.reader import Reader
# usage: python decode.py <path-to.bbl>   (writes logN.npz + logN_meta.json into the CWD)
# Each contiguous flight segment becomes its own logN (pauses from the BLACKBOX switch split a log).
if len(sys.argv) < 2:
    sys.exit("usage: python decode.py <path-to.bbl>   (run from inside a run folder)")
src = sys.argv[1]
nlogs = Reader(src, None).log_count   # log_index=None: just count logs, don't parse yet
print("logs in file:", nlogs)
out_i = 0
for li in range(1, nlogs+1):
    t0=time.time()
    try:   # a corrupt/near-empty log (e.g. a header-only aborted arm) can make orangebox fail here,
           # so give each log its own parser instead of set_log_index on a shared one
        p = Parser.load(src, li)
        names = list(p.field_names)
    except Exception as ex:
        print(f"log {li}: unreadable ({type(ex).__name__}: {ex}) - skipped"); continue
    nnum = names.index('flightModeFlags') if 'flightModeFlags' in names else len(names)
    rows=[]; flags=[]; state=[]
    trunc = None
    try:
        for f in p.frames():
            d = f.data
            rows.append(d[:nnum])
            flags.append(str(d[nnum]) if len(d)>nnum else '')
            state.append(str(d[nnum+1]) if len(d)>nnum+1 else '')
    except Exception as ex:   # orangebox raises on a truncated/corrupt tail; keep what we have
        trunc = f"{type(ex).__name__}: {ex}"
    if not rows:
        print(f"log {li}: NO FRAMES decoded ({trunc})"); continue
    arr = np.array(rows, dtype=np.float64); flags=np.array(flags); state=np.array(state)
    ev = []
    try:
        for e in p.events:
            ev.append({'name': str(getattr(e,'name',e)), 'data': str(getattr(e,'data','')), 'time': str(getattr(e,'time',''))})
    except Exception as ex:
        ev=[{'err':str(ex)}]
    hdr = {k:(v if isinstance(v,(int,float,str)) else str(v)) for k,v in p.headers.items()}
    tcol = names.index('time'); vbcol = names.index('vbatLatest'); thcol = names.index('rcCommand[3]')
    # drop leading pre-I-frame garbage (time==0 / vbat==0 rows orangebox emits before the first I frame)
    tt = arr[:,tcol]/1e6; dt = np.diff(tt)
    junk = np.where(((dt > 1.0) | (dt < 0) | (arr[:-1,tcol] <= 0) | (arr[:-1,vbcol] <= 0)) & (np.arange(len(dt)) < 2000))[0]
    start = int(junk[-1]) + 1 if len(junk) else 0
    if start >= len(arr): print(f"log {li}: no valid frames"); continue
    if start > 0: print(f"log {li}: dropping {start} leading garbage frames")
    arr=arr[start:]; flags=flags[start:]; state=state[start:]
    tt = arr[:,tcol]/1e6; dt = np.diff(tt)
    cuts = np.where((dt > 1.0) | (dt < 0))[0]
    starts = np.r_[0, cuts+1]; ends = np.r_[cuts+1, len(arr)]
    for k,(s,e) in enumerate(zip(starts,ends)):
        dur = tt[e-1]-tt[s]
        if dur < 5.0:
            if dur > 0.1: print(f"log {li} seg {k}: skipped ({e-s} frames, {dur:.2f}s)")
            continue
        out_i += 1
        a = arr[s:e]
        np.savez_compressed(f'log{out_i}.npz', data=a, names=np.array(names[:nnum]), flags=flags[s:e], state=state[s:e])
        json.dump({'headers':hdr,'events':ev,'truncated':trunc,'source_log':li,'segment':k,'nframes':int(e-s),'duration_s':float(dur)}, open(f'log{out_i}_meta.json','w'), indent=1)
        t = a[:,tcol]/1e6; d2 = np.diff(t); thr = a[:,thcol]; vb = a[:,vbcol]
        print(f"log{out_i} (file log {li} seg {k}): {len(a)} frames, {dur:.1f}s, median dt {np.median(d2)*1e6:.0f}us, max gap {d2.max()*1e3:.1f}ms, "
              f"thr {thr.min():.0f}-{thr.max():.0f}, vbat {vb.min()/100:.2f}-{vb.max()/100:.2f}V, decode {time.time()-t0:.1f}s"
              + (f"  [TRUNCATED TAIL: {trunc}]" if trunc and k==len(starts)-1 else ""))
    print("  events:", ev[:12])
    print("  flag values:", sorted(set(flags))[:10])
