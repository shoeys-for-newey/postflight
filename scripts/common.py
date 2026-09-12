import numpy as np, json, glob, re
from scipy import signal
LOGS = sorted(int(re.search(r'log(\d+)_meta', f).group(1)) for f in glob.glob('log*_meta.json'))
def load(li):
    z = np.load(f'log{li}.npz', allow_pickle=True)
    names = list(z['names']); d = z['data']
    c = {n: d[:, i] for i, n in enumerate(names)}
    meta = json.load(open(f'log{li}_meta.json'))
    t = c['time']/1e6; t = t - t[0]
    fs = 1.0/np.median(np.diff(t))
    return c, t, fs, meta['headers']
def pole_pairs(H):
    return int(H.get('motor_poles', 14))/2
def longest_log():
    best, bn = None, -1
    for li in LOGS:
        n = json.load(open(f'log{li}_meta.json')).get('nframes', 0)
        if n > bn: best, bn = li, n
    return best
def lp(x, fs, fc, order=2):
    b, a = signal.butter(order, fc/(fs/2)); return signal.filtfilt(b, a, x)
def quat_R22_and_aup(c):
    # imuQuaternion[0..2] = x,y,z * 0x7FFF, w = sqrt(1-|v|^2)
    x = c['imuQuaternion[0]']/32767.0; y = c['imuQuaternion[1]']/32767.0; zq = c['imuQuaternion[2]']/32767.0
    w = np.sqrt(np.clip(1 - x*x - y*y - zq*zq, 0, 1))
    ax = c['accSmooth[0]']/2048.0; ay = c['accSmooth[1]']/2048.0; az = c['accSmooth[2]']/2048.0
    r22 = 1 - 2*(x*x + y*y)
    aupA = 2*(x*zq - w*y)*ax + 2*(y*zq + w*x)*ay + r22*az   # R body->earth, row 3
    aupB = 2*(x*zq + w*y)*ax + 2*(y*zq - w*x)*ay + r22*az   # transpose
    return r22, aupA, aupB
