import math, json, itertools
exec(open('wind_core.py').read())
exec(open('pentagram_neuron.py').read().split('def run_tests')[0])

PHI = (1 + 5 ** 0.5) / 2
F, faces, SC = carrier()
XY = {q: tuple(float(t) for t in U.xy(q)) for q in F}

cells = [[round(XY[q][0], 4), round(XY[q][1], 4), F[q] % 2] for q in F]

# 五芒星30個：中心と、尖り・凹頂点の方位
stars = []
for a, c in faces:
    if abs(a - 2.9389) < 0.01:
        P = [tuple(float(t) for t in U.xy(p)) for p in c]
        cx = sum(p[0] for p in P) / 10
        cy = sum(p[1] for p in P) / 10
        rr = [math.hypot(p[0] - cx, p[1] - cy) for p in P]
        rmax, rmin = max(rr), min(rr)
        tip = min(math.degrees(math.atan2(p[1] - cy, p[0] - cx)) % 72
                  for p, r in zip(P, rr) if abs(r - rmax) < 1e-6)
        notch = min(math.degrees(math.atan2(p[1] - cy, p[0] - cx)) % 72
                    for p, r in zip(P, rr) if abs(r - rmin) < 1e-6)
        stars.append([round(cx, 4), round(cy, 4), round(tip, 3), round(notch, 3)])

# 正十二面体：面ひとつを下にして、面の頂点が方位 90+72k に来るように回す
V, E, FA, FD, D = dodeca()
f0 = FA[0]
c = [sum(V[i][k] for i in f0) / 5 for k in range(3)]
n = [x / math.dist((0, 0, 0), c) for x in c]
# n を −z に向ける回転
def rot_to(v, a, b):
    """a を b に重ねる回転を v に施す（ロドリゲス）。"""
    ax = [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
    s = math.dist((0,0,0), ax); co = sum(a[k]*b[k] for k in range(3))
    if s < 1e-12:
        return list(v) if co > 0 else [-x for x in v]
    k = [x/s for x in ax]; th = math.atan2(s, co)
    ct, st = math.cos(th), math.sin(th)
    d = sum(k[i]*v[i] for i in range(3))
    cr = [k[1]*v[2]-k[2]*v[1], k[2]*v[0]-k[0]*v[2], k[0]*v[1]-k[1]*v[0]]
    return [v[i]*ct + cr[i]*st + k[i]*d*(1-ct) for i in range(3)]

W = [rot_to(p, n, [0, 0, -1]) for p in V]
z0 = min(p[2] for p in W)
W = [[p[0], p[1], p[2] - z0] for p in W]           # 伏せた面を z=0 に
# 面の頂点の方位を 90+72k に揃える
fa = math.degrees(math.atan2(W[f0[0]][1], W[f0[0]][0]))
th = math.radians(90 - fa)
ct, st = math.cos(th), math.sin(th)
W = [[p[0]*ct - p[1]*st, p[0]*st + p[1]*ct, p[2]] for p in W]

rad = max(math.hypot(p[0], p[1]) for p in W)        # 投影の最大半径
hgt = max(p[2] for p in W)
facever = [sorted(f, key=lambda i: math.atan2(
    W[i][1] - sum(W[j][1] for j in f)/5, W[i][0] - sum(W[j][0] for j in f)/5))
    for f in FA]
print("辺長1 で 投影半径 %.6f（φ²/√(2+φ)=%.6f）／高さ %.6f" % (rad, PHI**2/math.sqrt(2+PHI), hgt))

json.dump({"cells": cells, "stars": stars,
           "dv": [[round(x, 5) for x in p] for p in W],
           "df": facever, "seat": 0,
           "rad": round(rad, 6), "hgt": round(hgt, 6),
           "face_r": round(max(math.hypot(W[i][0], W[i][1]) for i in f0), 6)},
          open('solid3d.json', 'w'), separators=(',', ':'))
print(len(cells), len(stars), "→ solid3d.json")
