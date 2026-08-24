#  検定FL19  一枚から二つの読みを取る（明るさ と 赤−青）
#
#      網膜はセルごとに数値一つあればよい。どの数値を渡すかは自由なつまみで、
#      画角も担体も変えずに済む。識別へ渡すときは元の画像をそのまま渡す。
#
#      OK なら：赤−青の読みが、明るさでは井戸にならない対象（マガモ）を捉える。
#               かつ色の偏らない明るい部分（上の帯＝反射に相当）を落とす
#      NG なら：どちらの読みでも同じ。分ける理由が無い
#
#      必ず落ちる設定：対象のいない場所から始めた視線が動かないこと。

import math, random, json
import numpy as np
from PIL import Image

W_IMG, H_IMG, NFR = 270, 480, 275
TAU2 = 2 * math.pi

def penrose_cells(spacing, w, h, seed=3):
    rng = random.Random(seed)
    e = [(math.cos(TAU2*j/5), math.sin(TAU2*j/5)) for j in range(5)]
    g = [rng.uniform(0.1, 0.9) for _ in range(4)]; g.append(-sum(g))
    R = math.hypot(w, h)/spacing*0.75; K = int(R)+2
    V = set()
    for j in range(5):
        for l in range(j+1, 5):
            det = e[j][0]*e[l][1] - e[j][1]*e[l][0]
            for kj in range(-K, K+1):
                for kl in range(-K, K+1):
                    a, b = kj-g[j], kl-g[l]
                    px = (a*e[l][1] - b*e[j][1])/det
                    py = (b*e[j][0] - a*e[l][0])/det
                    if px*px + py*py > R*R: continue
                    idx = []
                    for p in range(5):
                        if p == j:   idx.append(kj)
                        elif p == l: idx.append(kl)
                        else:        idx.append(math.ceil(px*e[p][0] + py*e[p][1] + g[p]))
                    V.add((round(sum(idx[p]*e[p][0] for p in range(5)), 5),
                           round(sum(idx[p]*e[p][1] for p in range(5)), 5)))
    return [(x*spacing+w/2, y*spacing+h/2) for x, y in V
            if 0 <= x*spacing+w/2 < w and 0 <= y*spacing+h/2 < h]

RGB = np.array([np.asarray(Image.open("fr/%04d.png" % t).convert("RGB")).astype(np.float32)
                for t in range(1, NFR+1)])
LUMI = 0.299*RGB[..., 0] + 0.587*RGB[..., 1] + 0.114*RGB[..., 2]
RMB  = RGB[..., 0] - RGB[..., 2]

READ = {"明るさ": (LUMI, 140.0), "赤−青": (RMB, 20.0),
        "両方": (np.minimum((LUMI-100.0)*0.5, RMB-20.0), 0.0)}

# ---- 照合（羽ごと）：赤−青で塊を拾う。装置とは無関係なものさし --------
def blobs(field, thr, minpx=30):
    m = field > thr
    seen = np.zeros(m.shape, bool); out = []
    ys, xs = np.nonzero(m)
    for y0, x0 in zip(ys, xs):
        if seen[y0, x0]: continue
        st, pts = [(y0, x0)], []; seen[y0, x0] = True
        while st:
            y, x = st.pop(); pts.append((y, x))
            for dy in range(-3, 4):
                for dx in range(-3, 4):
                    yy, xx = y+dy, x+dx
                    if 0 <= yy < m.shape[0] and 0 <= xx < m.shape[1] \
                       and m[yy, xx] and not seen[yy, xx]:
                        seen[yy, xx] = True; st.append((yy, xx))
        if len(pts) >= minpx:
            P = np.array(pts); out.append((float(P[:,1].mean()), float(P[:,0].mean()), len(pts)))
    return out

BL = [blobs(RMB[t], 22.0) for t in range(NFR)]
def chain(seed):
    tr, cur = [], seed
    for t in range(NFR):
        if not BL[t]: tr.append(None); continue
        b = min(BL[t], key=lambda q: (q[0]-cur[0])**2 + (q[1]-cur[1])**2)
        if math.hypot(b[0]-cur[0], b[1]-cur[1]) > 40: tr.append(None); continue
        cur = (b[0], b[1]); tr.append(cur)
    return tr

print("枚1の赤−青の塊：", [(round(x), round(y), n) for x, y, n in BL[0]])
MALL = chain((215, 315))
DUCK = chain((91, 232))
for nm, tr in (("白いアヒル", DUCK), ("マガモ", MALL)):
    ok = [p for p in tr if p]
    print("  照合 %s：%d/%d 枚  (%.0f,%.0f) → (%.0f,%.0f)"
          % (nm, len(ok), NFR, ok[0][0], ok[0][1], ok[-1][0], ok[-1][1]))

# 明るさで見たときのマガモ／アヒルの明るさ
for nm, tr in (("白いアヒル", DUCK), ("マガモ", MALL)):
    v = [LUMI[t][int(tr[t][1])-5:int(tr[t][1])+5, int(tr[t][0])-5:int(tr[t][0])+5].mean()
         for t in range(0, NFR, 25) if tr[t]]
    print("  %s の明るさ 平均 %.0f（水は 44〜55、しきい値は140）" % (nm, sum(v)/len(v)))

C = np.array(penrose_cells(5.0, W_IMG, H_IMG)); NC = len(C)
CXi, CYi = C[:,0].astype(int), C[:,1].astype(int)

def run(EXC, tau, gamma, start, sig=12.0, k=200.0, sub=32):
    S = np.zeros(NC); x, y = start; vx = vy = 0.0
    dt = 1.0/sub; inv = 1.0/(sig*sig); cut = 3.0*sig
    f = math.exp(-gamma*dt); path = []
    for t in range(NFR):
        L = EXC[t].astype(np.float64); Wt = L - S
        sel = np.abs(Wt) > 0.05; pw, pp = Wt[sel], C[sel]
        for _ in range(sub):
            dx = x-pp[:,0]; dy = y-pp[:,1]
            r2 = dx*dx+dy*dy; m = r2 < cut*cut
            if m.any():
                w = pw[m]*np.exp(-0.5*r2[m]*inv)*inv
                gx, gy = float((w*dx[m]).sum()), float((w*dy[m]).sum())
            else: gx = gy = 0.0
            vx = vx*f - k*gx*dt; vy = vy*f - k*gy*dt
            x = min(max(x+dt*vx, 0.0), W_IMG-1.0); y = min(max(y+dt*vy, 0.0), H_IMG-1.0)
        path.append((x, y)); S += (L-S)/tau
    return path

def sc(p, tr):
    d = [math.hypot(a[0]-b[0], a[1]-b[1]) for a, b in zip(p, tr) if b]
    return sum(d)/len(d), max(d), sum(1 for v in d if v < 20)/len(d)


# 二羽が離れている区間だけで測る
SEP = [t for t in range(NFR) if DUCK[t] and MALL[t]
       and math.hypot(DUCK[t][0]-MALL[t][0], DUCK[t][1]-MALL[t][1]) > 45]
END = max(SEP) + 1 if SEP else NFR
print("二羽が45px以上離れているのは 枚1〜%d（以降は合流するので照合が使えない）" % END)

def onobj(F, thr, tr, r=22):
    """その読みで、対象の上に乗ったセル数（平均）"""
    E = F[:, CYi, CXi] > thr
    n = []
    for t in range(0, END, 5):
        if not tr[t]: continue
        d = (C[:,0]-tr[t][0])**2 + (C[:,1]-tr[t][1])**2 < r*r
        n.append(int((E[t] & d).sum()))
    return sum(n)/len(n)

band = (C[:,1] < 45)
print()
print("%-7s %12s %10s %10s %12s" % ("読み", "励起セル/枚", "アヒルに", "マガモに", "上の帯に"))
for rn, (F, thr) in READ.items():
    E = F[:, CYi, CXi] > thr
    print("%-7s %12.0f %10.1f %10.1f %12.1f"
          % (rn, E.sum(1).mean(), onobj(F, thr, DUCK), onobj(F, thr, MALL),
             float((E[:, band]).sum(1).mean())))

def sc2(p, tr):
    d = [math.hypot(p[t][0]-tr[t][0], p[t][1]-tr[t][1]) for t in range(END) if tr[t]]
    return sum(d)/len(d), max(d), sum(1 for v in d if v < 20)/len(d)

print()
print("%-7s %-10s %-7s %9s %9s %10s" % ("読み", "対象", "τ", "平均ずれ", "最大ずれ", "20px以内"))
for rn, (F, thr) in READ.items():
    EXC = F[:, CYi, CXi] > thr
    for on, tr in (("白いアヒル", DUCK), ("マガモ", MALL)):
        for tau, nm in ((100.0, "100"), (1e9, "止め")):
            m, mx, fr = sc2(run(EXC, tau, 30.0, tr[0]), tr)
            print("%-7s %-10s %-7s %9.1f %9.1f %9.0f%%" % (rn, on, nm, m, mx, 100*fr))
    print()
