#  網の沈みを読む
#      各セル q の沈み w(q) ＝ q から半径 r 以内で点いているセルの枚数（数えるだけ）。
#      いちばん深く沈んだところへ視線を移す。閾値は置かない。
#      沈みの深さ＝対象がいるかどうか／沈んだ範囲＝大きさ／沈みの前後の偏り＝向きと速さ
exec(open('ballistic.py').read().split('print("■ 投げられたものを追う')[0])

def line(T, sp, ang=0.6):
    return [(O[0]-0.5*T*sp*math.cos(ang)+t*sp*math.cos(ang),
             O[1]-0.5*T*sp*math.sin(ang)+t*sp*math.sin(ang)) for t in range(T)]

def tri_at(x, y, R):
    V = [(x+R*math.cos(math.radians(90+120*k)), y+R*math.sin(math.radians(90+120*k))) for k in range(3)]
    out = set()
    for k in range(3):
        x0,y0 = V[k]; x1,y1 = V[(k+1)%3]
        for i in range(8):
            t = i/8; px, py = x0+t*(x1-x0), y0+t*(y1-y0)
            out.add(min(CELLS, key=lambda q: (XYC[q][0]-px)**2 + (XYC[q][1]-py)**2))
    return out

def scene(pts, sd, R=6.0):
    back = set(random.Random(sd).sample(CELLS, 100))
    return [back | set(random.Random(50+7*t).sample(CELLS, 8)) | tri_at(x, y, R)
            for t, (x, y) in enumerate(pts)]

def sink(E, G, Rnet, r=2.7):
    """網の中の各セルの沈み。返すのは (最大の沈み, 深いところの重心, 深いセルの枚数)"""
    P = [XYC[q] for q in CELLS
         if math.hypot(XYC[q][0]-G[0], XYC[q][1]-G[1]) <= Rnet and q in E]
    if not P:
        return 0, None, 0
    w = []
    for a in P:
        c = sum(1 for b in P if (a[0]-b[0])**2 + (a[1]-b[1])**2 <= r*r)
        w.append(c)
    mx = max(w)
    deep = [P[i] for i in range(len(P)) if w[i] >= mx]
    n_half = sum(1 for x in w if x >= mx/2)
    return mx, (sum(a[0] for a in deep)/len(deep), sum(a[1] for a in deep)/len(deep)), n_half

def follow(Es, pts, Rnet=6.47):
    G = pts[0]; d = []; deep = []
    for t, E in enumerate(Es):
        mx, c, _ = sink(E, G, Rnet)
        if c is not None:
            G = c
        deep.append(mx)
        d.append(math.hypot(G[0]-pts[t][0], G[1]-pts[t][1]))
    return d, sum(deep)/len(deep)

print("■ 沈んだところへ移る（網の半径6.47・総移動距離30に揃える・種3本）\n")
print("   一歩   捉えている枚  平均のずれ  対象ありの沈み  対象なしの沈み")
for i in (1, 5, 8, 12, 16, 24, 40):
    sp = 0.25*i; T = max(10, int(30/sp))
    c = 0.0; m = 0.0; dp = 0.0; dp0 = 0.0
    for sd in (1, 2, 3):
        pts = line(T, sp)
        d, deep = follow(scene(pts, sd), pts)
        c += sum(1 for x in d if x <= 4.236)/len(d)/3; m += sum(d)/len(d)/3; dp += deep/3
        back = set(random.Random(sd).sample(CELLS, 100))
        E0 = [back | set(random.Random(50+7*t).sample(CELLS, 8)) for t in range(T)]
        dp0 += sum(sink(E, O, 6.47)[0] for E in E0)/T/3
    print("   %4.2f      %.3f        %5.1f        %4.1f            %4.1f" % (sp, c, m, dp, dp0))

print("\n■ 沈んだ範囲は大きさを語るか（一歩1.0・網13.0・沈みの半分以上のセル枚数）")
print("   対象の外接半径   沈んだ範囲の枚数   最大の沈み")
for R in (4.0, 6.0, 9.0, 13.0):
    n = 0.0; mx = 0.0
    for sd in (1, 2, 3):
        pts = line(20, 1.0)
        Es = scene(pts, sd, R)
        for t, E in enumerate(Es):
            a, c, nh = sink(E, pts[t], 13.0)
            n += nh/len(Es)/3; mx += a/len(Es)/3
    print("      %4.1f            %5.1f            %4.1f" % (R, n, mx))
