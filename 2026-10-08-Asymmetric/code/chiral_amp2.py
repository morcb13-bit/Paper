# 留保2点の解決（整数と加算・比較だけ）
# 基準は chiral_amp.py と同じ：
#   F1: ある T 以降ずっと |D| >= 2|D(0)|
#   F2: 最終で 20*同側 >= 19*占有、かつ D != 0
#   対照: 対称初期で D=0、足し込みなしで F1 NG

def judge(D0, hist):
    a0 = D0 if D0 >= 0 else -D0
    Ds = [h[2] for h in hist]
    f1 = False
    if a0 > 0:
        for t in range(len(Ds)):
            if all((d if d >= 0 else -d) >= a0 + a0 for d in Ds[t:]):
                f1 = True; break
    R, L, D = hist[-1]
    occ = R + L; same = R if D >= 0 else L
    f2 = occ > 0 and 20 * same >= 19 * occ and D != 0
    return f1, f2

# ---------- 留保1：器で、組を除く量を一歩あたり cap 組に制限 ----------
def mixed_cap(R, L, food, cap, grow=True, steps=60000):
    hist = []
    for _ in range(steps):
        if grow:
            gr = gl = 0
            # 1個ずつ交互に配る。多い側から先に配る（型を見ず数の比較だけ。鏡像で同じ結果になる）
            first_R = R >= L
            while food > 0 and (gr < R or gl < L):
                if first_R:
                    if gr < R and food > 0: gr += 1; food += -1
                    if gl < L and food > 0: gl += 1; food += -1
                else:
                    if gl < L and food > 0: gl += 1; food += -1
                    if gr < R and food > 0: gr += 1; food += -1
            R, L = R + gr, L + gl
        m = R if R < L else L
        if m > cap: m = cap
        R, L = R + -m, L + -m
        hist.append((R, L, R + -L))
    return hist

print("=== 留保1：組を除く量を制限した器 ===")
for cap in (1, 10, 100, 1000):
    for (R, L, tag) in ((50, 50, "対称"), (51, 50, "+1"), (50, 51, "-1")):
        h = mixed_cap(R, L, 100000, cap)
        f1, f2 = judge(R + -L, h)
        # 飽和に達した最初の一歩
        tsat = next((t for t, (r, l, d) in enumerate(h)
                     if (r + l) > 0 and d != 0 and 20 * (r if d > 0 else l) >= 19 * (r + l)), None)
        print(f"cap={cap:5d} {tag:4s} F1={'OK' if f1 else 'NG'} F2={'OK' if f2 else 'NG'} "
              f"飽和の一歩={tsat} 最終={h[-1]}")
# 負の対照：増えない器で組だけ除く
for (R, L, tag) in ((51, 50, "+1"),):
    h = mixed_cap(R, L, 100000, 10, grow=False)
    f1, f2 = judge(R + -L, h)
    print(f"足し込みなし cap=10 {tag} F1={'OK' if f1 else 'NG'} F2={'OK' if f2 else 'NG'} 最終={h[-1]}")

# ---------- 留保2：隣接だけの輪に、型を見ない並べ替え（かき混ぜ）を入れる ----------
N = 60
def stir(x, k):
    # 位置 i の中身を位置 i+k+k+... へ：番地を k 飛びに読み直す（k と N は互いに素）
    y = [0] * N; j = 0
    for i in range(N):
        y[i] = x[j]
        j += k
        if j >= N: j += -N
    return y

def local(x, pair=True):
    y = list(x)
    for i in range(N):
        if x[i] == 0:
            s = x[i - 1] + x[(i + 1) % N]
            if s > 0: y[i] = 1
            elif s < 0: y[i] = -1
    if pair:
        z = list(y)
        for i in range(N):
            if y[i] != 0 and (y[i - 1] == -y[i] or y[(i + 1) % N] == -y[i]):
                z[i] = 0
        y = z
    return y

def seeds(extra):
    x = [0] * N
    for k in range(6):
        x[10 * k] = 1 if k % 2 == 0 else -1
    if extra: x[5] = extra
    return x

def cnt(x):
    p = sum(1 for v in x if v > 0); m = sum(1 for v in x if v < 0)
    return p, m

def ring(extra, k, pair=True, copy=True, steps=400):
    x = seeds(extra); p, m = cnt(x); D0 = p + -m
    hist = []
    for _ in range(steps):
        if copy: x = local(x, pair)
        if k: x = stir(x, k)
        p, m = cnt(x); hist.append((p, m, p + -m))
    return D0, hist

print("\n=== 留保2：かき混ぜを入れた輪 N=60 ===")
for k in (0, 7, 11, 13, 17):
    for extra, tag in ((0, "対称"), (1, "+1"), (-1, "-1")):
        D0, h = ring(extra, k)
        f1, f2 = judge(D0, h)
        print(f"k={k:2d} {tag:4s} F1={'OK' if f1 else 'NG'} F2={'OK' if f2 else 'NG'} 最終={h[-1]}")
# 負の対照：かき混ぜだけ（写さない・除かない）
for extra, tag in ((1, "+1"),):
    D0, h = ring(extra, 7, copy=False)
    f1, f2 = judge(D0, h)
    print(f"かき混ぜのみ k=7 {tag} F1={'OK' if f1 else 'NG'} F2={'OK' if f2 else 'NG'} 最終={h[-1]}")
# 対照：かき混ぜ＋写すだけ（組を除かない）
for extra, tag in ((1, "+1"), (-1, "-1")):
    D0, h = ring(extra, 7, pair=False)
    f1, f2 = judge(D0, h)
    print(f"かき混ぜ＋写すのみ k=7 {tag} F1={'OK' if f1 else 'NG'} F2={'OK' if f2 else 'NG'} 最終={h[-1]}")
