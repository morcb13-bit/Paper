# 偏りの増幅検定（整数と加算・比較だけ）
# 判定基準は走らせる前にここで固定する
#   F1: ある T で |D(T)| >= 2|D(0)| となり、以後それを下回らない
#   F2: 最終状態で 20*N_same >= 19*N_occ（占有のうち95%以上が同じ側）が保たれる
#   対照: 完全対称の初期状態では D=0 が続く（実装確認）
#         足し込みのない規則では F1 が NG になる（負の対照）
#   鏡像: +を1多くした場合と -を1多くした場合で D の符号が反転し、絶対値が一致

# ---------- (1) よく混ぜた器：個数だけを持つ ----------
def mixed(R, L, food, rule, steps=60):
    hist = []
    for _ in range(steps):
        if rule in ("self", "self+pair"):
            # 自己足し込み：各側が自分の数を足す（食料が尽きるまで）
            need = R + L
            if need <= food:
                food += -need; R, L = R + R, L + L
            else:
                # 食料を比率どおりに割らず、1個ずつ交互に配る（割り算を使わない）
                gr = gl = 0
                while food > 0 and (gr < R or gl < L):
                    if gr < R and food > 0: gr += 1; food += -1
                    if gl < L and food > 0: gl += 1; food += -1
                R, L = R + gr, L + gl
        if rule == "self+pair":
            # 右型と左型の組をつくって除く（組んだものは働かない）
            m = R if R < L else L
            R, L = R + -m, L + -m
        if rule == "none":
            pass  # 足し込みなし：何も増えない
        hist.append((R, L, R + -L))
    return hist

# ---------- (2) 隣接だけを見る輪 ----------
def ring_step(x, rule):
    n = len(x); y = list(x)
    for i in range(n):
        a, b = x[i - 1], x[(i + 1) % n]
        if rule in ("copy", "copy+pair") and x[i] == 0:
            s = a + b
            if s > 0: y[i] = 1
            elif s < 0: y[i] = -1
        if rule == "diffuse":
            # 足し込みなしの負の対照：左隣の値を受け取るだけ（総和保存）
            y[i] = a
    if rule == "copy+pair":
        z = list(y)
        for i in range(n):
            if y[i] != 0 and (y[i - 1] == -y[i] or y[(i + 1) % n] == -y[i]):
                z[i] = 0
        y = z
    return y

def counts(x):
    p = sum(1 for v in x if v > 0); m = sum(1 for v in x if v < 0)
    return p, m

def ring_run(x, rule, steps=200):
    hist = []
    for _ in range(steps):
        x = ring_step(x, rule)
        p, m = counts(x); hist.append((p, m, p + -m))
    return hist, x

def judge(D0, hist):
    Ds = [h[2] for h in hist]
    a0 = D0 if D0 >= 0 else -D0
    f1 = False
    for t in range(len(Ds)):
        if all((d if d >= 0 else -d) >= a0 + a0 for d in Ds[t:]) and a0 > 0:
            f1 = True; break
    R, L, D = hist[-1]
    occ = R + L; same = R if D >= 0 else L
    f2 = occ > 0 and 20 * same >= 19 * occ and D != 0
    return f1, f2, hist[-1]

print("=== (1) よく混ぜた器 ===")
for rule in ("none", "self", "self+pair"):
    for (R, L, tag) in ((50, 50, "対称"), (51, 50, "+1"), (50, 51, "-1")):
        h = mixed(R, L, 100000, rule)
        f1, f2, last = judge(R + -L, h)
        print(f"{rule:10s} {tag:4s} F1={'OK' if f1 else 'NG'} F2={'OK' if f2 else 'NG'} 最終(R,L,D)={last}")

print("\n=== (2) 隣接だけの輪 N=60 ===")
N = 60
def seeds(extra):
    # 10おきに +,- を交互に置く（左右反転と符号反転で自分に戻る配置）
    x = [0] * N
    for k in range(6):
        x[10 * k] = 1 if k % 2 == 0 else -1
    if extra:
        x[5] = extra  # 隙間の中央に1個足す
    return x
for rule in ("diffuse", "copy", "copy+pair"):
    for extra, tag in ((0, "対称"), (1, "+1"), (-1, "-1")):
        x = seeds(extra)
        p, m = counts(x)
        h, xf = ring_run(x, rule)
        f1, f2, last = judge(p + -m, h)
        print(f"{rule:10s} {tag:4s} F1={'OK' if f1 else 'NG'} F2={'OK' if f2 else 'NG'} 最終(N+,N-,D)={last}")
