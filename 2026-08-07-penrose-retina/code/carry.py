#  検定CY1  記憶は跳びの大きさによってどれだけ引き継げるか
#      設計：視線を番地 G から G+d へ跳ばしたとき、層の記憶（セルの番地の集合）に d を
#            足して引き継ぐ。ペンローズは平行移動で自分に重ならないので、ずらした先が
#            担体のセルであるとは限らない。その割合＝引き継ぎ率を、跳びの実長と残渣で見る。
#      測る量：引き継ぎ率（内側3層のセルのうち、d だけずらした先もセルである割合）
#      OK なら：引き継ぎ率が残渣 σ(d) で決まり、実長の長い跳びほど高い
#              ── 「跳ぶなら残渣の小さい番地へ」という指針が構造から出る
#      NG なら：実長とも残渣とも関係しない ── 引き継ぎは使えない
#      必ず落ちる設定：残渣が窓（正五角形）の差集合の外に出る d では引き継ぎ率が 0 に近いこと
exec(open('wind_core.py').read())
import math, itertools

XYC = {q: U.xy(q) for q in CELLS}
CS = set(CELLS)
O = U.xy(CTR)
GAPL = 1.618034

# 残渣（共役側）。b13_chain_units の xy は実空間。共役は φ → 1-φ の入れ替えで作る
def sigma(q):
    # (a,b,c,d) の Z[ζ10] 表現に対し、ζ → ζ^3 のガロア共役
    z = complex(0, 0)
    for k in range(4):
        z += q[k] * complex(math.cos(2 * math.pi * 3 * k / 10), math.sin(2 * math.pi * 3 * k / 10))
    return abs(z)

def real_len(q):
    x, y = U.xy(q)
    return math.hypot(x, y)

# 内側3層のセル（原点＝担体の中心）
INNER = []
for n in (1, 2, 3):
    lo, hi = n * GAPL, (n + 1) * GAPL
    INNER += [q for q in CELLS if lo <= math.hypot(XYC[q][0] - O[0], XYC[q][1] - O[1]) < hi]
print("内側3層のセル %d 枚\n" % len(INNER))

# 跳びの候補：小さい係数の格子ベクトル
cands = []
for a in range(-3, 4):
    for b in range(-3, 4):
        for c in range(-3, 4):
            for e in range(-3, 4):
                d = (a, b, c, e)
                if d == (0, 0, 0, 0):
                    continue
                L = real_len(d)
                if 0.5 <= L <= 25:
                    cands.append((L, sigma(d), d))
cands.sort()

rows = []
for L, S, d in cands:
    keep = sum(1 for q in INNER if U.zadd(q, d) in CS)
    rows.append((L, S, keep / len(INNER), d))

print("  跳びの実長  残渣σ(d)   実長×残渣  引き継ぎ率")
band = [(0.5, 1.5), (1.5, 2.5), (2.5, 3.5), (3.5, 4.5), (4.5, 6), (6, 8), (8, 12), (12, 25)]
for lo, hi in band:
    sel = [r for r in rows if lo <= r[0] < hi]
    if not sel:
        continue
    best = max(sel, key=lambda r: r[2])
    print("  %5.2f〜%-5.2f  最良 %.3f   %5.2f     %.3f  (%s)"
          % (lo, hi, best[1], best[0] * best[1], best[2], best[3]))

print("\n  ── 個別に見る（引き継ぎ率の高い順・上位10） ──")
for L, S, r, d in sorted(rows, key=lambda r: -r[2])[:10]:
    print("   実長 %6.3f  残渣 %.4f  積 %6.3f  引き継ぎ %.3f  d=%s" % (L, S, L * S, r, d))
print("\n  ── 引き継ぎ率の低い順（下位5） ──")
for L, S, r, d in sorted(rows, key=lambda r: r[2])[:5]:
    print("   実長 %6.3f  残渣 %.4f  積 %6.3f  引き継ぎ %.3f  d=%s" % (L, S, L * S, r, d))
