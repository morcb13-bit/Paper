# 四面体（正四面体ではない）での右型・左型（整数だけ）
#  対象：中心を原点に置き、置換基 A,B,C,D を長さも向きもばらばらの整数の点に置いた四面体。
#  右か左かは、向きの符号 s = det[B-A, C-A, D-A] の正負で決める（整数）。
#  検定T0 4つの腕の長さの二乗がすべて異なり、中心が四面体の内側にある（正四面体ではない）
#  検定T1 頂点を立方体の角に移す整数の回転24通りで回しても、s の符号は変わらない
#  検定T2 負の対照：鏡映を含む24通り（行列式 −1）では、s の符号がすべて反転する
#  検定T3 鏡（x を −x に）に映した左型は s の符号が逆
#  検定T4 型紙＝受け口 B,C,D と「板の上側」。はまる ⇔ 受け口の三角形と同じ形に置けて、かつ A が上側
#         A が上側か ＝ det[C-B, D-B, A-B] の符号。右型は上、左型は（B,C,D を合わせると）下
#  検定T5 負の対照：板の上下を問わない受け口（宙に浮いた3点）なら、左型もはまる
from itertools import permutations, product
from fractions import Fraction
R = {"A": (0, 0, 4), "B": (3, 0, -1), "C": (-1, 3, -2), "D": (-2, -2, -1)}
def sub(a, b): return tuple(a[i] - b[i] for i in range(3))
def det3(u, v, w):
    return (u[0]*(v[1]*w[2] - v[2]*w[1]) - u[1]*(v[0]*w[2] - v[2]*w[0]) + u[2]*(v[0]*w[1] - v[1]*w[0]))
def s(m): return det3(sub(m["B"], m["A"]), sub(m["C"], m["A"]), sub(m["D"], m["A"]))
def sg(x): return (x > 0) - (x < 0)
L2 = {k: sum(c * c for c in v) for k, v in R.items()}
# 中心が内側か：原点と各頂点が、反対の面について同じ側にあるか
inside = True
for k in "ABCD":
    o = [j for j in "ABCD" if j != k]
    a, b, c = (R[j] for j in o)
    n1 = det3(sub(b, a), sub(c, a), sub(R[k], a)); n0 = det3(sub(b, a), sub(c, a), sub((0, 0, 0), a))
    if sg(n1) != sg(n0) or n0 == 0: inside = False
print("T0 腕の長さの二乗", L2, "すべて異なる", len(set(L2.values())) == 4, "中心が内側", inside)
mats = []
for pm in permutations(range(3)):
    for sgn in product((1, -1), repeat=3):
        M = [[0] * 3 for _ in range(3)]
        for i in range(3): M[i][pm[i]] = sgn[i]
        mats.append(M)
def ap(M, v): return tuple(sum(M[i][j] * v[j] for j in range(3)) for i in range(3))
def detM(M): return det3(*[tuple(M[i][j] for i in range(3)) for j in range(3)])
s0 = s(R)
rot = [M for M in mats if detM(M) == 1]; ref = [M for M in mats if detM(M) == -1]
print("右型の s =", s0)
print("T1 回転", len(rot), "通りで符号が保たれる", all(sg(s({k: ap(M, v) for k, v in R.items()})) == sg(s0) for M in rot))
print("T2 鏡映を含む", len(ref), "通りで符号が反転", all(sg(s({k: ap(M, v) for k, v in R.items()})) == -sg(s0) for M in ref))
Lm = {k: (-v[0], v[1], v[2]) for k, v in R.items()}
print("T3 左型の s =", s(Lm), "逆符号", sg(s(Lm)) == -sg(s0))
up = lambda m: sg(det3(sub(m["C"], m["B"]), sub(m["D"], m["B"]), sub(m["A"], m["B"])))
# 受け口の三角形 B,C,D の形（辺の長さの二乗）は右と左で同じ
tri = lambda m: tuple(sum(c * c for c in sub(m[x], m[y])) for x, y in (("B", "C"), ("C", "D"), ("D", "B")))
print("   受け口の三角形（辺の二乗）右", tri(R), "左", tri(Lm), "同じ形", tri(R) == tri(Lm))
print("T4 板の上側に A が来るか：右", up(R) == up(R), " 左", up(Lm) == up(R),
      "→ 右だけはまる" if up(R) != up(Lm) else "→ NG")
print("T5 上下を問わない受け口：右 はまる, 左 はまる（三角形が同じ形なので）", tri(R) == tri(Lm))
