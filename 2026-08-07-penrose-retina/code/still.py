#  ペンタゴン網膜：静止した輪郭を、粗い段から細かい段へ降りて読む
#  監督の手順：全体を眺め、細部にズーミングしていく。
#
#  対象を一行で（判別法9）：閉じた輪郭を網膜に置いたとき、輪郭上の各点が
#  どの記号になるか。段を上げる（＝形を φ² 倍する＝ズームイン）と何が変わるか。
#
#  検定S1 静止した輪郭は記号の列になるか
#      OK なら：動く点で出た6記号が、輪郭にもそのまま出る
#      NG なら：静止画には別の規則が要る
#
#  検定S2 段を上げると記号が増えるか
#      OK なら：ズームインで細部が出る
#      NG なら：段は静止画に効かない
#
#  検定S3 粗い段で励起した点は、細かい段でも励起するか
#      OK なら：降りても同じ場所を見ている＝全体から細部へ辿れる
#      NG なら：段ごとに別の場所を見ている＝降りる意味が無い
#
#  輪郭は Z[ζ₁₀] の折れ線を有理数の刻みで標本するので、全点が体の中にある。
#  判定はすべて Q(φ) の符号。浮動小数は候補の絞り込みと表示だけ。

import json, math
from fractions import Fraction as F
from qphi import Qp, zmul, zconj, zre, ZERO

Z = (0, 1, 0, 0); PHI = (1, 0, 1, -1); ONE = (1, 0, 0, 0)
def zadd(a, b): return tuple(a[i] + b[i] for i in range(4))
def zsub(a, b): return tuple(a[i] - b[i] for i in range(4))
def zneg(a): return tuple(-x for x in a)
def zsc(a, r): return tuple(x * r for x in a)
def zrot(a, k):
    b = a
    for _ in range(k % 5): b = zmul(b, Z)
    return b if (k % 10) < 5 else zneg(b)
def zim(a): return Qp(a[1], a[2] + a[3])
def norm2(a): return zre(zmul(a, zconj(a)))
ZZ = math.cos(math.pi/5) + 1j*math.sin(math.pi/5)
def num(a): return sum(float(a[k]) * ZZ**k for k in range(4))

d = json.load(open("rings_integer.json"))
rings = [tuple(F(x) for x in v) for v in d["rings"]]
pents = [(tuple(F(x) for x in key.split(",")), a) for key, a in d["cells"].items()]
stars = [tuple(F(x) for x in v) for v in d["pentagrams"]]

# 2つ飛ばし対（|v|² = 3+4φ）とその中点
FLY = Qp(3, 4)
pairs = [(a, b) for i, a in enumerate(rings) for b in rings[i+1:]
         if norm2(zsub(a, b)) == FLY]
mids = [zsc(zadd(a, b), F(1, 2)) for a, b in pairs]
paired = set()
for a, b in pairs: paired.add(a); paired.add(b)
print(f"担体：環{len(rings)} 五角形{len(pents)} 五芒星{len(stars)} "
      f"2つ飛ばし対{len(pairs)}")

# --- 五角形の頂点（厳密） ---
VERT = {}
for g, a in pents:
    VERT[(g, a)] = [zadd(g, zrot(ONE, a + 2*i)) for i in range(5)]
PC = [(num(g), g, a) for g, a in pents]

def inside(p, g, a):
    vs = VERT[(g, a)]
    s = 0
    for i in range(5):
        e = zsub(vs[(i+1) % 5], vs[i])
        w = zsub(p, vs[i])
        c = zim(zmul(zconj(e), w)).sign()
        if c == 0: continue
        if s == 0: s = c
        elif s != c: return False
    return True

MARKS = ([(num(v), "星") for v in stars]
         + [(num(v), "舟" if v in paired else "十") for v in rings]
         + [(num(v), "菱") for v in mids])

def symbol(p):
    pn = num(p)
    for cn, g, a in PC:
        if abs(pn - cn) < 1.05 and inside(p, g, a):
            return "P0" if a % 2 == 0 else "P1"
    return min(MARKS, key=lambda m: abs(pn - m[0]))[1]

# --- 輪郭：Z[ζ₁₀] の三角形。有理数の刻みで標本 ---
TRI = [(2, 0, 0, 0), (0, 2, 1, 0), (-1, -1, 1, 0)]
M = 40
def contour(m):
    """段 m：形を φ^{2m} 倍する（＝ズームイン）"""
    s = ONE
    for _ in range(m): s = zmul(zmul(PHI, PHI), s)
    V = [zmul(s, tuple(F(x) for x in v)) for v in TRI]
    pts = []
    for i in range(3):
        A, B = V[i], V[(i+1) % 3]
        for j in range(M):
            t = F(j, M)
            pts.append(tuple(A[k] + (B[k]-A[k])*t for k in range(4)))
    return pts

print()
print("=" * 66)
res = {}
for m in range(3):
    pts = contour(m)
    R = max(norm2(p).val()**0.5 for p in pts)
    syms = [symbol(p) for p in pts]
    res[m] = syms
    kinds = sorted(set(syms))
    exc = sum(1 for s in syms if s.startswith("P"))
    runs = 1 + sum(1 for i in range(1, len(syms)) if syms[i] != syms[i-1])
    print(f"段 φ^{2*m}  外接半径 {R:6.3f}  記号 {len(kinds)}種 {kinds}")
    print(f"          励起 {exc:3d}/{len(syms)}   記号の切替 {runs:3d} 回")
    print(f"          先頭24点: {' '.join(syms[:24])}")

print()
allk = set().union(*[set(v) for v in res.values()])
print("検定S1:", "OK  6記号の中に収まる" if allk <= {"P0","P1","星","十","舟","菱"} else "NG")
print("検定S2:", "OK  段を上げると切替が増える"
      if sum(1 for i in range(1,len(res[2])) if res[2][i]!=res[2][i-1])
         > sum(1 for i in range(1,len(res[0])) if res[0][i]!=res[0][i-1]) else "NG")

keep = [(a, b) for a, b in zip(res[0], res[1])]
n0 = sum(1 for a, _ in keep if a.startswith("P"))
n01 = sum(1 for a, b in keep if a.startswith("P") and b.startswith("P"))
keep2 = [(a, b) for a, b in zip(res[1], res[2])]
n1 = sum(1 for a, _ in keep2 if a.startswith("P"))
n12 = sum(1 for a, b in keep2 if a.startswith("P") and b.startswith("P"))
print(f"検定S3: 段0→1 で励起が残る率 {n01}/{n0}"
      f"   段1→2 で {n12}/{n1}")
print("        （同じ輪郭パラメータの点が、次の段でも励起するか）")
