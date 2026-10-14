# m15: 正四面体を 頂点・辺・面 で共有させる（単・二重・三重結合の模型）
# 座標は立方体の一つおきの頂点を3倍したもの。整数と足し算・掛け算だけ。
import itertools

A = [(3, 3, 3), (3, -3, -3), (-3, 3, -3), (-3, -3, 3)]
O = (0, 0, 0)

def add(a, b): return tuple(x + y for x, y in zip(a, b))
def sub(a, b): return tuple(x - y for x, y in zip(a, b))
def dot(a, b): return sum(x * y for x, y in zip(a, b))
def sc(k, a): return tuple(k * x for x in a)

def point_reflect(p, m2):   # 中心 m の点対称（m2 = 2m を渡す）
    return sub(m2, p)

def plane_reflect(p, n, d):  # 平面 n·x = d での鏡映。n·n が割り切る場合だけ使う
    t = 2 * (dot(n, p) - d)
    nn = dot(n, n)
    assert t % nn == 0
    return sub(p, sc(t // nn, n))

def pair(c, p, q):  # 中心 c から見た p, q の組（内積, 長さ²）。長さ²が等しいときだけ角の組になる
    u, v = sub(p, c), sub(q, c)
    return dot(u, v), dot(u, u), dot(v, v)

def report(name, B, cB):
    shared = [p for p in A if p in B]
    freeA = [p for p in A if p not in B]
    freeB = [p for p in B if p not in A]
    print(f"== {name}: 共有 {len(shared)} 点  中心どうしの長さ² {dot(sub(cB,O),sub(cB,O))}")
    print("   相手の中心", cB, " A の空き", freeA, " B の空き", freeB)
    for p in freeA:
        print("   空き–C–C の組", pair(O, p, cB))
    for p, q in itertools.combinations(freeA, 2):
        print("   空き–C–空き の組", pair(O, p, q))
    pts = freeA + freeB
    if len(pts) >= 4:
        base = pts[0]
        vs = [sub(p, base) for p in pts[1:]]
        cr = (vs[0][1]*vs[1][2]-vs[0][2]*vs[1][1], vs[0][2]*vs[1][0]-vs[0][0]*vs[1][2], vs[0][0]*vs[1][1]-vs[0][1]*vs[1][0])
        print("   空き4点が同一平面か:", all(dot(cr, v) == 0 for v in vs), " 法線", cr)
    if len(pts) == 2:
        d1, d2 = sub(pts[0], O), sub(cB, O)
        cr = (d1[1]*d2[2]-d1[2]*d2[1], d1[2]*d2[0]-d1[0]*d2[2], d1[0]*d2[1]-d1[1]*d2[0])
        e = sub(pts[1], cB)
        cr2 = (e[1]*d2[2]-e[2]*d2[1], e[2]*d2[0]-e[0]*d2[2], e[0]*d2[1]-e[1]*d2[0])
        print("   空き2点が軸の上か:", cr == (0, 0, 0) and cr2 == (0, 0, 0))
    return B, cB

# 頂点共有（単結合）：共有点 v で点対称
v = A[0]
report("頂点共有（単結合 C–C）", [point_reflect(p, sc(2, v)) for p in A], sc(2, v))

# 辺共有（二重結合 C=C）：辺の中点で点対称
m2 = add(A[0], A[1])           # 2m = (6,0,0)
report("辺共有（二重結合 C=C）", [point_reflect(p, m2) for p in A], m2)

# 面共有（三重結合 C≡C）：面の平面で鏡映。面 A0,A1,A2 は n=(1,1,-1), n·x=3
n, d = (1, 1, -1), 3
assert all(dot(n, p) == d for p in A[:3])
report("面共有（三重結合 C≡C）", [plane_reflect(p, n, d) for p in A], plane_reflect(O, n, d))

# CO2：中心 C の向かい合う二辺をそれぞれ別の正四面体と共有
print("== CO2：向かい合う二辺 (A0A1) と (A2A3)")
mO1 = add(A[0], A[1]); mO2 = add(A[2], A[3])
print("   O の中心", mO1, mO2, " 和", add(mO1, mO2), "（0 なら O=C=O が一直線）")
e1, e2 = sub(A[0], A[1]), sub(A[2], A[3])
print("   二つの辺の向き", e1, e2, " 内積", dot(e1, e2), "（0 なら二つの面が直交）")
print("   長さ²の階段 頂点:辺:面 =", dot(sc(2, v), sc(2, v)), ":", dot(m2, m2), ":", dot(plane_reflect(O, n, d), plane_reflect(O, n, d)))
