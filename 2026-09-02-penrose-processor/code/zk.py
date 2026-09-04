import json, math, pickle
from collections import defaultdict, Counter
import b13_chain_units_big as U
phi = (1 + 5 ** .5) / 2
SC = pickle.load(open("rscan_0.pkl", "rb"))

# Z[φ] の道具
# (a,b) は a + bφ
def val(x): return x[0] + x[1] * phi
def divphi(x): return (x[1] - x[0], x[0])          # ×φ⁻¹
def mulphi(x): return (x[1], x[0] + x[1])          # ×φ
def reduce_phi(x):
    if x == (0, 0): return 0, (0, 0)
    k = 0
    while val(x) >= phi * phi - 1e-9: x = divphi(x); k += 1
    while val(x) < 1 - 1e-9: x = mulphi(x); k -= 1
    return k, x

FIB = [1, 2]
while FIB[-1] < 10 ** 6: FIB.append(FIB[-1] + FIB[-2])
def zeck(n):
    out = []
    for f in reversed(FIB):
        if f <= n: out.append(f); n -= f
    return out
FSET = set()
a, b = 1, 1
while b < 10 ** 6: FSET.add(b); a, b = b, a + b
LSET = set()
a, b = 2, 1
while b < 10 ** 6: LSET.add(b); a, b = b, a + b

rows = []
for c, R, lim in SC:
    r = math.hypot(*U.xy(c)) / 10.0
    p, q = U.norm2(c)
    if p % 100 or q % 100: continue
    p, q = p // 100, q // 100
    k, x0 = reduce_phi((p, q))
    rows.append(dict(R=round(R, 4), r=r, p=p, q=q, N=p*p + p*q - q*q, k=k, x0=x0,
                     zp=len(zeck(p)) if p else 0, zq=len(zeck(q)),
                     fib=(p in FSET and q in FSET), luc=(p in LSET and q in LSET)))
print("対象の星 %d 個" % len(rows))
print()
byR = defaultdict(list)
for x in rows: byR[x["R"]].append(x)
for R in sorted(byR):
    g = byR[R]
    print("R = %8.4f   %4d 個" % (R, len(g)))
    print("   k の分布      %s" % dict(sorted(Counter(x["k"] for x in g).items())))
    x0s = Counter(x["x0"] for x in g)
    print("   x₀ の種類 %d   上位 %s" % (len(x0s), x0s.most_common(4)))
    print("   N の種類 %d    上位 %s" % (len(set(x["N"] for x in g)),
                                          Counter(x["N"] for x in g).most_common(4)))
    print("   ゼッケンドルフ長 (p,q) の種類 %d  上位 %s"
          % (len(set((x["zp"], x["zq"]) for x in g)),
             Counter((x["zp"], x["zq"]) for x in g).most_common(3)))
    print("   フィボナッチ対 %d   リュカ対 %d" % (sum(x["fib"] for x in g), sum(x["luc"] for x in g)))
    print()

print("判定：R の類の中で一定になった量")
for name, key in (("k", lambda x: x["k"]), ("x₀", lambda x: x["x0"]),
                  ("N", lambda x: x["N"]),
                  ("ゼッケンドルフ長", lambda x: (x["zp"], x["zq"]))):
    ok = all(len(set(key(x) for x in byR[R])) == 1 for R in byR)
    dif = len(set(key(byR[R][0]) for R in byR)) == len(byR)
    print("  %-16s 類の中で一定 %-5s   類のあいだで違う %s" % (name, ok, dif))
