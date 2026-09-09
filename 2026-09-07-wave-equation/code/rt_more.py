#  rt_more.py ── rt_walk.py の続き。因数の中身・檻の素因数・壁の置き場所
exec(open('rt_walk.py').read().split('# ================================================================ 実行')[0])
verts, n_ico, edges, elen2 = rt_graph()
faces, vadj = rt_faces(verts, n_ico, edges)
fidx = {}
fedges = {}
for i, f in enumerate(faces):
    for k in range(4):
        fedges.setdefault(frozenset((f[k], f[(k+1)%4])), []).append(i)
nbr = {i: [] for i in range(len(faces))}
for e, fs in fedges.items():
    if len(fs) == 2:
        nbr[fs[0]].append(fs[1]); nbr[fs[1]].append(fs[0])

def to_arc_factor(g, t):
    d = len(g) - 1
    out = [0] * (2 * d + 1)
    for j, c in enumerate(g):
        if not c: continue
        term = [1]
        for _ in range(j): term = poly_mul(term, [t * t, 0, 1])
        for _ in range(d - j): term = poly_mul(term, [0, 1])
        for i, v in enumerate(term): out[i] += c * v
    return poly_trim(out)

A = [[0]*30 for _ in range(30)]
for i in nbr:
    for j in nbr[i]: A[i][j] = 1
cpA = charpoly(A)
facsA, restA = factor_small(cpA)
from collections import Counter
cntA = Counter(tuple(f) for f in facsA)
t = 2
print("面の網の固有値ごとの、辺の側の因数（x^2 - nu x + 4）")
arcfacs = Counter()
for f, k in sorted(cntA.items(), key=lambda kv:(len(kv[0]), kv[0])):
    af = to_arc_factor(list(f), t)
    sub, rest = factor_small(af)
    if len(rest) > 1: sub = sub + [rest]
    print(f"  nu の因数 {pstr(list(f)):<12} 重複{k:>2}  →  {pstr(af):<26} " +
          ("既約" if len(sub)==1 else "可約: " + " × ".join(pstr(x) for x in sub)))
    for x in sub: arcfacs[tuple(x)] += k
# x^2 - 4 の分（辺60 - 面30 = 30重）
for x in ([-2,1],[2,1]): arcfacs[tuple(x)] += 30
print()
print("相異なる因数と次数の和（クリロフ次元の上限）")
tot = 0
for f, k in sorted(arcfacs.items(), key=lambda kv:(len(kv[0]), kv[0])):
    print(f"  {pstr(list(f)):<26} 重複 {k}")
    tot += len(f)-1
print(f"  次数の和 = {tot}   （有向辺 120 に対して）")

print()
print("檻の素因数分解")
pred = [1]
for f, k in cntA.items():
    af = to_arc_factor(list(f), t)
    for _ in range(k): pred = poly_mul(pred, af)
for _ in range(30): pred = poly_mul(pred, [-4,0,1])
pred = poly_trim(pred)
for m in (5,7,11,13,17,19,23):
    k = cage(pred, m)
    fac = factorize(k)
    ss = " · ".join(f"{p}^{e}" if e>1 else f"{p}" for p,e in sorted(fac.items()))
    print(f"  m={m:<3} {k}  =  {ss}")

print()
print("壁を置ける場所（面の網の部分担体の次数）")
def direction(u, v):
    d = vsub(verts[u], verts[v])
    for comp in d:
        s = psign(comp)
        if s:
            if s < 0: d = tuple(pneg(c) for c in d)
            break
    return d
dirs = {}
for a,b in edges: dirs.setdefault(direction(a,b), []).append(frozenset((a,b)))
zones = []
for dvec, es in dirs.items():
    fs = set()
    for e in es: fs.update(fedges[e])
    zones.append(sorted(fs))
z = set(zones[0])
rest = [f for f in range(30) if f not in z]
# 帯の外の20面が2つの笠に割れるか
seen=set(); caps=[]
for f in rest:
    if f in seen: continue
    comp=[f]; seen.add(f); st=[f]
    while st:
        x=st.pop()
        for y in nbr[x]:
            if y in rest and y not in seen: seen.add(y); comp.append(y); st.append(y)
    caps.append(sorted(comp))
print(f"  帯10面を外すと残りは {[len(c) for c in caps]} 面に割れる")
for name, sub in [("帯 10面", sorted(z)), ("笠 %d面"%len(caps[0]), caps[0]),
                  ("帯+笠 %d面"%(10+len(caps[0])), sorted(z|set(caps[0])))]:
    S=set(sub); dd=Counter(len([g for g in nbr[f] if g in S]) for f in sub)
    ne=sum(dd[k]*k for k in dd)//2
    print(f"  {name:<12} 次数分布 {dict(sorted(dd.items()))}  辺 {ne}")
