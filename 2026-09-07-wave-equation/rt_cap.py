#  rt_cap.py ── 帯が空だったので、壁のある部分担体（笠10面）で走らせる
#  検定C1 笠の担体
#      OK なら：10面、壁（次数2）が5枚・内側（次数4）が5枚、辺15
#  検定C2 一歩
#      OK なら：dc = 2t が全次数で書ける最小は dc=4（t=2）。T^T T = 4I
#      NG なら：壁のある担体では一歩が書けない
#  検定C3 井戸になっているか（1環との対比）
#      判別式に 5 が出る なら：phi が時間の側に乗る
exec(open('rt_walk.py').read().split('# ================================================================ 実行')[0])
from collections import Counter
verts, n_ico, edges, elen2 = rt_graph()
faces, vadj = rt_faces(verts, n_ico, edges)
fedges = {}
for i, f in enumerate(faces):
    for k in range(4): fedges.setdefault(frozenset((f[k], f[(k+1)%4])), []).append(i)
nbr = {i: [] for i in range(len(faces))}
for e, fs in fedges.items():
    if len(fs)==2: nbr[fs[0]].append(fs[1]); nbr[fs[1]].append(fs[0])
def direction(u,v):
    d = vsub(verts[u], verts[v])
    for c in d:
        s = psign(c)
        if s:
            if s<0: d = tuple(pneg(x) for x in d)
            break
    return d
dirs={}
for a,b in edges: dirs.setdefault(direction(a,b),[]).append(frozenset((a,b)))
zones=[]
for dv,es in dirs.items():
    fs=set()
    for e in es: fs.update(fedges[e])
    zones.append(sorted(fs))
z=set(zones[0]); restf=[f for f in range(30) if f not in z]
seen=set(); caps=[]
for f in restf:
    if f in seen: continue
    comp=[f]; seen.add(f); st=[f]
    while st:
        x=st.pop()
        for y in nbr[x]:
            if y in restf and y not in seen: seen.add(y); comp.append(y); st.append(y)
    caps.append(sorted(comp))
cap=caps[0]; S=set(cap)
nb={f:[g for g in nbr[f] if g in S] for f in cap}
dd=Counter(len(nb[f]) for f in cap)
print("検定C1 笠の担体")
print(f"  面 {len(cap)}  次数分布 {dict(sorted(dd.items()))}  辺 {sum(dd[k]*k for k in dd)//2}")
print(f"  → {'OK' if dict(dd)=={2:5,4:5} else 'NG'}")

def to_arc_factor(g,t):
    d=len(g)-1; out=[0]*(2*d+1)
    for j,c in enumerate(g):
        if not c: continue
        term=[1]
        for _ in range(j): term=poly_mul(term,[t*t,0,1])
        for _ in range(d-j): term=poly_mul(term,[0,1])
        for i,v in enumerate(term): out[i]+=c*v
    return poly_trim(out)

for t in (1,2):
    try:
        T,arcs,_=arc_operator(nb,t)
        n=len(arcs)
        TT=mat_mul([[T[j][i] for j in range(n)] for i in range(n)],T)
        ok=all(TT[i][j]==(t*t if i==j else 0) for i in range(n) for j in range(n))
        print(f"\n検定C2 一歩  t={t}  有向辺 {n}   T^T T = {t*t} I : {ok}")
    except AssertionError as e:
        print(f"\n検定C2 一歩  t={t} → 落ちる: {e}")

t=2; T,arcs,_=arc_operator(nb,t); n=len(arcs)
cp=charpoly(T)
fs,rest=factor_small(cp)
cnt=Counter(tuple(x) for x in fs)
print("\n笠の上の一歩の特性多項式（120ではなく %d次元）"%n)
tot=0
for f,k in sorted(cnt.items(), key=lambda kv:(len(kv[0]),kv[0])):
    print(f"  {pstr(list(f)):<26} 重複 {k}")
    tot+=len(f)-1
if len(rest)>1:
    print(f"  残り {pstr(rest)}（{len(rest)-1}次）")
    tot+=len(rest)-1
print(f"  相異なる因数の次数の和 = {tot}  （有向辺 {n} に対して）")

print("\n檻（T^k = I の最小 k）")
for m in (5,7,11,13,17,19,23):
    k=cage(cp,m)
    if k is None: print(f"  m={m:<3} 可逆でない"); continue
    fac=factorize(k); ss=" · ".join(f"{p}^{e}" if e>1 else f"{p}" for p,e in sorted(fac.items()))
    print(f"  m={m:<3} {k} = {ss}")

# 壁で跳ね返るか（1枚に置いた状態を追う）
print("\n壁の効き（次数2の面に立てた状態を1歩ずつ）")
wall=[f for f in cap if len(nb[f])==2][0]
idx={a:i for i,a in enumerate(arcs)}
v=[0]*n; v[idx[(wall, nb[wall][0])]]=1
for step in range(1,7):
    v=[sum(T[i][j]*v[j] for j in range(n)) for i in range(n)]
    occ=Counter()
    for (u,w),i in idx.items():
        if v[i]: occ[u]+=abs(v[i])
    print(f"  {step}歩  二乗和 {sum(x*x for x in v):>8}   居る面 {sorted(occ)}")
