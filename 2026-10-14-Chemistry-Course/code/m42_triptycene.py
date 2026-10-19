# トリプチセン：ベンゼン環 3 つ（星型正八面体の帯）を、橋の炭素 C9・C10（正四面体の中心）でつなぐ。整数だけ
# 決めごと：芳香環の C–C は辺を分け合う。橋の C–C（単結合）は、相手の正四面体の頂点に自分の中心が入る（ダイヤモンドと同じ、長さ² 27）
from itertools import product, combinations
add=lambda a,b:tuple(x+y for x,y in zip(a,b)); sub=lambda a,b:tuple(x-y for x,y in zip(a,b)); neg=lambda a:tuple(-x for x in a)
dot=lambda a,b:sum(x*y for x,y in zip(a,b)); d2=lambda a,b:dot(sub(a,b),sub(a,b))
half=lambda a:tuple(x//2 for x in a)
DA=[(3,3,3),(3,-3,-3),(-3,3,-3),(-3,-3,3)]          # A の正四面体（中心からの 4 本）
def tet(c,t): return [add(c,d) for d in (DA if t=="A" else [neg(d) for d in DA])]
def neighbor(T,tip,edge):
    """T の底の面（tip 以外）の辺 edge を挟んで、星型正八面体の隣の面に立つ正四面体"""
    base=[v for v in T if v!=tip]; r=[v for v in base if v not in edge][0]; p,q=edge
    s=add(add(p,q),r); o2=sub(s,tip)                      # 2o = p+q+r−t
    r2=sub(o2,r)                                           # r' = 2o − r
    s2=add(add(p,q),r2); t2=sub(s2,o2)                     # t' = 3g' − 2o = (p+q+r') − 2o
    T2=[p,q,r2,t2]; c=tuple(sum(v[i] for v in T2) for i in range(3))
    assert all(x%4==0 for x in c); return T2,t2,tuple(x//4 for x in c)
def ring_from(c0,T0,tip0,cap_edge):
    """炭素 c0 から、ふたと分け合う辺を cap_edge にして帯を一周する。戻り値：中心・正四面体・とがった先の列"""
    base=[v for v in T0 if v!=tip0]; es=[e for e in combinations(base,2) if set(e)!=set(cap_edge)]
    out=[(c0,T0,tip0)]; T,tip,came=T0,tip0,None
    e=es[0]
    for _ in range(5):
        T2,t2,c2=neighbor(T,tip,e); out.append((c2,T2,t2))
        # 次の辺：共有した辺 e でもなく、ふたへの辺でもない辺（ふたへの辺は A/B で上下が交互）
        b2=[v for v in T2 if v!=t2]; cand=[f for f in combinations(b2,2) if set(f)!=set(e)]
        # ふたとの辺は、帯で一つおきに同じふたへ向く。隣の正四面体を作って、既にある帯の正四面体と重ならない方を選ぶ
        nxt=None
        for f in cand:
            T3,_,c3=neighbor(T2,t2,f)
            if c3 not in [o[0] for o in out] or (len(out)==6 and c3==c0):
                # 帯の向き：正八面体の中心から見て、ふたでない（軸のまわりを回る）側
                nxt=f if nxt is None else nxt
        T,tip,e=T2,t2,nxt
    return out
C9=(0,0,0); H9=(3,3,3); T9=tet(C9,"A")
arms=[v for v in T9 if v!=H9]                              # 環の炭素 3 個の居場所
print("C9 の正四面体",T9,"環の炭素の居場所",arms)
sols=[]
choices=[]
for a in arms:
    Ta=tet(a,"B"); tip=C9                                   # a のとがった先が C9
    assert tip in Ta
    base=[v for v in Ta if v!=tip]; opts=[]
    for cap in combinations(base,2):
        for e in [f for f in combinations(base,2) if set(f)!=set(cap)]:
            T2,t2,c2=neighbor(Ta,tip,e); opts.append({"a":a,"cap":cap,"e":e,"b":c2,"Tb":T2,"tipb":t2})
    choices.append(opts)
for combo in product(*choices):
    C10s={o["tipb"] for o in combo}
    if len(C10s)==1:
        C10=C10s.pop(); bs=[o["b"] for o in combo]
        T10=tet(C10,"B") if all(b in tet(C10,"B") for b in bs) else (tet(C10,"A") if all(b in tet(C10,"A") for b in bs) else None)
        sols.append((C10,combo,T10))
print("3 つの環の b のとがった先が 1 点にそろう選び方",len(sols))
seen=set()
for C10,combo,T10 in sols:
    k=(C10,tuple(o["b"] for o in combo))
    if k in seen: continue
    seen.add(k)
    print("C10",C10,"b",[o["b"] for o in combo],"C10 の正四面体に b が 3 個とも入る",T10 is not None,"C9–C10 の長さ²",d2(C9,C10))
