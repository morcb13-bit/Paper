# 二層のシェルピンスキー（星型正八面体）の中の六角形：立方体の対角線 4 本のそれぞれから見て正六角形になる輪。整数だけ
from itertools import product, combinations, permutations
dot=lambda a,b:sum(x*y for x,y in zip(a,b)); sub=lambda a,b:tuple(x-y for x,y in zip(a,b)); neg=lambda a:tuple(-x for x in a)
CORN=[tuple(3*x for x in s) for s in product((1,-1),repeat=3)]
edge=lambda p,q:dot(sub(p,q),sub(p,q))==36
AX=[(1,1,1),(-1,-1,1),(-1,1,-1),(1,-1,-1)]                  # 対角線 4 本（向き）
HEX={}
for ax in AX:
    ends=[c for c in CORN if dot(c,ax) in (9,-9)]
    rest=[c for c in CORN if c not in ends]
    # 輪の順に並べる
    ring=[rest[0]]
    while len(ring)<6:
        ring.append([c for c in rest if c not in ring and edge(c,ring[-1])][0])
    hs=[dot(c,ax) for c in ring]
    w=[tuple(3*dot(c,c)*0+3*x-dot(c,ax)*y for x,y in zip(c,ax)) for c in ring]   # 軸に垂直な成分×3
    ok=len({dot(v,v) for v in w})==1 and all(2*dot(w[i],w[(i+1)%6])==dot(w[i],w[i]) for i in range(6))
    HEX[ax]=ring
    print("軸",ax,"両端（ふた）",ends,"輪",ring,"軸の向きの高さ",hs,"軸から見て正六角形",ok)
E=[frozenset((p,q)) for p,q in combinations(CORN,2) if edge(p,q)]
cnt={e:sum(1 for r in HEX.values() if any(frozenset((r[i],r[(i+1)%6]))==e for i in range(6))) for e in E}
print("立方体の辺",len(E),"どの辺も入る六角形の数",set(cnt.values()))
print("どの角も入る六角形の数",{c:sum(1 for r in HEX.values() if c in r) for c in CORN}.values())
P,Q=(3,3,-3),(-3,-3,3)
thru=[ax for ax,r in HEX.items() if P in r and Q in r]
print("C9（P）と C10（Q）を両方含む六角形",len(thru),"その軸",thru,"含まない六角形の軸",[ax for ax in AX if ax not in thru])
for ax in thru:
    r=HEX[ax]; i,j=r.index(P),r.index(Q); print("  軸",ax,"で P と Q は輪の上で",(j-i)%6,"つ離れる（向かい合う＝3）")
# トリプチセンの芯の 3 本の道（m47 の 1 通り目）と、3 つの六角形の関係
PATHS=[[P,(3,3,3),(3,-3,3),Q],[P,(3,-3,-3),(-3,-3,-3),Q],[P,(-3,3,-3),(-3,3,3),Q]]
for k1,k2 in combinations(range(3),2):
    cyc=PATHS[k1]+PATHS[k2][::-1][1:-1]
    m=[ax for ax,r in HEX.items() if set(r)==set(cyc)]
    print("道",k1+1,"と",k2+1,"でできる 6 員環",cyc,"は六角形の軸",m)
print("---- トリプチセンの道は、どの六角形の半分か")
for k,p in enumerate(PATHS):
    es={frozenset((p[i],p[i+1])) for i in range(3)}
    m=[ax for ax,r in HEX.items() if es<={frozenset((r[i],r[(i+1)%6])) for i in range(6)}]
    print("道",k+1,p,"→ 六角形の軸",m)
