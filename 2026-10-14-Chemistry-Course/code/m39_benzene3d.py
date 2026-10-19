# ベンゼン環を立体で：各炭素の三角形（結合 3 本の先）に 4 本目の頂点を戻して正四面体にする。整数だけ
from itertools import combinations
dot=lambda a,b:sum(x*y for x,y in zip(a,b))
add=lambda a,b:tuple(x+y for x,y in zip(a,b))
sub=lambda a,b:tuple(x-y for x,y in zip(a,b))
neg=lambda a:tuple(-x for x in a)
d2=lambda a,b:dot(sub(a,b),sub(a,b))
GA=[(4,-2,-2),(-2,4,-2),(-2,-2,4)]; GB=[neg(d) for d in GA]
C=[(4,-2,-2),(2,2,-4),(-2,4,-2),(-4,2,2),(-2,-2,4),(2,-4,2)]   # m38 と同じ（A,B,A,B,A,B）
isA=lambda i:i%2==0
H=[]
for i,p in enumerate(C):
    for d in (GA if isA(i) else GB):
        q=add(p,d)
        if q not in C: H.append(q)
# 4 本目の頂点：三角形の重心（＝炭素）から (1,1,1) の向きへ、辺² 72 になる高さ
for k in range(1,9):
    if all(d2(add(GA[0],(0,0,0)),(k,k,k))==72 for _ in [0]): print("4 本目の頂点 (k,k,k) の k",k)
TOP=(4,4,4)
TET=[]
for i,p in enumerate(C):
    base=[add(p,d) for d in (GA if isA(i) else GB)]
    apex=add(p,TOP) if isA(i) else sub(p,TOP)
    T=base+[apex]; TET.append(T)
    assert all(d2(a,b)==72 for a,b in combinations(T,2))
print("6 個とも正四面体（辺² 72）")
print("頂点の内訳：A の四面体 =",["B,B,H,上" for _ in range(3)][0],"／ B の四面体 = A,A,H,下")
# 大きな四面体（シェルピンスキー 2 段目）：A の 3 個＋上の 1 個
AP=[TET[i][3] for i in (0,2,4)]; BP=[TET[i][3] for i in (1,3,5)]
for name,P,cs,sg in [("A（上）",AP,[H[i] for i in (0,2,4)],1),("B（下）",BP,[H[i] for i in (1,3,5)],-1)]:
    apx=tuple(sum(x)//3 for x in zip(*P)); s=sum(apx[i]-sum(c[i] for c in cs)//3 for i in range(1))
    # 大きな頂点：角 3 つの重心から上へ。上の小四面体＝4 本目の頂点 3 つ＋大きな頂点
    for k in range(1,31):
        k=sg*k
        top=add(tuple(sum(x)//3 for x in zip(*cs)),(k,k,k))
        if all(d2(top,c)==288 for c in cs): break
    small_top=P+[top]
    print(name,"大きな三角形の角",cs,"角どうし辺²",[d2(a,b) for a,b in combinations(cs,2)])
    print("   大きな頂点",top,"角との辺²",[d2(top,c) for c in cs],"／ 4 本目の頂点 3 つの辺²",[d2(a,b) for a,b in combinations(P,2)])
    print("   4 本目の頂点は大きな辺の中点か",all(add(c,top)==tuple(2*x for x in m) for c,m in zip(cs,P)) or sorted(tuple(2*x for x in m) for m in P)==sorted(add(c,top) for c in cs))
    print("   上の小四面体（4 本目の頂点 3＋大きな頂点）辺²",sorted(set(d2(a,b) for a,b in combinations(small_top,2))))
print("B は A の反転（中心について）か", sorted(neg(x) for T in TET[0::2] for x in T)==sorted(x for T in TET[1::2] for x in T))
# 番地の格子：辺 6 の立方体の格子に乗るか
mod=lambda p:tuple(x%6 for x in p)
print("A の四面体の頂点 mod 6",set(mod(x) for T in TET[0::2] for x in T),"／ B",set(mod(x) for T in TET[1::2] for x in T))
print("A の炭素 mod 6",set(mod(C[i]) for i in (0,2,4)),"／ B の炭素",set(mod(C[i]) for i in (1,3,5)))
# 隣どうしの四面体の接し方
for i in range(6):
    j=(i+1)%6; sh=set(TET[i])&set(TET[j])
    print("四面体",i,"と",j,"共有する頂点",len(sh),"／ 炭素",C[i],"は相手の頂点か",C[i] in TET[j])
import json; json.dump({"C":C,"H":H,"TET":TET},open("benzene3d.json","w"))
