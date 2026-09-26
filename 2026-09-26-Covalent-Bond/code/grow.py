# 一歩先：一つの正四面体のすべての辺で「互いの中心を通る」相手を置く（辺のまわりに ±120°）
import sympy as sp, itertools, numpy as np
T=[sp.Matrix(v) for v in [(1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1)]]
def rot(P,a,b,th):
    ax=(b-a)/sp.sqrt((b-a).dot(b-a)); m=(a+b)/2; c,s=sp.cos(th),sp.sin(th)
    return [sp.simplify(m+(p-m)*c+ax.cross(p-m)*s+ax*(ax.dot(p-m))*(1-c)) for p in P]
nbr=[]
for i,j in itertools.combinations(range(4),2):
    for th in (2*sp.pi/3,-2*sp.pi/3):
        nbr.append(rot(T,T[i],T[j],th))
cent=[sp.simplify(sum(P,sp.zeros(3,1))/4) for P in nbr]
d2=lambda p,q:sp.simplify((p-q).dot(p-q))
print("相手の数",len(cent),"  中心までの距離²",sorted({d2(c,sp.zeros(3,1)) for c in cent}))
dd={}
for a,b in itertools.combinations(range(12),2):
    k=d2(cent[a],cent[b]); dd[k]=dd.get(k,0)+1
print("12個の中心どうしの距離²の分布:",{str(k):v for k,v in sorted(dd.items(),key=lambda t:float(sp.N(t[0])))})
print("参照 立方八面体（半径²=3）：辺²3が24組, 6が24, 9が12, 12が6 ／ 正二十面体：辺²は 3(… ) で30組")
# 相手どうしの重なり（浮動小数、分離軸で判定）
def tet(P): return np.array([[float(sp.N(x)) for x in p] for p in P])
Ps=[tet(T)]+[tet(P) for P in nbr]
def faces(P):
    out=[]
    for f in itertools.combinations(range(4),3):
        n=np.cross(P[f[1]]-P[f[0]],P[f[2]]-P[f[0]]); out.append(n)
    return out
def overlap(P,Q):
    axes=faces(P)+faces(Q)
    eP=[P[j]-P[i] for i,j in itertools.combinations(range(4),2)]; eQ=[Q[j]-Q[i] for i,j in itertools.combinations(range(4),2)]
    axes+= [np.cross(a,b) for a in eP for b in eQ]
    for n in axes:
        if np.linalg.norm(n)<1e-12: continue
        a,b=P@n,Q@n
        if a.max()<=b.min()+1e-9 or b.max()<=a.min()+1e-9: return False
    return True
bad=[(i,j) for i,j in itertools.combinations(range(13),2) if overlap(Ps[i],Ps[j])]
print("重なる組（0＝真ん中）:",len(bad),bad[:20])
# 真ん中の頂点のまわり：各頂点に何個の正四面体が集まるか
V0=Ps[0]
for k,v in enumerate(V0):
    n=sum(any(np.linalg.norm(Q-v,axis=1)<1e-9) for Q in Ps)
    print(f"  真ん中の頂点{k}に集まる正四面体の数 {n}")
np.save('cent12.npy',np.array([[float(sp.N(x)) for x in c] for c in cent]))

# 重ならない相手の選び方をすべて数え、辺のまわりが閉じる（±120°の両方が入る）辺が最大いくつかを見る
edges=list(itertools.combinations(range(4),2))
badset={frozenset((i-1,j-1)) for i,j in bad}           # 相手どうしの重なり（真ん中とは重ならない）
print("真ん中と重なる相手:",[j for i,j in bad if i==0])
best={}
for mask in range(1<<12):
    S=[k for k in range(12) if mask>>k&1]
    if any(frozenset((a,b)) in badset for a,b in itertools.combinations(S,2)): continue
    closed=tuple(e for n,e in enumerate(edges) if 2*n in S and 2*n+1 in S)
    key=len(closed); best.setdefault(key,set()).add((len(S),closed))
for k in sorted(best):
    ex=sorted(best[k],key=lambda t:-t[0])
    print(f"閉じる辺 {k} 本：置ける相手の最大 {ex[0][0]}  例 {sorted({c for _,c in ex})[:4]}")
# 重なる3つ組は、真ん中のどの面の辺か
for tri in [(1,4,7),(2,5,10),(3,6,11),(8,9,12)]:
    es=[edges[(k-1)//2] for k in tri]; vs=set(itertools.chain(*es))
    print("重なる3つ組",tri,"→ 辺",es," 面",sorted(vs))
