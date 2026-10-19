# トリプチセンを閉じる：環を蝶番（分け合う辺）で 4 次元へ折ってよいとする。整数だけ
# 座標は 3 次元の 2 倍の単位：正四面体の中心→頂点の長さ² 108、頂点どうしの内積 −36、辺² 288
import sys
from itertools import product, combinations
from collections import defaultdict
dot=lambda a,b:sum(x*y for x,y in zip(a,b)); add=lambda a,b:tuple(x+y for x,y in zip(a,b)); sub=lambda a,b:tuple(x-y for x,y in zip(a,b))
neg=lambda a:tuple(-x for x in a)
DIM4 = sys.argv[1]=="4d" if len(sys.argv)>1 else True
r=10
W=[w for w in product(range(-r,r+1),repeat=4) if dot(w,w)==108 and (DIM4 or w[3]==0)]
print("中心→頂点の候補",len(W))
Wset=set(W)
def tets_with_vertex(c,v):
    """中心 c、頂点 v を持つ正四面体（頂点の集まり）"""
    w1=sub(v,c); out=set()
    if w1 not in Wset: return out
    cand=[w for w in W if dot(w,w1)==-36]
    for w2,w3 in combinations(cand,2):
        if dot(w2,w3)!=-36: continue
        w4=neg(add(add(w1,w2),w3))
        if w4 in Wset: out.add(frozenset(add(c,w) for w in (w1,w2,w3,w4)))
    return out
def tets_with_edge(p,q):
    """辺 p–q を持つ正四面体（中心と頂点）"""
    out=[]
    for w in W:
        b=sub(p,w); wq=sub(q,b)
        if wq in Wset and dot(w,wq)==-36:
            for T in tets_with_vertex(b,p):
                if q in T: out.append((b,T))
    return out
C9=(0,0,0,0)
U=[(6,-6,-6,0),(-6,6,-6,0),(-6,-6,6,0)]; H9=(6,6,6,0)
# 環ごとに：a の正四面体（頂点に C9）→ 底の辺 → b の正四面体 → C10 は b の、分けた辺以外の頂点
chains=[]
for k,a in enumerate(U):
    L=defaultdict(list)
    for Ta in tets_with_vertex(a,C9):
        base=[v for v in Ta if v!=C9]
        for p,q in combinations(base,2):
            for b,Tb in tets_with_edge(p,q):
                if b==a or Tb==Ta: continue
                for C10 in Tb:
                    if C10 in (p,q): continue
                    L[C10].append((a,Ta,(p,q),b,Tb))
    chains.append(L); print("環",k+1,"：C10 の候補の場所",len(L),"鎖",sum(len(v) for v in L.values()))
common=set(chains[0])&set(chains[1])&set(chains[2])
print("3 つの環に共通の C10",len(common))
sols=[]
for C10 in common:
    for x in chains[0][C10]:
        for y in chains[1][C10]:
            for z in chains[2][C10]:
                bs=[x[3],y[3],z[3]]
                ws=[sub(b,C10) for b in bs]
                if not all(w in Wset for w in ws): continue
                if any(dot(p,q)!=-36 for p,q in combinations(ws,2)): continue
                if len({x[3],y[3],z[3]})<3: continue
                sols.append((C10,x,y,z))
print("段 1：C10 のまわりも正四面体になる選び方",len(sols))
if sols:
    from collections import Counter
    print("  C9–C10 の長さ²（2 倍の単位）",sorted(Counter(dot(s[0],s[0]) for s in sols).items()))
    print("  C10 の第 4 成分",sorted(Counter(s[0][3] for s in sols).items()))
    import pickle; pickle.dump(sols,open("hinge_sols_%s.pkl"%("4d" if DIM4 else "3d"),"wb"))
