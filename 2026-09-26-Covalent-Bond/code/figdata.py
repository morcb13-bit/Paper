import numpy as np, itertools, json, pickle
exec(open('grow2.py').read().split("tets={key(T0):T0}")[0])
# 枝分かれ：世代4までを、重なりの印つきで
tets=[(T0,0,False)]; keys={key(T0)}; closed=set()
def close(P,i,j,gen):
    ek=edge_key(P,i,j)
    if ek in closed: return []
    closed.add(ek); return [(rot(P,P[i],P[j],th),ek,gen) for th in (2*np.pi/3,-2*np.pi/3)]
cand=close(T0,0,1,1)+close(T0,2,3,1)
for g in range(1,5):
    nxt=[]
    for Q,ek,gen in cand:
        if key(Q) in keys: continue
        hit=any(overlap(P,Q) for P,_,_ in tets)
        tets.append((Q,gen,hit)); keys.add(key(Q))
        if not hit:
            for (i,j),(a,b) in opp.items():
                if edge_key(Q,i,j)==ek: nxt+=close(Q,a,b,gen+1); break
    cand=nxt
print("枝分かれ：",[(g,sum(1 for _,gg,h in tets if gg==g and not h),sum(1 for _,gg,h in tets if gg==g and h)) for g in range(5)])
tree=[{"v":P.tolist(),"g":g,"hit":h} for P,g,h in tets]
# ねじの一本道：一歩 S（T0 を 次の辺が 0,1 に来るよう読み替えた T1 へ）
T1=rot(T0,T0[0],T0[1],2*np.pi/3); T1r=T1[[2,3,0,1]]
A=np.stack([T0[i]-T0[0] for i in (1,2,3)],1); B=np.stack([T1r[i]-T1r[0] for i in (1,2,3)],1)
R=B@np.linalg.inv(A); t=T1r[0]-R@T0[0]
print("一歩の cos",round((np.trace(R)-1)/2,12))
chain=[T0]
for k in range(29): chain.append(np.array([R@p+t for p in chain[-1]]))
ov=[(a,b) for a,b in itertools.combinations(range(len(chain)),2) if overlap(chain[a],chain[b])]
print("一本道 30 個の中の重なり",len(ov),ov[:5])
w,v=np.linalg.eig(R); axis=np.real(v[:,np.argmin(abs(w-1))]); axis/=np.linalg.norm(axis)
json.dump({"tree":tree,"chain":[P.tolist() for P in chain],"axis":axis.tolist()},open('fig2.json','w'))
