# 規則：各正四面体は「すでに閉じている辺」と向かい合う辺でも、±120°の二つを置いて一周を閉じる。
# これを繰り返し、重なり・一致（輪が閉じる）・広がり方を見る。（浮動小数、許容 1e-9）
import numpy as np, itertools
T0=np.array([(1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1)],float)
def rot(P,a,b,th):
    ax=(b-a)/np.linalg.norm(b-a); m=(a+b)/2; c,s=np.cos(th),np.sin(th)
    Q=P-m; return m+Q*c+np.cross(ax,Q)*s+np.outer(Q@ax,ax)*(1-c)
def key(P): return tuple(np.round(P.mean(0),6))
def overlap(P,Q):
    fs=lambda X:[np.cross(X[f[1]]-X[f[0]],X[f[2]]-X[f[0]]) for f in itertools.combinations(range(4),3)]
    eP=[P[j]-P[i] for i,j in itertools.combinations(range(4),2)]; eQ=[Q[j]-Q[i] for i,j in itertools.combinations(range(4),2)]
    for n in fs(P)+fs(Q)+[np.cross(a,b) for a in eP for b in eQ]:
        if np.linalg.norm(n)<1e-12: continue
        a,b=P@n,Q@n
        if a.max()<=b.min()+1e-9 or b.max()<=a.min()+1e-9: return False
    return True
def edge_key(P,i,j): return tuple(sorted([tuple(np.round(P[i],6)),tuple(np.round(P[j],6))]))
opp={(0,1):(2,3),(2,3):(0,1),(0,2):(1,3),(1,3):(0,2),(0,3):(1,2),(1,2):(0,3)}
tets={key(T0):T0}
closed_edges=set()
# 種：真ん中の辺 0-1 を閉じる
queue=[]
def close(P,i,j,gen):
    ek=edge_key(P,i,j)
    if ek in closed_edges: return []
    closed_edges.add(ek); new=[]
    for th in (2*np.pi/3,-2*np.pi/3):
        Q=rot(P,P[i],P[j],th); k=key(Q)
        if k in tets: continue
        tets[k]=Q; new.append((Q,ek,gen))
    return new
queue+=close(T0,0,1,1)
# 真ん中自身も、向かい合う辺 2-3 を閉じる
queue+=close(T0,2,3,1)
log=[]
gen_count={0:1}
while queue:
    Q,ek,gen=queue.pop(0)
    if len(tets)>400: break
    gen_count[gen]=gen_count.get(gen,0)+1
    # この四面体がどの辺で閉じたか → 向かい合う辺を閉じる
    for (i,j),(a,b) in opp.items():
        if edge_key(Q,i,j)==ek:
            queue+=close(Q,a,b,gen+1); break
allT=list(tets.values())
bad=[(x,y) for x,y in itertools.combinations(range(len(allT)),2) if overlap(allT[x],allT[y])]
cent=np.array([P.mean(0) for P in allT])
print("正四面体の数",len(allT),"  閉じた辺",len(closed_edges),"  世代ごとの数",gen_count)
print("重なる組",len(bad))
r=np.linalg.norm(cent,axis=1); print("中心からの距離の最大",r.max().round(4))
np.save('grow_cent.npy',cent)
