exec(open('grow2.py').read().split("tets={key(T0):T0}")[0])
# 世代ごとに置き、置くたびに重なり・一致を調べる。重なったら止める。
tets=[(T0,0)]; keys={key(T0)}; closed=set()
frontier=[]
def close(P,i,j,gen):
    ek=edge_key(P,i,j)
    if ek in closed: return []
    closed.add(ek); out=[]
    for th in (2*np.pi/3,-2*np.pi/3):
        Q=rot(P,P[i],P[j],th); out.append((Q,ek,gen))
    return out
cand=close(T0,0,1,1)+close(T0,2,3,1)
g=1
while cand and g<=6:
    placed=0; coinc=0; ov=[]
    nxt=[]
    for Q,ek,gen in cand:
        k=key(Q)
        if k in keys: coinc+=1; continue
        hits=[t for t,(P,gg) in enumerate(tets) if overlap(P,Q)]
        if hits: ov.append((gen,[tets[h][1] for h in hits])); continue
        tets.append((Q,gen)); keys.add(k); placed+=1
        for (i,j),(a,b) in opp.items():
            if edge_key(Q,i,j)==ek: nxt+=close(Q,a,b,gen+1); break
    print(f"世代{g}: 置けた {placed}  既存と一致 {coinc}  重なって置けない {len(ov)}  （ぶつかった相手の世代 {sorted({x for _,h in ov for x in h})}）")
    cand=nxt; g+=1
cent=np.array([P.mean(0) for P,_ in tets])
print("置けた総数",len(tets))
np.save('grow_cent.npy',cent); import pickle; pickle.dump(tets,open('tets.pkl','wb'))
