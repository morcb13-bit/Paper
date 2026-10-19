# 段 2：段 1 の解（C10≠C9）ごとに、各環の残り 4 個（c,d,e,f）が辺の共有で 6 員環に閉じるか
# 帯の決めごと：環の正四面体は、底の面（とがった先の反対の面）の辺 2 本で前後の隣と辺を分け合う。とがった先は誰とも分けない（H、a・b では C9・C10）
import pickle, sys
from itertools import combinations
exec(open("m45_hinge.py").read().split("C9=(0,0,0,0)")[0].replace("DIM4 = sys.argv[1]==\"4d\" if len(sys.argv)>1 else True","DIM4=True"))
C9=(0,0,0,0)
exec(open("hinge_fast.py").read())
print("辺の向き",len(D))
sols=[s for s in pickle.load(open("hinge_sols_4d.pkl","rb")) if dot(s[0],s[0])>0]
print("段 1 の解（C10≠C9）",len(sols))
E=lambda T:{frozenset(e) for e in combinations(T,2)}
cache={}
def twe(p,q):
    k=frozenset((p,q))
    if k not in cache: cache[k]=tets_with_edge(p,q)
    return cache[k]
def complete_ring(Ta,tipa,Tb,tipb,shared):
    """a,b を含む 6 員環の残り 4 個を探す。戻り値：見つかった環（正四面体の列）"""
    baseA=[v for v in Ta if v!=tipa]; baseB=[v for v in Tb if v!=tipb]
    endA=[frozenset(e) for e in combinations(baseA,2) if frozenset(e)!=shared]   # a が f と分ける辺
    startB=[frozenset(e) for e in combinations(baseB,2) if frozenset(e)!=shared] # b が c と分ける辺
    found=[]
    def step(path,T,tip,used_edge,depth):
        base=[v for v in T if v!=tip]
        if depth==4:
            return
        for e in combinations(base,2):
            e=frozenset(e)
            if e==used_edge: continue
            if not (e & used_edge): continue          # 底の 2 本は 1 頂点を共有する
            p,q=tuple(e)
            for c,Tc in twe(p,q):
                if Tc in [x[0] for x in path] or Tc in (Ta,Tb): continue
                tipc=[v for v in Tc if v not in e]
                for t in tipc:                          # とがった先の選び方：e と同じ面にない方は後で決まる
                    nxt=path+[(Tc,t,e)]
                    if depth==3:
                        # f は a と endA の辺を分ける
                        for ea in endA:
                            if ea <= Tc and t not in ea and (ea & e):
                                found.append(nxt)
                    else:
                        step(nxt,Tc,t,e,depth+1)
    for eb in startB:
        if not (eb & shared): continue
        p,q=tuple(eb)
        for c,Tc in twe(p,q):
            if Tc in (Ta,Tb): continue
            for t in [v for v in Tc if v not in eb]:
                step([(Tc,t,eb)],Tc,t,eb,1)
    return found
ok=0; rec=[]
for n,(C10,x,y,z) in enumerate(sols):
    res=[]
    for ch in (x,y,z):
        a,Ta,pq,b,Tb=ch
        rings=complete_ring(Ta,C9,Tb,C10,frozenset(pq))
        res.append(len(rings))
    if all(res): ok+=1; rec.append((C10,res))
    print("  解",n,res,flush=True) if n%10==0 else None
    if n<3 or all(res): print("解",n,"C10",C10,"長さ²",dot(C10,C10),"環ごとに閉じる道の数",res)
print("段 2：3 つの環がすべて 6 員環に閉じる解",ok,"／",len(sols))
pickle.dump(rec,open("ring_close.pkl","wb"))
