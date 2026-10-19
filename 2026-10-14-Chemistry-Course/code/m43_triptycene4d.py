# トリプチセンを 4 次元で閉じる：星型正八面体のベンゼン環（m41）を、24 胞体の向きで 4 次元に置く。整数だけ
# 世界座標は 3 次元の座標の 2 倍の単位（長さ² は 4 倍）。正四面体：中心から頂点の長さ² 108、頂点どうしの内積 −36
import sys
from itertools import product, permutations, combinations
from collections import defaultdict
dot=lambda a,b:sum(x*y for x,y in zip(a,b)); add=lambda a,b:tuple(x+y for x,y in zip(a,b)); sub=lambda a,b:tuple(x-y for x,y in zip(a,b))
neg=lambda a:tuple(-x for x in a)
# ---- 星型正八面体のベンゼン（m41 と同じ、3 次元の座標） ----
SG=list(product((1,-1),repeat=3))
TET={s:[tuple(6*s[k] if i==k else 0 for i in range(3)) for k in range(3)]+[tuple(6*x for x in s)] for s in SG}
RING=[(1,1,-1),(1,-1,-1),(1,-1,1),(-1,-1,1),(-1,1,1),(-1,1,-1)]
CEN={s:tuple(3*x for x in s) for s in SG}; TIP={s:tuple(6*x for x in s) for s in SG}
ATOMS3=[CEN[s] for s in RING]
# ---- 4 次元の向き：24 胞体の 24 本（長さ² 4）から直交する 3 本 ----
V24=[v for v in product((-2,-1,0,1,2),repeat=4) if dot(v,v)==4 and (all(abs(x)==1 for x in v) or sorted(map(abs,v))==[0,0,0,2])]
assert len(V24)==24
def frames(only3d):
    vs=[v for v in V24 if (not only3d) or (v[3]==0 and all(abs(x) in (0,2) for x in v))]
    for f1 in vs:
        for f2 in vs:
            if dot(f1,f2): continue
            for f3 in vs:
                if dot(f1,f3) or dot(f2,f3): continue
                yield (f1,f2,f3)
R=lambda F,x:tuple(x[0]*F[0][i]+x[1]*F[1][i]+x[2]*F[2][i] for i in range(4))   # 3 次元 → 4 次元（2 倍）
def run(only3d):
    # 環 1 個の置き方：C9 を原点に、帯の炭素 c のとがった先が C9。隣 b のとがった先が C10
    table=defaultdict(list)          # C10 → [(u=a, b, 環の原子の集まり, 枠)]
    for F in frames(only3d):
        for i,c in enumerate(RING):
            T=neg(R(F,TIP[c]))                                  # 環の平行移動
            a=add(T,R(F,CEN[c]))
            for j in ((i+1)%6,(i-1)%6):
                b=add(T,R(F,CEN[RING[j]])); C10=add(T,R(F,TIP[RING[j]]))
                atoms=frozenset(add(T,R(F,p)) for p in ATOMS3)
                hs=frozenset(add(T,R(F,TIP[s])) for k,s in enumerate(RING) if k not in (i,j))
                table[C10].append((a,b,atoms,hs))
    sols=[]
    for C10,L in table.items():
        for x,y,z in combinations(L,3):
            us=[x[0],y[0],z[0]]; bs=[x[1],y[1],z[1]]
            if any(dot(p,q)!=-36 for p,q in combinations(us,2)): continue          # C9 の正四面体
            ws=[sub(b,C10) for b in bs]
            if any(dot(p,q)!=-36 for p,q in combinations(ws,2)): continue          # C10 の正四面体
            if len(x[2]|y[2]|z[2])!=18: continue                                    # 炭素が重ならない
            H9=neg(add(add(us[0],us[1]),us[2])); H10=sub(C10,add(add(ws[0],ws[1]),ws[2]))
            allp=list(x[2]|y[2]|z[2])+list(x[3]|y[3]|z[3])+[(0,0,0,0),C10,H9,H10]
            if len(set(allp))!=len(allp): continue                                  # どの原子も重ならない
            sols.append((C10,(x,y,z),H9,H10))
    return table,sols
for name,only3d in [("Q0 3 次元の向きだけ",True),("Q1 4 次元の向き（24 胞体）",False)]:
    table,sols=run(only3d)
    c10=set(s[0] for s in sols)
    print(name,"：環 1 個の置き方",sum(len(v) for v in table.values()),"／ 3 環が閉じる選び方",len(sols),"／ C10 の場所",len(c10))
    if sols:
        C10,(x,y,z),H9,H10=sols[0]
        print("   例：C10",C10,"C9–C10 の長さ²",dot(C10,C10),"（3 次元の単位で",dot(C10,C10)//4,"）")
        print("   a",[t[0] for t in (x,y,z)],"b",[t[1] for t in (x,y,z)],"H9",H9,"H10",H10)
        print("   C10 の長さ²の種類",sorted({dot(s[0],s[0]) for s in sols}))
        import json; json.dump([{"C10":s[0],"rings":[{"a":r[0],"b":r[1],"C":sorted(r[2]),"H":sorted(r[3])} for r in s[1]],"H9":s[2],"H10":s[3]} for s in sols[:200]],open("trip4d.json","w"))
