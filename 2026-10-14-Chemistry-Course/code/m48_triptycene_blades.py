# トリプチセンの本当の姿：芯（星型正八面体の 8 個）に、羽根（別の星型正八面体のベンゼン）を 3 枚付ける。整数だけ
# 羽根は芯の a・b の正四面体を、帯の隣どうしとして共有する。羽根の正八面体は a・b の正四面体の面から決まる
from itertools import product, combinations
dot=lambda a,b:sum(x*y for x,y in zip(a,b)); add=lambda a,b:tuple(x+y for x,y in zip(a,b)); sub=lambda a,b:tuple(x-y for x,y in zip(a,b))
neg=lambda a:tuple(-x for x in a); d2=lambda a,b:dot(sub(a,b),sub(a,b))
cen=lambda T:tuple(sum(v[i] for v in T)//4 for i in range(3))
# 芯：m41 の星型正八面体（正八面体の中心 0）
SG=list(product((1,-1),repeat=3))
CORE={s:frozenset([tuple(6*s[k] if i==k else 0 for i in range(3)) for k in range(3)]+[tuple(6*x for x in s)]) for s in SG}
P,Q=(1,1,-1),(-1,-1,1)                       # C9・C10（への字の頂点と底）
PATHS=[((1,1,1),(1,-1,1)),((1,-1,-1),(-1,-1,-1)),((-1,1,-1),(-1,1,1))]   # a→b（m47 の 1 通り目）
def stella_on_face(F,T_on):
    """面 F の上に T_on が立っている正八面体（T_on の反対側）と、その 8 面の正四面体"""
    t=[v for v in T_on if v not in F][0]
    g3=tuple(sum(v[i] for v in F) for i in range(3))           # 3g
    o2=sub(g3,t)                                                # 2o = 3g − t
    if any(x%2 for x in o2): return None
    o=tuple(x//2 for x in o2)
    OCT=set(F)|{sub(o2,v) for v in F}
    if len(OCT)!=6 or any(d2(a,b) not in (72,144) for a,b in combinations(OCT,2)): return None
    faces=[f for f in combinations(OCT,3) if all(d2(a,b)==72 for a,b in combinations(f,2))]
    tets=[]
    for f in faces:
        g3f=tuple(sum(v[i] for v in f) for i in range(3)); tip=sub(g3f,o2)   # 3g − 2o
        tets.append(frozenset(f)|{tip})
    return o,OCT,tets
results=[]
for k,(sa,sb) in enumerate(PATHS):
    Ta,Tb=CORE[sa],CORE[sb]; e=Ta&Tb
    assert len(e)==2
    opts=[]
    fa_list=[frozenset(f) for f in combinations(Ta,3) if e<=frozenset(f)]
    fb_list=[frozenset(f) for f in combinations(Tb,3) if e<=frozenset(f)]
    for Fa in fa_list:
        r=stella_on_face(Fa,Ta)
        if r is None: continue
        o,OCT,tets=r
        if Ta not in tets or Tb not in tets: continue
        if set(OCT)==set(v for s in SG for v in CORE[s] if sum(abs(x) for x in v)==6): continue   # 芯そのもの
        opts.append((o,OCT,tets))
    print("羽根",k+1,"：a",tuple(3*x for x in sa),"b",tuple(3*x for x in sb),"芯以外の正八面体",len(opts),[o for o,_,_ in opts])
    results.append(opts)
# 羽根の帯：a・b を隣どうしに含む 3 回軸（立方体の対角線のうち a・b を通らないもの）
def belts(o,tets,Ta,Tb):
    C={T:cen(T) for T in tets}; out=[]
    for T1,T2 in combinations(tets,2):
        if sub(C[T1],o)!=neg(sub(C[T2],o)): continue          # 対角線の両端＝ふた
        if Ta in (T1,T2) or Tb in (T1,T2): continue
        out.append(((T1,T2),[T for T in tets if T not in (T1,T2)]))
    return out
core_tets=set(CORE.values())
best=[]
for combo in product(*results):
    ok=True; used=set(core_tets); blades=[]
    for k,(o,OCT,tets) in enumerate(combo):
        sa,sb=PATHS[k]; new=[T for T in tets if T not in (CORE[sa],CORE[sb])]
        if any(T in used for T in new): ok=False;break
        used|=set(new); blades.append((o,OCT,tets))
    if ok: best.append(blades)
print("3 枚の羽根が芯とも互いとも正四面体を重ねずに付く組",len(best))
for blades in best:
    print("  羽根の正八面体の中心",[b[0] for b in blades])
    for k,(o,OCT,tets) in enumerate(blades):
        sa,sb=PATHS[k]; B=belts(o,tets,CORE[sa],CORE[sb])
        print("   羽根",k+1,"帯の取り方",len(B))
import pickle; pickle.dump(best,open("blades.pkl","wb"))
