# トリプチセンとペンチプセン：正八面体と正四面体で空間を埋める並び（その4の結晶と同じ）の上で組む。整数だけ
# 星型正八面体＝正八面体の中心 o と、その 8 面の正四面体。正四面体の中心＝炭素は o+3s、とがった先は o+6s（s は ±1 の 3 つ組）
# 炭素どうしの結合＝正四面体が辺を分け合う（中心どうしの長さ² 36）
from itertools import product, combinations
dot=lambda a,b:sum(x*y for x,y in zip(a,b)); add=lambda a,b:tuple(x+y for x,y in zip(a,b)); sub=lambda a,b:tuple(x-y for x,y in zip(a,b))
neg=lambda a:tuple(-x for x in a); d2=lambda a,b:dot(sub(a,b),sub(a,b))
SG=list(product((1,-1),repeat=3))
def tet(o,s): return frozenset([tuple(o[i]+(6*s[k] if i==k else 0) for i in range(3)) for k in range(3)]+[tuple(o[i]+6*s[i] for i in range(3))])
def stella(o): return {add(o,tuple(3*x for x in s)):(tet(o,s),add(o,tuple(6*x for x in s))) for s in SG}   # 炭素 → (正四面体, とがった先)
def ring(o,axis):
    """中心 o の星型正八面体で、軸 axis（±1 の 3 つ組）の両端を外した 6 個を輪の順に"""
    rest=[c for c in stella(o) if abs(dot(sub(c,o),axis))!=9]; r=[rest[0]]
    while len(r)<6: r.append([c for c in rest if c not in r and d2(c,r[-1])==36][0])
    return r
def other_oct(o,c1,c2):
    """炭素 c1・c2 の正四面体（隣どうし）を、o とは別の正八面体の面に持つ星型正八面体の中心"""
    T1=stella(o)[c1][0]; T2=stella(o)[c2][0]; e=T1&T2
    F=[v for v in T1 if v in e]+[stella(o)[c1][1]]          # 分けた辺＋c1 のとがった先
    t=[v for v in T1 if v not in F][0]
    o2=sub(tuple(sum(v[i] for v in F) for i in range(3)),t); o3=tuple(x//2 for x in o2)
    assert c1 in stella(o3) and c2 in stella(o3) and stella(o3)[c1][0]==T1 and stella(o3)[c2][0]==T2
    return o3
class Mol:
    def __init__(s): s.C={}; s.bonds=set(); s.H={}
    def tet_of(s,c): return s.C[c]
    def add_C(s,c,T,tip):
        if c in s.C: assert s.C[c][0]==T,"同じ場所に別の正四面体"; return
        s.C[c]=(T,tip)
    def bond(s,a,b):
        assert len(s.C[a][0]&s.C[b][0])==2,"辺を分け合っていない"; s.bonds.add(frozenset((a,b)))
def core(m,o,P,Q,pairs_hint=None):
    """中心 o の星型正八面体で、P を C9・Q を C10 とする芯。重ならない 3 本の道（2 通りのうち、hint があればそれに合う方）"""
    S=stella(o); NP=[c for c in S if d2(c,P)==36]; NQ=[c for c in S if d2(c,Q)==36]
    opts=[]
    for perm in __import__("itertools").permutations(NQ):
        if all(d2(n,q)==36 for n,q in zip(NP,perm)): opts.append(list(zip(NP,perm)))
    if pairs_hint: opts=[op for op in opts if all(any(set(h)=={x,y} for x,y in op) for h in pairs_hint)]
    op=opts[0]
    for c in S: m.add_C(c,*S[c])
    m.Hsite=getattr(m,"Hsite",set())|{P,Q}
    for n,q in op: m.bond(P,n); m.bond(n,q); m.bond(q,Q)
    return op
def blade(m,o_core,a,b,avoid_axis_ends):
    """a–b を隣どうしに持つ別の星型正八面体のベンゼン（帯）。帯の取り方は重ならない方"""
    o=other_oct(o_core,a,b); S=stella(o)
    for ax in SG[:4]:
        r=ring(o,ax)
        if a in r and b in r and abs(r.index(a)-r.index(b)) in (1,5):
            new=[c for c in r if c not in (a,b)]
            if any(c in m.C and m.C[c][0]!=S[c][0] for c in new): continue
            if any(c in m.C for c in new): continue
            for c in r: m.add_C(c,*S[c])
            for i in range(6): m.bond(r[i],r[(i+1)%6])
            m.Hsite|=set(new); return o,r
    raise Exception("羽根が付かない")
def finish(m,name):
    deg={c:sum(1 for e in m.bonds if c in e) for c in m.C}
    H=[m.C[c][1] for c in m.C if c in m.Hsite]
    clash=set(H)&set(m.C)
    n_tet=len({m.C[c][0] for c in m.C})
    print(name,"：C",len(m.C),"H",len(H),"結合",len(m.bonds),"／ 炭素の結合の数",sorted(__import__("collections").Counter(deg.values()).items()),
          "／ H の重なり",len(H)-len(set(H)),"H と炭素の重なり",len(clash),"正四面体",n_tet)
    return H
# ---- トリプチセン ----
m=Mol(); O=(0,0,0); P,Q=(3,3,-3),(-3,-3,3)
op=core(m,O,P,Q)
bl=[blade(m,O,n,q,None) for n,q in op]
print("トリプチセンの羽根の正八面体",[b[0] for b in bl])
finish(m,"トリプチセン（C20H14 になるはず）")
# ---- ペンチプセン：羽根 1 を真ん中の環にして、その向かい側の結合に 2 つ目の芯 ----
o1,r1=bl[0]; a,b=op[0]; i,j=r1.index(a),r1.index(b)
d,e=r1[(i+3)%6],r1[(j+3)%6]
print("真ん中の環",r1,"芯 1 の結合",(a,b),"向かい側の結合",(d,e))
o2=other_oct(o1,d,e); S2=stella(o2)
print("2 つ目の芯の正八面体の中心",o2)
# 芯 2：d・e が「C9 の隣・C10 の隣」になる対角線を探す
found=None
for P2 in S2:
    for Q2 in S2:
        if dot(sub(P2,o2),sub(Q2,o2))!=-27: continue
        if d2(P2,d)==36 and d2(Q2,e)==36 and P2 not in m.C and Q2 not in m.C: found=(P2,Q2); break
    if found: break
print("芯 2 の C9・C10",found)
if found:
    P2,Q2=found
    newc=[c for c in S2 if c not in (d,e)]
    print("芯 2 の新しい炭素のうち既にある炭素と重なるもの",[c for c in newc if c in m.C])
    op2=core(m,o2,P2,Q2,pairs_hint=[(d,e)]); m.Hsite-={d,e}
    bl2=[blade(m,o2,n,q,None) for n,q in op2 if {n,q}!={d,e}]
    print("ペンチプセンの外側の羽根の正八面体",[b[0] for b in bl2])
    finish(m,"ペンチプセン（C34H22 になるはず）")
    import json; json.dump({"C":[list(c) for c in m.C],"bonds":[[list(x) for x in e] for e in m.bonds],
        "H":[[list(c),list(m.C[c][1])] for c in m.C if c in m.Hsite],"octs":[O,o1,bl[1][0],bl[2][0],o2]+[b[0] for b in bl2],
        "cores":[[list(P),list(Q)],[list(P2),list(Q2)]],"T":{",".join(map(str,c)):[list(v) for v in m.C[c][0]] for c in m.C}},open("pentiptycene.json","w"))
