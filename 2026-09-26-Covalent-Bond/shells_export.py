exec(open('h2walk.py').read().split("# 照合")[0])
import json
def sq(w): return dot(w,w)
def vsub(u,v): return tuple(sub(u[i],v[i]) for i in range(3))
O=(Z,Z,Z); OB=refl(O)
# 頂点：A の20個＋B の新しい15個
allv=list(V); owner=["A"]*20
for v in VB:
    if v not in V: allv.append(v); owner.append("B")
for i,v in enumerate(allv):
    if v in [V[k] for k in face]: owner[i]="S"
N=len(allv)
key=lambda v: (sq(vsub(v,O)), sq(vsub(v,OB)))
# 殻：{d1², d2²} の組（順序なし）で決まる。値で並べるため実数は並べ替えにだけ使う
phi=(1+5**0.5)/2; val=lambda x:x[0]+x[1]*phi
shellkey=lambda v: tuple(sorted([sq(vsub(v,O)),sq(vsub(v,OB))]))
keys=sorted(set(shellkey(v) for v in allv), key=lambda k: val(k[1]))
shell=[keys.index(shellkey(v)) for v in allv]
print("殻の数",len(keys),"各殻の頂点数",[shell.count(s) for s in range(len(keys))])
idx={v:i for i,v in enumerate(allv)}
# 鏡（共有面）：refl
mir=[idx[refl(v)] for v in allv]
assert all(shell[mir[i]]==shell[i] for i in range(N))
# 各原子の中心について反対（B は自分の中心 OB について）
def anti(i):
    v=allv[i]
    if owner[i]=="B": return idx[vsub(OB,v)] if False else idx[tuple(sub(OB[k],v[k]) if False else sub(OB[k],sub(v[k],OB[k])) for k in range(3))]
    return idx[tuple(neg(c) for c in v)]
antiA=[idx[tuple(neg(c) for c in v)] if owner[i]!="B" else None for i,v in enumerate(allv)]
antiB=[]
for i,v in enumerate(allv):
    w=tuple(sub(sc2,c) for sc2,c in zip(tuple(add(OB[k],OB[k]) for k in range(3)),v))
    antiB.append(idx.get(w))
# 検査：A 側（共有面含む）の対蹠は A 内、B 側（共有面含む）の対蹠は B 内。殻は 0↔3, 1↔2
ok=True
for i in range(N):
    if owner[i] in "AS":
        j=antiA[i]; ok&= shell[i]+shell[j]==3
    if owner[i] in "BS":
        j=antiB[i]; ok&= (j is not None) and shell[i]+shell[j]==3
print("中心の対が殻1↔4・2↔3に分かれる:", ok)
# 辺・共有面
E=set()
for i in range(N):
    for j in range(i+1,N):
        if sq(vsub(allv[i],allv[j]))==EDGE2: E.add((i,j))
print("辺",len(E))
fl=lambda x:(x[0]+x[1]*phi)/5
coords=[[fl(c) for c in v] for v in allv]
s_vals=[]
for k in keys:
    a,b=k; s_vals.append((val(a)**0.5+val(b)**0.5)/5)
json.dump({"coords":coords,"edges":sorted(E),"owner":owner,"shell":shell,"mir":mir,
           "antiA":antiA,"antiB":antiB,"OB":[fl(c) for c in OB],"sums":s_vals},open('shells.json','w'))
