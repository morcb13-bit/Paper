exec(open('h2walk.py').read().split("# 出発点")[0])
import json, math
phi=(1+5**0.5)/2
fl=lambda x:(x[0]+x[1]*phi)/5
adj,region,ne=build(set(face))
N=len(adj)
# 座標（表示用の浮動小数）
coords=[None]*N
for i in range(20): coords[i]=[fl(c) for c in V[i]]
idB={}
for k,v in enumerate(VB):
    if v in V and V.index(v) in face: idB[k]=V.index(v)
# build と同じ順で B 専用頂点に番号を振り直す
nxt=20
for k,v in enumerate(VB):
    if k not in idB:
        idB[k]=nxt; coords[nxt]=[fl(c) for c in v]; nxt+=1
E=sorted({tuple(sorted(e)) for v in adj for e in [(v,u) for u in adj[v]]})
far=[i for i,v in enumerate(V) if dot(v,n)==neg(PHI2)]
sA=far[0]; kB=V.index(V[sA]); sB=idB[kB]   # A の出発点の鏡像（B 側）
reg=[region[i] for i in range(N)]
# 検査：build の region と合っているか
assert sum(r=="B" for r in reg)==15 and len(E)==55 and reg[sB]=="B"
json.dump({"coords":coords,"edges":E,"region":reg,"sA":sA,"sB":sB},open('h2data.json','w'))
print(N,len(E),sA,sB)
# 鏡像の対応：A の頂点 k の鏡像は B 側の idB[k]、共有面の頂点は自分自身
mir=[None]*N
for k in range(20):
    mir[k]=idB[k]; mir[idB[k]]=k
assert all(m is not None for m in mir)
assert all(mir[mir[i]]==i for i in range(N))                      # 二回で元に戻る
assert all((mir[i]==i)==(reg[i]=="S") for i in range(N))          # 動かないのは共有面だけ
assert all(tuple(sorted((mir[i],mir[j]))) in set(E) for i,j in E)  # 辺を辺へ写す
assert mir[sA]==sB
dd=json.load(open('h2data.json')); dd["mirror"]=mir; json.dump(dd,open('h2data.json','w'))
print("鏡像の対応 OK")
