# ベンゼン環：黒鉛の層（m31）から六角形を 1 つ切り出す。整数だけ
# B0 環の位置  B1 ケクレの数え上げ  B2 外した 4 本目とシクロヘキサン（椅子形）
from itertools import product
dot=lambda a,b:sum(x*y for x,y in zip(a,b))
add=lambda a,b:tuple(x+y for x,y in zip(a,b))
sub=lambda a,b:tuple(x-y for x,y in zip(a,b))
neg=lambda a:tuple(-x for x in a)
s=lambda a:sum(a)                                   # 高さ（(1,1,1) の向きの目盛り）
DIA=[(3,-3,-3),(-3,3,-3),(-3,-3,3)]; UPD=(3,3,3)   # ダイヤモンドの 4 本（4 本目が UPD）
GA=[add(d,(1,1,1)) for d in DIA]; GB=[neg(d) for d in GA]

# ---- B0：原点の A から 6 歩で戻る道（同じ点を通らない） ----
def nb(p,isA): return [(add(p,d),not isA) for d in (GA if isA else GB)]
rings=set()
def walk(path,isA):
    if len(path)==6:
        if (0,0,0) in [q for q,_ in nb(path[-1],isA)]: rings.add(frozenset(path))
        return
    for q,b in nb(path[-1],isA):
        if q not in path: walk(path+[q],b)
walk([(0,0,0)],True)
print("B0 原点の炭素を通る六角形",len(rings))
R0=sorted(rings,key=lambda r:sorted(r))[0]
ctr=tuple(sum(p[i] for p in R0) for i in range(3)); assert all(c%6==0 for c in ctr)
ctr=tuple(c//6 for c in ctr); print("  六角形の中心",ctr)
# 中心を原点に、A から回る順に並べる
start=(0,0,0); ring=[start]; isA=True
while len(ring)<6:
    for q,b in nb(ring[-1],isA):
        if q in R0 and q not in ring: ring.append(q); isA=b; break
C=[sub(p,ctr) for p in ring]; KIND=["A","B"]*3
print("  炭素（中心を原点）",C,"高さ",[s(p) for p in C])
print("  中心からの距離²",[dot(p,p) for p in C],"隣どうし",[dot(sub(C[i],C[(i+1)%6]),sub(C[i],C[(i+1)%6])) for i in range(6)])
Bv=[C[i] for i in range(1,6,2)]
print("  B の 3 個がつくる三角形の辺²",[dot(sub(Bv[i],Bv[(i+1)%3]),sub(Bv[i],Bv[(i+1)%3])) for i in range(3)],"（正四面体の面 72）")
# A は三角形の中心。その三角形の頂点
for i in range(0,6,2):
    T=[add(C[i],d) for d in GA]; inring=[t for t in T if t in C]
    print("  A",C[i],"の三角形の頂点",T,"環に入るもの",len(inring))
corners=[add(C[i],d) for i in range(0,6,2) for d in GA if add(C[i],d) not in C]
print("  A の外向きの頂点（大きな三角形の角）",corners,"角どうしの辺²",dot(sub(corners[0],corners[1]),sub(corners[0],corners[1])),"＝ 72×4")
# ---- H：各炭素の、環に入らない 1 本の先 ----
H=[]
for i,p in enumerate(C):
    for d in (GA if KIND[i]=="A" else GB):
        q=add(p,d)
        if q not in C: H.append(q)
print("  H",H,"中心からの距離²",[dot(h,h) for h in H],"H の数",len(H))
print("  C–H の長さ²",[dot(sub(H[i],C[i]),sub(H[i],C[i])) for i in range(6)],"（模型では C–C と同じ 24）")

# ---- B1：ケクレ。6 本の結合から、どの炭素にもちょうど 1 本 ----
E=[(i,(i+1)%6) for i in range(6)]
kek=[m for m in product([0,1],repeat=6) if all(sum(m[j] for j,e in enumerate(E) if i in e)==1 for i in range(6))]
print("B1 64 通りのうち条件を満たすもの",len(kek),kek)
# A の三角形ごとに見る：A は環の結合 2 本のうち 1 本を選ぶ
for m in kek:
    pick=[ "後" if m[i]==1 else "前" for i in range(0,6,2)]
    print("  ケクレ",m,"A ごとの選び方（後＝次の B へ、前＝一つ前の B へ）",pick)
# 鏡映（x と y の入れ替え）・回転（座標の巡回）で 2 通りがどう移るか
def mapping(f):
    idx=[C.index(f(p)) for p in C]          # 炭素の行き先
    def img(m):
        out=[0]*6
        for j,(a,b) in enumerate(E):
            if m[j]:
                ia,ib=idx[a],idx[b]; k=E.index((ia,ib)) if (ia,ib) in E else E.index((ib,ia)); out[k]=1
        return tuple(out)
    return idx,img
for name,f in [("鏡映 x↔y",lambda p:(p[1],p[0],p[2])),("鏡映 y↔z",lambda p:(p[0],p[2],p[1])),
               ("120 度（座標の巡回）",lambda p:(p[1],p[2],p[0])),("60 度（巡回して反転）",lambda p:neg((p[2],p[0],p[1])))]:
    try:
        idx,img=mapping(f); print("  ",name,"炭素",idx,"ケクレ 1 →",img(kek[0])==kek[0] and "そのまま" or (img(kek[0])==kek[1] and "もう一方" or "外れる"))
    except ValueError: print("  ",name,"六角形が自分に移らない")

# ---- B2：外した 4 本目とシクロヘキサン ----
AX=[UPD if k=="A" else neg(UPD) for k in KIND]
print("B2 外した 4 本目の向き（高さの増え）",[s(a) for a in AX],"上",sum(1 for a in AX if s(a)>0),"下",sum(1 for a in AX if s(a)<0))
# 椅子形：A の炭素と B の H を (1,1,1) 上げる（m31 で下ろした分を戻す）
CH=[add(p,(1,1,1)) if k=="A" else p for p,k in zip(C,KIND)]
HE=[h if k=="A" else add(h,(1,1,1)) for h,k in zip(H,KIND)]
HX=[add(p,a) for p,a in zip(CH,AX)]
print("  椅子形の炭素の高さ",[s(p) for p in CH])
print("  C–C の長さ²",[dot(sub(CH[i],CH[(i+1)%6]),sub(CH[i],CH[(i+1)%6])) for i in range(6)],"（ダイヤモンド 27）")
print("  C–C の向きはダイヤモンドの 4 本か",all(sub(CH[(i+1)%6],CH[i]) in (DIA+[UPD] if KIND[i]=="A" else [neg(d) for d in DIA+[UPD]]) for i in range(6)))
print("  横の H の高さ",[s(h) for h in HE],"縦の H の高さ",[s(h) for h in HX])
print("  C–H（横）の長さ²",[dot(sub(HE[i],CH[i]),sub(HE[i],CH[i])) for i in range(6)],"（縦）",[dot(sub(HX[i],CH[i]),sub(HX[i],CH[i])) for i in range(6)])
print("  H の数：ベンゼン",len(H),"シクロヘキサン",len(HE)+len(HX),"差",len(HX),"＝ 外した 4 本目の数")
# 黒鉛の中で：A の三角形 3 個と真ん中の穴（シェルピンスキー 2 段目）
import json
json.dump({"C":C,"KIND":KIND,"H":H,"CH":CH,"HE":HE,"HX":HX,"kek":kek,"corners":corners},open("benzene.json","w"))
