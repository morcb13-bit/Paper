# 金属結合の 3 つの型をシェルピンスキーの枠（辺 12、小立方体 6）で数える
from itertools import product, combinations
B=list(product([0,1],repeat=3)); odd=lambda n:sum(n)%2==1
def tet(o,e,w): return [tuple(o[i]+b[i]*e for i in range(3)) for b in B if odd(b)==w]
dot=lambda a,b:sum(x*y for x,y in zip(a,b)); sub=lambda a,b:tuple(a[i]-b[i] for i in range(3))
def grow(seed,nbr,K):
    gen={seed:0}; front=[seed]; sh=[1]
    for k in range(1,K+1):
        nxt=[]
        for p in front:
            for d in nbr(p):
                q=tuple(p[i]+d[i] for i in range(3))
                if q not in gen: gen[q]=k; nxt.append(q)
        sh.append(len(nxt)); front=nxt
    return gen,sh
print("== 面心立方 ==")
FCC=lambda p:((p[0]+p[1]+p[2])//6)%2==0 and all(x%6==0 for x in p)
N12=[d for d in product([-6,0,6],repeat=3) if dot(d,d)==72]
g,sh=grow((0,0,0),lambda p:N12,5); print("最近接",len(N12),"長さ² 72（正四面体の辺）","世代",sh)
SUB=[tet((6*n[0],6*n[1],6*n[2]),6,True) for n in B if odd(n)]; V=set(v for T in SUB for v in T)
print("シェルピンスキー2段目の頂点10 はすべて原子",all(FCC(v) for v in V),"／角8もすべて原子",all(FCC(c) for c in product([0,12],repeat=3)),"／中心(6,6,6)は",FCC((6,6,6)))
# 小立方体ごとの正四面体（偶の角 4 つ）
cubes=list(product(range(-2,4),repeat=3))
T={n:[c for c in [tuple(6*n[i]+6*b[i] for i in range(3)) for b in B] if FCC(c)] for n in cubes}
print("どの小立方体にも原子の正四面体が 1 個",set(len(v) for v in T.values()),"向き：奇の小立方体は",{('A' if sorted(T[n])==sorted(tet(tuple(6*x for x in n),6,True)) else 'B') for n in cubes if odd(n)},"偶は",{('A' if sorted(T[n])==sorted(tet(tuple(6*x for x in n),6,True)) else 'B') for n in cubes if not odd(n)})
a=(1,0,0);b=(0,0,0); print("隣の小立方体の正四面体どうし（奇と偶）の共有頂点",len(set(T[a])&set(T[b])),"→辺")
p=(0,0,0); print("原子 1 個が属する正四面体",sum(1 for n in cubes if p in T[n]),"正八面体",sum(1 for d in [(6,0,0),(-6,0,0),(0,6,0),(0,-6,0),(0,0,6),(0,0,-6)]))
print("単位格子（辺12）の体積 1728＝正四面体 8×72＋正八面体 4×288 ＝",8*72+4*288,"／原子",sum(1 for q in product(range(0,12,6),repeat=3) if FCC(q)))
print("\n== 体心立方 ==")
BCC=lambda p:all(x%12==0 for x in p) or all(x%12==6 for x in p)
N8=[d for d in product([-6,6],repeat=3)]
g,sh=grow((0,0,0),lambda p:N8,5); print("最近接",len(N8),"長さ²",dot(N8[0],N8[0]),"＝27×4（ダイヤモンドの結合の 2 倍の長さ）","世代",sh)
c=(6,6,6); dirs=[sub(q,c) for q in product([0,12],repeat=3)]
A=[d for d in dirs if sum(x>0 for x in d)%2==1]; Bd=[d for d in dirs if sum(x>0 for x in d)%2==0]
print("中心から 8 本：A の向き",A,"B の向き",Bd)
print("角 8 と中心は原子",all(BCC(q) for q in product([0,12],repeat=3)) and BCC(c),"面の中心（シェルピンスキーの共有点）は",[BCC(v) for v in V if v not in product([0,12],repeat=3)])
print("単位格子あたり原子",sum(1 for q in product(range(0,12,6),repeat=3) if BCC(q)))
print("2 番目に近い",sum(1 for d in product([-12,0,12],repeat=3) if dot(d,d)==144),"長さ² 144")
print("\n== 六方最密 ==")
e1=(6,-6,0); e2=(0,6,-6)
L=[tuple(i*e1[t]+j*e2[t] for t in range(3)) for i in range(-6,7) for j in range(-6,7)]
UPF=[(6,6,0),(6,0,6),(0,6,6)]          # 面心立方：B 層から上へ（C 層）
UPH=[(2,2,8),(2,8,2),(8,2,2)]          # 六方最密：B 層から上へ（A 層の真上）
print("上への 3 本の長さ²",[dot(d,d) for d in UPF],[dot(d,d) for d in UPH])
lay=lambda base:set(tuple(base[t]+l[t] for t in range(3)) for l in L)
A0=lay((0,0,0)); B1=lay((6,6,0)); A2=lay((8,8,8)); C2=lay((12,12,0))
print("六方最密の 3 層目は A の真上（(8,8,8) は (1,1,1) の向き）","面心立方の 3 層目 C は (12,12,0)")
b=(6,6,0)
print("B 層の原子の隣：同じ層",sum(1 for q in B1 if dot(sub(q,b),sub(q,b))==72),"下",sum(1 for q in A0 if dot(sub(q,b),sub(q,b))==72),"上（六方最密）",sum(1 for q in A2 if dot(sub(q,b),sub(q,b))==72),"上（面心立方）",sum(1 for q in C2 if dot(sub(q,b),sub(q,b))==72))
# 層をまたいで面を分け合う正四面体の組：A0 の三角形 + 上の B、上の B の三角形 + その上
def tets_between(lo,hi):
    out=[]
    for p in hi:
        if max(abs(x) for x in p)>30: continue
        nb=[q for q in lo if dot(sub(q,p),sub(q,p))==72]
        if len(nb)==3: out.append((frozenset(nb),p))
    return out
lowH=tets_between(A0,B1); upH=tets_between(B1,A2)
lowB=set(); 
# 下向きの正四面体：B 層の三角形 + 下の A 層の原子
def tets_down(hi_layer,lo_layer):
    res=[]
    for p in lo_layer:
        if max(abs(x) for x in p)>30: continue
        nb=[q for q in hi_layer if dot(sub(q,p),sub(q,p))==72]
        if len(nb)==3: res.append((frozenset(nb),p))
    return res
downAB=tets_down(B1,A0)   # B 層の三角形を面に、下の A 層に頂点
upBA=tets_between(B1,A2); upBC=tets_between(B1,C2)
fH=sum(1 for f,_ in downAB for g,_ in upBA if f==g); fF=sum(1 for f,_ in downAB for g,_ in upBC if f==g)
print("B 層の三角形を挟んで上下の正四面体が面を分け合う組：六方最密",fH,"面心立方",fF)
