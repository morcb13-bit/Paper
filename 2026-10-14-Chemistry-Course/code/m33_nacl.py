# NaCl（岩塩型）：辺 6 の格子点に Na と Cl が番地の偶奇で交互に並ぶ。シェルピンスキーとの対応
from itertools import product, combinations
B=list(product([0,1],repeat=3)); odd=lambda n:sum(n)%2==1
def tet(o,e,w): return [tuple(o[i]+b[i]*e for i in range(3)) for b in B if odd(b)==w]
ion=lambda p:"Na" if (sum(p)//6)%2==0 else "Cl"     # 番地 p/6 の和が偶数なら Na
SUB=[tet((6*n[0],6*n[1],6*n[2]),6,True) for n in B if odd(n)]
V=set(v for T in SUB for v in T)
print("シェルピンスキー2段目の頂点10 はすべて Na か",all(ion(v)=="Na" for v in V))
OCT=[(12,6,6),(0,6,6),(6,12,6),(6,0,6),(6,6,12),(6,6,0)]
print("真ん中の穴（正八面体）の頂点＝共有点6",sorted(OCT)==sorted(v for v in V if sorted(x in(0,12) for x in v).count(True)==1),"中心 (6,6,6) は",ion((6,6,6)))
print("小さな正四面体の中心 (9,3,3) 型は格子点でない（空き）")
# 小立方体 1 個の 8 角：Na の正四面体と Cl の正四面体
for n in [(1,0,0),(0,0,0)]:
    o=tuple(6*x for x in n); C=[tuple(o[i]+b[i]*6 for i in range(3)) for b in B]
    print(" 小立方体",n,"Na",[c for c in C if ion(c)=="Na"],"→正四面体の向き",'A' if sorted(c for c in C if ion(c)=='Na')==sorted(tet(o,6,True)) else 'B')
# 結合（最近接 Na–Cl）は立方体の辺、正四面体の辺（面の対角線）は同じイオンどうし
d2=lambda p,q:sum((p[i]-q[i])**2 for i in range(3))
print("最近接 長さ²",36,"＝立方体の辺。正四面体の辺 72 は Na–Na／Cl–Cl")
# 育つ：Cl 1 個（正八面体の中心）から最近接をたどる
seed=(6,6,6); gen={seed:0}; front=[seed]; sh=[1]
D=[(6,0,0),(-6,0,0),(0,6,0),(0,-6,0),(0,0,6),(0,0,-6)]
for k in range(1,6):
    nxt=[]
    for p in front:
        for d in D:
            q=tuple(p[i]+d[i] for i in range(3))
            if q not in gen: gen[q]=k; nxt.append(q)
    sh.append(len(nxt)); front=nxt
print("世代ごとのイオン",sh,"偶数世代は",set(ion(p) for p in gen if gen[p]%2==0),"奇数世代は",set(ion(p) for p in gen if gen[p]%2==1))
cum=[];c=0
for k in range(0,6,2): c+=sh[k]; cum.append(c)
print("正八面体（Cl が中心）の累計",cum)
# 正八面体どうしの共有：隣の Cl (12,12,6)
o1=set(tuple(seed[i]+d[i] for i in range(3)) for d in D); s2=(12,12,6); o2=set(tuple(s2[i]+d[i] for i in range(3)) for d in D)
print("隣の正八面体と共有する頂点",sorted(o1&o2),"→ 辺を 1 本共有")
# 単位格子（辺 12）あたり
cell=[p for p in product(range(0,12,6),repeat=3)]   # 片側だけ含む格子 [0,12) で数える（整数）
print("単位格子あたり Na",sum(1 for p in cell if ion(p)=="Na"),"Cl",sum(1 for p in cell if ion(p)=="Cl"))
print("\n== 頂点に集まる正四面体の数 ==")
# シェルピンスキー2段目の内側の共有点
print("シェルピンスキー2段目：共有点ごと",set(sum(1 for T in SUB if v in T) for v in V if sum(1 for T in SUB if v in T)>1))
# ダイヤモンド：A の正四面体は奇の小立方体すべて
DT=[tet((6*n[0],6*n[1],6*n[2]),6,True) for n in product(range(-2,4),repeat=3) if odd(n)]
print("ダイヤモンド：頂点 (12,6,6) に集まる A の正四面体",sum(1 for T in DT if (12,6,6) in T))
GA=[(4,-2,-2),(-2,4,-2),(-2,-2,4)]
print("黒鉛：B の炭素に集まる三角形 3（結合 3 本の先がそれぞれ三角形の中心）")
