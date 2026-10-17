# 黒鉛：ダイヤモンドの正四面体から 4 本目（(3,3,3) の向き）を外し、残る 3 本を面 x+y+z=一定 に寝かせる
from itertools import combinations
dot=lambda a,b:sum(x*y for x,y in zip(a,b))
DIA=[(3,-3,-3),(-3,3,-3),(-3,-3,3)]
# (1,1,1) の向きの成分を引く（3 倍してから引くと整数）：d − (d·w)/3 w,  d·w=−3 なので d+(1,1,1)
# (3,-3,-3)+(1,1,1)=(4,-2,-2)
GA=[tuple(d[i]+1 for i in range(3)) for d in DIA]
print("寝かせた 3 本",GA,"和",tuple(sum(v[i] for v in GA) for i in range(3)))
print("ダイヤモンド：長さ²",dot(DIA[0],DIA[0]),"内積",dot(DIA[0],DIA[1]),"／黒鉛：長さ²",dot(GA[0],GA[0]),"内積",dot(GA[0],GA[1]))
e=[tuple(a[i]-b[i] for i in range(3)) for a,b in combinations(GA,2)]; print("三角形の辺の長さ²",[dot(x,x) for x in e],"（正四面体の辺 72）")
GB=[tuple(-x for x in d) for d in GA]
key=lambda p:p
SEED=(0,0,0); gen={SEED:0}; side={SEED:"A"}; front=[SEED]; sh=[1]
for k in range(1,9):
    nxt=[]
    for p in front:
        for d in (GA if side[p]=="A" else GB):
            q=tuple(p[i]+d[i] for i in range(3))
            if q not in gen: gen[q]=k; side[q]="B" if side[p]=="A" else "A"; nxt.append(q)
    sh.append(len(nxt)); front=nxt
print("世代ごとの炭素",sh)
cA=[];cB=[];a=b=0
for k,n in enumerate(sh):
    if k%2==0: a+=n; cA.append(a)
    else: b+=n; cB.append(b)
print("A の三角形の累計（偶数世代）",cA,"／B の炭素の累計",cB)
allp=set(gen); print("内側の結合数",set(sum(1 for d in (GA if side[p]=="A" else GB) if tuple(p[i]+d[i] for i in range(3)) in allp) for p in allp if gen[p]<=6))
# 平面のシェルピンスキー：B の炭素の三角格子 B=i*e1+j*e2、A は上向きの三角形の中心
e1=(6,-6,0); e2=(0,6,-6)
Bs=set(p for p in gen if side[p]=="B"); As=set(p for p in gen if side[p]=="A")
def Bpt(i,j): return tuple(-GA[0][t]+i*e1[t]+j*e2[t] for t in range(3))
# 上向きの三角形 (i,j) の頂点と中心
def upc(i,j):
    V=[Bpt(i,j),Bpt(i+1,j),Bpt(i+1,j+1)]
    return V,tuple(sum(v[t] for v in V)//3 for t in range(3)) if all(sum(v[t] for v in V)%3==0 for t in range(3)) else None
V,c=upc(0,0); print("三角形 (0,0) の頂点",V,"中心",c,"A の炭素か",c in As)
for n in [2,4]:
    ups=[(i,j) for i in range(n) for j in range(n) if j<=i]
    cs=[upc(i,j)[1] for i,j in ups]
    print("辺",n,"の大きな三角形：上向きの小三角形",len(ups),"すべて A の炭素",all(x in As for x in cs))
print("平面のシェルピンスキー：2段目 3、3段目 9（辺 4 の中は 10＝9＋穴の中の 1）")
