# 道具2（その二）：繰り上がりを、共通の中心 C のまわりの相似で離れた桁へ一度に渡す
#  各桁は二つの記憶 A（直前の F(k) 桁ぶんの繰り上がりの関数）と B（その前の F(k−1) 桁ぶん）を持つ
#  一段：相手＝自分の型紙の目印を C のまわりに φ^(−3F(k)) 倍した所にある桁。A ← A∘(相手の B)、B ← 旧 A
#   相似の比は一段ごとに「前の二段の比の積」（φ^(3F(k+1)) = φ^(3F(k))·φ^(3F(k−1))）で、外から与えない
#   相手が見つからない（桁0より下）＝下から入る繰り上がり 0
#  関数は各桁の表 (x,y,c) → c_out の三つの値 (c=−1,0,+1)。合成は表を引き直すだけ
#  事前の基準：OK＝貫く入力で和が一致し、段数が桁数に依らないか log で伸びる。NG＝段数が 10・20・40
#   対照1：止める桁（x+y=0）を途中に挟んでも和が一致（繰り上がりがそこで止まる）
#   対照2：比を一つずらす（φ^(3F(k)+1)）と相手が見つからず和が落ちる
import random
exec(open('st1.py').read().split('random.seed(6)\nfor')[0])
PHI3i=U.zmul(U.zmul(U.zsub(U.PHI,U.ONE),U.zsub(U.PHI,U.ONE)),U.zsub(U.PHI,U.ONE))   # φ^−3（単数）
def upw(base,n):
    p=U.ONE
    for _ in range(n): p=U.zmul(p,base)
    return p
Z0=tuple(z0); D0=U.zsub(U.zadd(U.zadd(tuple(TPL['sa']),Z0),U.zadd(tuple(TPL['sa']),Z0)),C2)
def ref(i): return U.zadd(C2,U.zmul(upw(PHI3,i),D0))      # 桁 i の型紙の目印（スリット a、2倍の座標）
IDX={ref(i):i for i in range(41)}

def fastadd(xd,yd,T,wrong=False):
    n=len(xd)
    f=[tuple(T[i][(xd[i],yd[i],c)][1] for c in (-1,0,1)) for i in range(n)]
    A=f[:]; P=f[:]; a,b=1,1
    ra=rb=PHI3i                                  # φ^(−3a), φ^(−3b)
    steps=0; first_const=None
    while a<n:
        s=U.zmul(ra,U.zsub(U.PHI,U.ONE)) if wrong else ra
        newA=[]
        for i in range(n):
            j=IDX.get(U.zadd(C2,U.zmul(s,U.zsub(ref(i),C2))))
            q=P[j] if (j is not None and j<n) else (0,0,0)
            if j is None and i-a>=0: q=(-1,0,1)   # 相手が描けていない：何も渡らない（恒等）
            newA.append(tuple(A[i][q[c]+1] for c in range(3)))
        P=A; A=newA; a,b=a+b,a; ra,rb=U.zmul(ra,rb),ra; steps+=1
        if first_const is None and all(len(set(x))==1 for x in A): first_const=steps
    c=[0]+[A[i][1] for i in range(n)]
    s=[T[i][(xd[i],yd[i],c[i])][0] for i in range(n)]; steps+=1
    return steps, sum(v*5**i for i,v in enumerate(s))+c[n]*5**n, first_const
def val(d): return sum(v*5**i for i,v in enumerate(d))
random.seed(6)
for nd in (10,20,40):
    T=tabs[:nd]; lim=(5**nd-1)//2
    xd=[2]+[1]*(nd-1); k,v,fc=fastadd(xd,xd,T); kr,vr=run(xd,xd,T)
    print(f"{nd}桁 貫く入力: 和 {v==2*val(xd)}  段数 {k}（相似の渡し {k-1}＋読み出し 1）   全桁が同時の順回し {kr}")
    ok=0; ks=set()
    for _ in range(100):
        X,Y=random.randint(-lim,lim),random.randint(-lim,lim); kk,vv,_=fastadd(todig_n(X,nd),todig_n(Y,nd),T); ok+=vv==X+Y; ks.add(kk)
    print(f"      無作為100組: 和 {ok}/100  段数 {sorted(ks)}")
    yd=xd[:]; yd[nd//2]=-1                       # 真ん中の桁を x+y=0（止める）に
    k2,v2,_=fastadd(xd,yd,T); print(f"      対照1 真ん中に止める桁: 和 {v2==val(xd)+val(yd)}  段数 {k2}")
    k3,v3,_=fastadd(xd,xd,T,wrong=True); print(f"      対照2 比を φ だけずらす: 和 {v3==2*val(xd)}")
