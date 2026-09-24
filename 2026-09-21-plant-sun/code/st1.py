# 道具2（その一）：描き足した型紙の一桁の表を、全桁が同時に更新するオートマトンとして回し、落ち着くまでの段数を数える
#  一段＝全桁が一斉に、隣（一つ下の桁）の今の繰り上がりを読んで (s, c_out) を引き直す。初めの繰り上がりは全部 0
#  事前の基準：和が整数の和と一致（NG を返せる：表を一つずらすと落ちる）
#             段数を、無作為の入力と、繰り上がりが全桁を貫く入力の二通りで数える
import random
src=open('vp2.py').read().split('random.seed(6)')[0]
exec(src)
def todig_n(X,n):
    d=todig(X); return d+[0]*(n-len(d))
def run(xd,yd,T):
    n=len(xd); c=[0]*(n+1); s=[0]*n; t=0
    while True:
        new=[None]*(n+1); new[0]=0
        ns=[0]*n
        for i in range(n):
            v=T[i][(xd[i],yd[i],c[i])]; ns[i]=v[0]; new[i+1]=v[1]
        t+=1
        if new==c and ns==s: return t-1, sum(v*5**i for i,v in enumerate(s))+c[n]*5**n
        c,s=new,ns
random.seed(6)
for nd in (10,20,40):
    lim=(5**nd-1)//2; T=tabs[:nd]
    ok=0; steps=[]
    for _ in range(100):
        X,Y=random.randint(-lim,lim),random.randint(-lim,lim)
        k,v=run(todig_n(X,nd),todig_n(Y,nd),T); ok+=(v==X+Y); steps.append(k)
    xd=[2]+[1]*(nd-1); yd=[2]+[1]*(nd-1)        # 桁0で繰り上がり、以後 1+1+1=3 で運び続ける
    X=sum(v*5**i for i,v in enumerate(xd)); k,v=run(xd,yd,T)
    from collections import Counter
    print(f"{nd}桁  無作為: 和の一致 {ok}/100  段数 最大{max(steps)} 分布{dict(sorted(Counter(steps).items()))}")
    print(f"      貫く入力: 和の一致 {v==2*X}  段数 {k}")
# 対照：表の一つを壊すと和が落ちるか
bad=[dict(t) for t in tabs[:10]]; bad[3][(1,1,1)]=(0,0)
lim=(5**10-1)//2; random.seed(1)
print("対照 表を一か所壊す 10桁:",sum(run(todig_n(X,10),todig_n(Y,10),bad)[1]==X+Y for X,Y in [(random.randint(-lim,lim),random.randint(-lim,lim)) for _ in range(300)]),"/300")
