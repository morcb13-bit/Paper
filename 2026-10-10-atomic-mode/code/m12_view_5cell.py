# 正5胞体（4次元目に支点）を固定し、3次元への見かけだけを変える
# 頂点：H 三本 v1,v2,v3、埋まった席 L、支点 F（整数座標 n=4、内積 -5、長さ² 20）
# 見る向き w = a·L + b·F を潰すと、H どうしの組は (-(W+5s²), 4W-5s²)、W=20a²-10ab+20b²、s=a+b
import math
def pair(a,b):
    W=20*a*a-10*a*b+20*b*b; s=a+b
    p,q=-(W+5*s*s),4*W-5*s*s; g=math.gcd(abs(p),q); return p//g,q//g
# 整数の確認（座標を実際に作って潰す）
d=5; V=[[(d if k==i else 0)-1 for k in range(d)] for i in range(d)]
H,L,F=V[:3],V[3],V[4]
dot=lambda u,v:sum(x*y for x,y in zip(u,v))
for a,b in [(1,0),(0,1),(1,-1),(1,1),(2,1),(3,1)]:
    w=[a*x+b*y for x,y in zip(L,F)]; W=dot(w,w)
    proj=lambda v:[W*x-dot(v,w)*y for x,y in zip(v,w)]   # W 倍して整数のまま
    P=[proj(h) for h in H]; p,q=dot(P[0],P[1]),dot(P[0],P[0]); g=math.gcd(abs(p),q)
    print(f"w={a}L+{b}F  座標で {(p//g,q//g)}  式で {pair(a,b)}  角 {math.degrees(math.acos(p/q)):.4f}")
print()
for nm,deg in [("NH3",106.769),("H3O+",111.95)]:
    best=[]
    for a in range(0,40):
        for b in range(-40,41):
            if (a,b)==(0,0) or math.gcd(a,abs(b))!=1: continue
            p,q=pair(a,b); ang=math.degrees(math.acos(p/q))
            if abs(ang-deg)<=0.02: best.append((abs(a)+abs(b),a,b,(p,q),round(ang,4)))
    print(nm, sorted(best)[:4])
