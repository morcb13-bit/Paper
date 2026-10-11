# 正三角形 n 枚を 120/n 度ずつ回して同じ中心に重ねたときの数（交点・辺・面）
# 注意：40度などは整数の格子に乗らないので、ここだけ50桁の数で同一判定する（観察用）
from mpmath import mp, mpf, cos, sin, pi
mp.dps=50
EPS=mpf(10)**-35
def count(n):
    segs=[]
    for k in range(n):
        th=2*pi*k/(3*n)
        P=[(cos(th+2*pi*j/3+pi/2),sin(th+2*pi*j/3+pi/2)) for j in range(3)]
        segs+= [(P[j],P[(j+1)%3]) for j in range(3)]
    pts=[]
    def key(p):
        for i,q in enumerate(pts):
            if abs(p[0]-q[0])<EPS and abs(p[1]-q[1])<EPS: return i
        pts.append(p); return len(pts)-1
    on=[set() for _ in segs]
    for i,(a,b) in enumerate(segs):
        on[i]|={key(a),key(b)}
    for i in range(len(segs)):
        for j in range(i+1,len(segs)):
            (a,b),(c,d)=segs[i],segs[j]
            r=(b[0]-a[0],b[1]-a[1]); s=(d[0]-c[0],d[1]-c[1])
            den=r[0]*s[1]-r[1]*s[0]
            if abs(den)<EPS: continue
            t=((c[0]-a[0])*s[1]-(c[1]-a[1])*s[0])/den
            u=((c[0]-a[0])*r[1]-(c[1]-a[1])*r[0])/den
            if -EPS<=t<=1+EPS and -EPS<=u<=1+EPS:
                k=key((a[0]+t*r[0],a[1]+t*r[1])); on[i].add(k); on[j].add(k)
    V=len(pts); E=sum(len(o)-1 for o in on); F=E-V+2-1   # 外側を除いた面
    tips=3*n
    return V,E,F,tips,V-tips
print('枚数 角度 点 辺 面(内側) 先端 交点')
for n in range(1,6):
    V,E,F,t,x=count(n); print(n, f'{120/n:g}度', V,E,F,t,x)
