# 手を整数の組 (a,b,b)（とその巡回）で書く：手どうしの内積 b(2a+b)、長さ² a²+2b²
# (1,-1,-1) 正四面体の手、(2,-1,-1) 平面 120 度、(1,0,0)/(1,-2,-2) 直角
# 実測の角 ±0.1 度に入る最小の (a,b) を探す（探索のための換算だけ浮動小数）
import math
DATA=[("BeH2",180.0),("BH3",120.0),("CH4",109.47),("NH3",106.67),("H2O",104.48),
      ("PH3",93.35),("H2S",92.11),("AsH3",92.08),("H2Se",90.90),("SbH3",91.56),("H2Te",90.26)]
def found(deg,tol=0.1,N=200):
    best=None
    for a in range(0,N):
        for b in range(-N,1):
            if math.gcd(a,abs(b))!=1: continue
            L=a*a+2*b*b
            if L==0: continue
            c=b*(2*a+b)/L
            if abs(math.degrees(math.acos(max(-1,min(1,c))))-deg)<=tol:
                k=(abs(a)+abs(b))
                if best is None or k<best[0]: best=(k,a,b,b*(2*a+b),L)
    return best
for nm,deg in DATA:
    r=found(deg)
    print(nm,deg, "(a,b)=",(r[1],r[2]) if r else None, "組=",(r[3],r[4]) if r else None, " b/a 側:", "直角(0)寄り" if r and abs(r[2])<r[1]/2 else "直角(-2)寄り" if r else "")
