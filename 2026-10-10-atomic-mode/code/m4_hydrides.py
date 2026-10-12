# L殻・M殻の水素化合物：席・手・手の組（整数）と実測の組の照合（検定2をM殻へ広げたもの）
# 判定は整数の交差掛け算だけ。浮動小数は実測角を組に換算する一回だけ。
import math
from itertools import combinations
A=[(1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1)]; B=[tuple(-x for x in v) for v in A]
dot=lambda u,v:sum(a*b for a,b in zip(u,v))
same=lambda a,b:a[0]*b[1]==b[0]*a[1]
def hands(e):
    fa=min(e,4); fb=max(0,e-4)
    if fa<4: return A[:fa]          # A 完成前：埋まった A の席が手
    return B[fb:]                   # A 完成後：B の空き席が手
def pairs(hs): return sorted({(dot(u,v),dot(u,u)) for u,v in combinations(hs,2)})
def to_pair(deg,tol=0.1,qmax=40):
    lo,hi=sorted([math.cos(math.radians(deg-tol)),math.cos(math.radians(min(deg+tol,180)))])
    for q in range(1,qmax+1):
        for p in range(-q,q+1):
            if lo<=p/q<=hi: return (p,q)
# 実測の角（[知識] 気相・単量体、CCCBDB 等）。角のないものは None
DATA=[("L","LiH",1,None),("L","BeH2",2,180.0),("L","BH3",3,120.0),("L","CH4",4,109.47),
      ("L","NH3",5,106.67),("L","H2O",6,104.48),("L","HF",7,None),("L","Ne",8,None),
      ("M","NaH",1,None),("M","MgH2",2,180.0),("M","AlH3",3,120.0),("M","SiH4",4,109.47),
      ("M","PH3",5,93.35),("M","H2S",6,92.11),("M","HCl",7,None),("M","Ar",8,None)]
print("殻 分子   電子 手  模型の組        実測角   実測の組  判定")
for sh,nm,e,deg in DATA:
    hs=hands(e); mp=pairs(hs)
    if deg is None or len(hs)<2:
        print(f"{sh}  {nm:5} {e:3} {len(hs):2}  {str(mp):15} {'―':>7}")
        continue
    m=to_pair(deg); ok=all(same(x,m) for x in mp)
    print(f"{sh}  {nm:5} {e:3} {len(hs):2}  {str(mp):15} {deg:7}  {str(m):9} {'OK' if ok else 'NG'}")
print("\nL : M の大きさ：フィボナッチの組 F3:F4 =",2,":",3," その2乗 =",4,":",9," ボーアの n²（n=2,3）=",4,":",9," 一致:",2*2==4 and 3*3==9)
