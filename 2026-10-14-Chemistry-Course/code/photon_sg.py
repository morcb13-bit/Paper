# 1〜3 段のシェルピンスキーを、銅原子 6 個の導線（3 段のシェルピンスキーを並べたもの）の上で +x へ転がす（整数だけ）
from itertools import product
mid=lambda a,b:tuple((a[i]+b[i])//2 for i in range(3))
split=lambda t:[[t[i]]+[mid(t[i],t[j]) for j in range(4) if j!=i] for i in range(4)]
def tet(o,s,par):   # 辺 s の立方体（原点 o）の頂点のうち、番地の和/s の偶奇が par のもの
    return [tuple(o[i]+s*c[i] for i in range(3)) for c in product((0,1),repeat=3) if (sum(o)//s+sum(c))%2==par]
def gasket(t,m):
    t=[t]
    for _ in range(m-1): t=[u for v in t for u in split(v)]
    return {p for u in t for p in u}
K=6; C=set()
for k in range(K): C|=gasket(tet((24*k-12,-12,-12),24,0),3)      # 導線の点
X0,X1=-12,24*K-12
for m in (1,2,3):
    s=6*2**(m-1)
    for par,name in ((0,"偶"),(1,"奇")):
        sit=best=tot=0
        for y in range(-12,12-s+1,s):
            for z in range(-12,12-s+1,s):
                for x0 in range(X0,X1-s+1,s):
                    tot+=1; sit+=gasket(tet((x0,y,z),s,par),m)<=C
                    run=0; x=x0
                    while x+s<=X1 and gasket(tet((x,y,z),s,par),m)<=C: run+=1; x+=s
                    best=max(best,run)
        if sit: print(f"{m} 段（辺 {s:2d}）{name}の族：導線に座れる升 {sit}／{tot}　続けて転がれた最長 {best}／{(X1-X0)//s}")
# 一辺 96 の 5 段のギャスケット（フォトン）の頂点が、銅の導線の点と重なる数（導線は原子 14 個）
C14=set()
for k in range(14): C14|=gasket(tet((24*k-12,-12,-12),24,0),3)
for off in (0,12,24,36):
    row=[]
    for k in range(3):
        c=(48+off+96*k,0,0)
        row.append(len(gasket(tet(tuple(x-48 for x in c),96,0),5)&C14))
    print(f"フォトンの位置のずれ {off:2d}：銅の頂点と重なる点（3 拍） {row}")
