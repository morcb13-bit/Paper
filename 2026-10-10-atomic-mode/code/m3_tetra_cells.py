# 殻は正四面体の複合体：整数だけの検算
from itertools import combinations
def sub(u,v): return tuple(a-b for a,b in zip(u,v))
def det(a,b,c): return a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])
def v6(t): p=t[0]; return abs(det(sub(t[1],p),sub(t[2],p),sub(t[3],p)))   # 6×体積
dot=lambda u,v: sum(a*b for a,b in zip(u,v))

print("== 1. 星型八面体 = 芯の正八面体 + 小正四面体8（座標は2倍して整数） ==")
A=[(2,2,2),(2,-2,-2),(-2,2,-2),(-2,-2,2)]; B=[tuple(-x for x in v) for v in A]
core=[(2,0,0),(-2,0,0),(0,2,0),(0,-2,0),(0,0,2),(0,0,-2)]
# 芯の八面体の6×体積：中心から8つの面への四面体の和
faces=[(s1,s2,s3) for s1 in (core[0],core[1]) for s2 in (core[2],core[3]) for s3 in (core[4],core[5])]
octa=sum(v6(((0,0,0),)+f) for f in faces)
def spike(tip):  # 先端 tip と、tip に向いた芯の面の3点
    f=[c for c in core if dot(c,tip)>0]; return [tip]+f
small=[v6(spike(t)) for t in A+B]
bigA=v6(A); cube=4*4*4*6
print(" 大きな正四面体 A:",bigA," 芯の正八面体:",octa," 小正四面体8個:",small)
print(" A = 芯 + 先端4 :",bigA,"=",octa,"+",4*small[0],"->",bigA==octa+4*small[0])
print(" 星 = 芯 + 先端8 :",octa+8*small[0]," 立方体:",cube," 星×2=立方体 ->",2*(octa+8*small[0])==cube)
print(" 先端8個の頂点 = 立方体の8頂点 = L殻の8席 ->",sorted(A+B)==sorted((x,y,z) for x in(2,-2) for y in(2,-2) for z in(2,-2)))
print(" 小正四面体の辺²:",{dot(sub(a,b),sub(a,b)) for a,b in combinations(spike(A[0]),2)}," 大きな正四面体の辺²:",{dot(sub(a,b),sub(a,b)) for a,b in combinations(A,2)})

print("\n== 2. 席の埋まり方（A から埋め、B の空き席が手） ==")
for name,e in [("C",4),("N",5),("O",6),("F",7),("Ne",8)]:
    fa=min(e,4); fb=e-fa; h=4-fb
    print(f" {name:2} 電子{e} A{fa} B{fb} 手{h}  B側={'そろう' if fb in(0,4) else '混ざる'}")

print("\n== 3. 手の組：模型 (-1,3) と実測（m2_hand_angles.py の換算） ==")
print(" 分子  実測(p,q)  交差の差  3p+q  そろえた長さ²L=3q  内積(模型→実測)  底辺²(模型→実測) 差  単位8の何個分  内積-2にそろえた長さ²")
for name,(p,q) in [("CH4",(-1,3)),("NH3",(-2,7)),("H2O",(-1,4))]:
    L=3*q//(3 if q==3 else 1) if q!=3 else 3
    dm=-L//3; dr=p*L//q
    em=2*L-2*dm; er=2*L-2*dr
    cross=p*3-(-1)*q
    k=-2*q//p  # 内積を-2にそろえた長さ²
    print(f" {name:4} ({p},{q})   {cross:3}      {3*p+q:3}   {L:3}   {dm}→{dr}   {em}→{er}  {er-em:3}   {em//8}   {k}")
