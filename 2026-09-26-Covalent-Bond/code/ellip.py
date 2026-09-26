exec(open('h2walk.py').read().split("# 照合")[0])   # V（5倍座標）, VB, refl, dot, face などを読む
def sq(w): return dot(w,w)
def vsub(u,v): return tuple(sub(u[i],v[i]) for i in range(3))
O=(Z,Z,Z); OB=refl(O)
allv=list(V)+[v for v in VB if v not in V]
print("頂点",len(allv))
# 判定1：焦点までの距離の二乗
d1=set(sq(vsub(v,O)) for v in V); d2A=sorted(set(sq(vsub(v,OB)) for v in V))
print("A側の O_A までの距離二乗の種類:",len(d1))
print("A側の O_B までの距離二乗の種類:",len(d2A), d2A)
# 判定2：回転楕円体（中心＝共有面の中心 M=OB/2、軸=n）
# 2倍座標で M を整数に：P=2v, 2M=OB
nn=dot(n,n)
pts=set()
for v in allv:
    w=vsub(tuple(sc(2,c) for c in v) if False else tuple((2*c[0],2*c[1]) for c in v),OB)
    t=dot(w,n); X=mul(t,t); Y=sub(mul(nn,sq(w)),X)
    pts.add((X,Y))
pts=sorted(pts); print("(X,Y) の種類:",len(pts))
def collinear(P):
    if len(P)<3: return True
    (x0,y0),(x1,y1)=P[0],P[1]
    for (x,y) in P[2:]:
        cr=sub(mul(sub(x1,x0),sub(y,y0)),mul(sub(y1,y0),sub(x,x0)))
        if cr!=(0,0): return False
    return True
print("全頂点が一枚の回転楕円体に乗る:", "OK" if collinear(pts) else "NG")
# 負でない対照：正十二面体一個（中心 O、軸 n）
P1=set()
for v in V:
    t=dot(v,n); X=mul(t,t); Y=sub(mul(nn,sq(v)),X); P1.add((X,Y))
print("対照 正十二面体一個が回転楕円体（球）に乗る:", "OK" if collinear(sorted(P1)) else "NG", " (X,Y) の種類",len(P1))
