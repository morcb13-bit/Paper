# 手が先端の小正四面体の辺を、先端 (c,-c,-c) から芯の正八面体の頂点 (0,-c,0) へ滑る
# 手 v1=(p,-c,-p) とその巡回 v2=(-p,p,-c), v3=(-c,-p,p)、p=c で先端、p=0 で芯の頂点
dot=lambda u,v:sum(a*b for a,b in zip(u,v))
c=6
for p in range(c,-1,-1):
    v1=(p,-c,-p); v2=(-p,p,-c); v3=(-c,-p,p)
    s=tuple(a+b+d for a,b,d in zip(v1,v2,v3))
    print(p, "組",(dot(v1,v2),dot(v1,v1)), "三つの和",s, "v2v3",dot(v2,v3),"v3v1",dot(v3,v1))
