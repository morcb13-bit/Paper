# 星型八面体と立方体の往復：辺の中点 e を (2N,2N,0) から (N,N,0) へ、6×体積を整数で
from itertools import product
N=16
C=[(x,y,z) for x in(-2*N,2*N) for y in(-2*N,2*N) for z in(-2*N,2*N)]
def tris(j):
    T=[]
    for a in range(3):
        b,c=[i for i in range(3) if i!=a]
        for sb,sc in product((-1,1),repeat=2):
            e=[0,0,0];e[b]=sb*(2*N-j);e[c]=sc*(2*N-j)
            f1=[0,0,0];f1[b]=sb*2*N
            f2=[0,0,0];f2[c]=sc*2*N
            for sa in(-1,1):
                cc=[0,0,0];cc[a]=sa*2*N;cc[b]=sb*2*N;cc[c]=sc*2*N
                T+= [(cc,e,f1),(cc,e,f2)]
    return T
det=lambda p,q,r:p[0]*(q[1]*r[2]-q[2]*r[1])-p[1]*(q[0]*r[2]-q[2]*r[0])+p[2]*(q[0]*r[1]-q[1]*r[0])
for j in (0,1,N):
    T=tris(j); print(j,len(T),sum(abs(det(*t)) for t in T), 6*64*N**3)
