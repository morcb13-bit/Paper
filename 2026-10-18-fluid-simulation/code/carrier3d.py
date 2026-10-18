# 立体の担体（6次元→3次元の切断投影）。番地＝6次元の整数。採否は Z[φ] の整数比較だけで決める。
import json, sys
from collections import deque
def zadd(x,y): return (x[0]+y[0],x[1]+y[1])
def zsub(x,y): return (x[0]-y[0],x[1]-y[1])
def zmul(x,y): a,b=x;c,d=y; return (a*c+b*d, a*d+b*c+b*d)   # φ²=φ+1
def zint(n,x): return (n*x[0],n*x[1])
def conj(x): return (x[0]+x[1],-x[1])                        # φ → 1−φ
def sgn(x):
    p,q=x; A=2*p+q; B=q
    if A>=0 and B>=0: return 0 if A==0 and B==0 else 1
    if A<=0 and B<=0: return -1
    return (1 if A*A>5*B*B else -1 if A*A<5*B*B else 0) if A>0 else (1 if 5*B*B>A*A else -1 if 5*B*B<A*A else 0)
def zabs(x): return x if sgn(x)>=0 else zint(-1,x)
O=(0,0); I=(1,0); P=(0,1)
E=[(I,P,O),(zint(-1,I),P,O),(O,I,P),(O,zint(-1,I),P),(P,O,I),(zint(-1,P),O,I)]   # 五回軸6本
Ep=[tuple(conj(c) for c in e) for e in E]
def vadd(u,v): return tuple(zadd(a,b) for a,b in zip(u,v))
def dot(u,v):
    s=O
    for a,b in zip(u,v): s=zadd(s,zmul(a,b))
    return s
def cross(u,v): return (zsub(zmul(u[1],v[2]),zmul(u[2],v[1])), zsub(zmul(u[2],v[0]),zmul(u[0],v[2])), zsub(zmul(u[0],v[1]),zmul(u[1],v[0])))
def proj(n,B):
    x=(O,O,O)
    for k in range(6): x=vadd(x,tuple(zint(n[k],c) for c in B[k]))
    return x
# 窓（ひし形30面体）の15面
D=37; g=(1,2,3,5,8,13)        # 窓のずらし γ=g/D（境界に格子点が乗らない一般の位置）
faces=[]
for i in range(6):
    for j in range(i+1,6):
        N=cross(Ep[i],Ep[j]); h=O; c=O
        for k in range(6):
            d=dot(N,Ep[k]); h=zadd(h,zabs(d)); c=zadd(c,zint(g[k],d))
        faces.append((N,zint(D,h),zint(2,c)))
ties=0
def accept(n):
    global ties
    x=proj(n,Ep)
    for N,h,c in faces:
        t=zsub(zint(2*D,dot(N,x)),c)
        s=sgn(zsub(h,zabs(t)))
        if s==0: ties+=1
        if s<=0: return False
    return True
R2=zint(int(sys.argv[1])**2 if len(sys.argv)>1 else 9, (2,1))    # 半径²（辺の長さ² = 2+φ の倍数）
def inside(n): x=proj(n,E); return sgn(zsub(R2,dot(x,x)))>=0
start=(0,)*6
assert accept(start), "原点が窓の外"
seen={start}; q=deque([start]); verts=[]
while q:
    n=q.popleft(); verts.append(n)
    for k in range(6):
        for s in (1,-1):
            m=list(n); m[k]+=s; m=tuple(m)
            if m not in seen and inside(m) and accept(m): seen.add(m); q.append(m)
idx={v:i for i,v in enumerate(verts)}
edges=[]; deg=[0]*len(verts)
for v in verts:
    for k in range(6):
        m=list(v); m[k]+=1; m=tuple(m)
        if m in idx: edges.append((idx[v],idx[m],k)); deg[idx[v]]+=1; deg[idx[m]]+=1
phi=(1+5**0.5)/2
def f(x): return x[0]+x[1]*phi
pos=[[round(f(c),4) for c in proj(v,E)] for v in verts]
hist={}
for d in deg: hist[d]=hist.get(d,0)+1
print('vertices',len(verts),'edges',len(edges),'ties',ties,'deg hist',sorted(hist.items()))
# 中心の頂点（原点）の渡れる向き
dirs=[]
for k in range(6):
    for s in (1,-1):
        m=list(start); m[k]+=s; dirs.append([k,s,tuple(m) in idx])
json.dump(dict(pos=pos,edges=edges,deg=deg,center=idx[start],dirs=dirs,addr=[list(v) for v in verts],hist=hist),open('carrier3d.json','w'))
