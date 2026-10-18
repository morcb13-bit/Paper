# 立体の担体（番地表）を配列でまとめて作る。Z[φ] の数は (a,b)=a+bφ の整数組。
import numpy as np
def zm(x,y):  # (…,2)×(…,2)
    a,b=x[...,0],x[...,1]; c,d=y[...,0],y[...,1]
    return np.stack([a*c+b*d, a*d+b*c+b*d],-1)
def conj(x): return np.stack([x[...,0]+x[...,1],-x[...,1]],-1)
def sgn(x):
    p,q=x[...,0],x[...,1]; A=2*p+q; B=q
    s=np.zeros(A.shape,np.int64)
    s=np.where((A>=0)&(B>=0)&((A!=0)|(B!=0)),1,s)
    s=np.where((A<=0)&(B<=0)&((A!=0)|(B!=0)),-1,s)
    s=np.where((A>0)&(B<0),np.sign(A*A-5*B*B),s)
    s=np.where((A<0)&(B>0),np.sign(5*B*B-A*A),s)
    return s
def zabs(x): return x*np.where(sgn(x)<0,-1,1)[...,None]
O,I,P=[0,0],[1,0],[0,1]
E=np.array([[I,P,O],[[-1,0],P,O],[O,I,P],[O,[-1,0],P],[P,O,I],[[0,-1],O,I]],np.int64)  # 6×3×2
Ep=conj(E)
def dot(u,v): return zm(u,v).sum(-2)
def cross(u,v):
    return np.stack([zm(u[...,1,:],v[...,2,:])-zm(u[...,2,:],v[...,1,:]),
                     zm(u[...,2,:],v[...,0,:])-zm(u[...,0,:],v[...,2,:]),
                     zm(u[...,0,:],v[...,1,:])-zm(u[...,1,:],v[...,0,:])],-2)
import os
D=int(os.environ.get('WIN_D','37')); g=np.array([int(x) for x in os.environ.get('WIN_G','1,2,3,5,8,13').split(',')])
F=[]; 
for i in range(6):
    for j in range(i+1,6): F.append(cross(Ep[i],Ep[j]))
F=np.array(F)                                    # 15×3×2
A=np.stack([dot(F,Ep[k][None]) for k in range(6)],1)   # 15×6×2   N_f·Ep_k
H=zabs(A).sum(1)                                  # 15×2
C=(A*g[None,:,None]).sum(1)                       # 15×2
def accept(n):                                    # n: M×6
    s=(n[:,None,:,None]*A[None]).sum(2)           # M×15×2
    t=2*D*s-2*C[None]
    return np.all(sgn(D*H[None]-zabs(t))>0,1), np.any(sgn(D*H[None]-zabs(t))==0)
def par(n): return (n[:,:,None,None]*E[None]).sum(1)   # M×3×2
def build(region, maxit=400):
    """region(xpar)->bool。原点から ±単位の一歩で広げ、窓と region の両方を満たす点を集める"""
    def key(n): return ((n+64)*(128**np.arange(6))).sum(1)
    start=np.zeros((1,6),np.int64)
    ok,_=accept(start); assert ok[0]
    allv=[start]; seen=set(key(start).tolist()); fr=start; ties=0
    steps=np.concatenate([np.eye(6,dtype=np.int64),-np.eye(6,dtype=np.int64)])
    for it in range(maxit):
        cand=(fr[:,None,:]+steps[None]).reshape(-1,6)
        k=key(cand); _,ui=np.unique(k,return_index=True); cand=cand[ui]; k=k[ui]
        new=np.array([x not in seen for x in k.tolist()],bool)
        cand=cand[new]; k=k[new]
        if len(cand)==0: break
        a,tie=accept(cand); ties+=int(tie)
        r=region(par(cand)); m=a&r
        seen.update(k.tolist())        # 調べた点は二度調べない
        fr=cand[m]
        if len(fr)==0: break
        allv.append(fr)
    V=np.concatenate(allv)
    return V, ties
