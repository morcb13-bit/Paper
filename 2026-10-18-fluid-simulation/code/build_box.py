import numpy as np, sys
from carrier3d_np import *
e0=E[0]; u=np.array([O,O,I]); v=np.array([P,[-1,0],O])
Lf,Lt,M=int(sys.argv[1]),int(sys.argv[2]),3       # 箱の半分の長さ（流れ方向・横）と縁の帯
def region(x):
    a=dot(x,e0[None]); b=dot(x,u[None]); c=dot(x,v[None])
    def lt(z,L): return (sgn(np.array([L,0])-z)>0)&(sgn(np.array([L,0])+z)>0)
    return lt(a,2*(Lf+M))&lt(b,Lt+M)&lt(c,2*(Lt+M))
import time; t=time.time()
V,ties=build(region)
x=par(V); a=dot(x,e0[None]); b=dot(x,u[None]); c=dot(x,v[None])
def lt(z,L): return (sgn(np.array([L,0])-z)>0)&(sgn(np.array([L,0])+z)>0)
inner=lt(a,2*Lf)&lt(b,Lt)&lt(c,2*Lt)
# 帯の番号（流れ方向に3等分）: a を 2Lf の三分割で
band=(sgn(a-np.array([-2*Lf//3,0]))>0).astype(int)+(sgn(a-np.array([2*Lf//3,0]))>0).astype(int)
order=np.concatenate([np.where(inner)[0],np.where(~inner)[0]])
V=V[order]; inner=inner[order]; band=band[order]; x=x[order]
key=((V+64)*(128**np.arange(6))).sum(1); idx={k:i for i,k in enumerate(key.tolist())}
nb=-np.ones((len(V),12),np.int64)
for k in range(12):
    st=np.zeros(6,np.int64); st[k%6]=1 if k<6 else -1
    kk=(((V+st)+64)*(128**np.arange(6))).sum(1)
    nb[:,k]=[idx.get(q,-1) for q in kk.tolist()]
import os
np.savez_compressed(os.environ.get('OUT','box3d_%d_%d.npz'%(Lf,Lt)),addr=V,nb=nb,inner=inner,band=band,xpar=x)
deg=(nb[inner]>=0).sum(1)
print('verts',len(V),'inner',inner.sum(),'ties',ties,'time',round(time.time()-t,1),'deg(inner)',np.bincount(deg))
