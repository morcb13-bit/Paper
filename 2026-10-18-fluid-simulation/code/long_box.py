import numpy as np
from carrier3d_np import build, par, sgn
LX,LY=120,30
def inbox(x):
    ok=(sgn(np.array([LX,0])-x[:,0])>0)&(sgn(np.array([LX,0])+x[:,0])>0)
    for a in (1,2): ok&=(sgn(np.array([LY,0])-x[:,a])>0)&(sgn(np.array([LY,0])+x[:,a])>0)
    return ok
V,ties=build(inbox,maxit=2000)
x=par(V)
key=((V+64)*(128**np.arange(6))).sum(1); idx={k:i for i,k in enumerate(key.tolist())}
nb=-np.ones((len(V),12),np.int64)
for k in range(12):
    st=np.zeros(6,np.int64); st[k%6]=1 if k<6 else -1
    kk=(((V+st)+64)*(128**np.arange(6))).sum(1); nb[:,k]=[idx.get(q,-1) for q in kk.tolist()]
near=np.zeros((len(V),6),bool)
for a,L in ((0,LX),(1,LY),(2,LY)):
    near[:,2*a]=sgn(np.array([L-2,0])-x[:,a])<=0; near[:,2*a+1]=sgn(np.array([L-2,0])+x[:,a])<=0
np.savez_compressed('longbox.npz',addr=V,nb=nb,near=near,xpar=x)
print('verts',len(V),'ties',ties)
