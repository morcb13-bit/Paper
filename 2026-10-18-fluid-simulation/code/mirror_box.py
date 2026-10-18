import numpy as np, sys, os
from carrier3d_np import build, par, sgn, E
from tunnel3d import LAB, ADJ, OPP
L=int(sys.argv[1]); out=sys.argv[2]
def inbox(x,M):
    ok=np.ones(len(x),bool)
    for a in range(3): ok&=(sgn(np.array([L+M,0])-x[:,a])>0)&(sgn(np.array([L+M,0])+x[:,a])>0)
    return ok
V,ties=build(lambda x: inbox(x,0))
x=par(V)
key=((V+64)*(128**np.arange(6))).sum(1); idx={k:i for i,k in enumerate(key.tolist())}
nb=-np.ones((len(V),12),np.int64)
for k in range(12):
    st=np.zeros(6,np.int64); st[k%6]=1 if k<6 else -1
    kk=(((V+st)+64)*(128**np.arange(6))).sum(1); nb[:,k]=[idx.get(q,-1) for q in kk.tolist()]
# 壁に近い（2未満）か：6面
near=np.zeros((len(V),6),bool)
for a in range(3):
    near[:,2*a]=sgn(np.array([L-2,0])-x[:,a])<=0      # +側の壁
    near[:,2*a+1]=sgn(np.array([L-2,0])+x[:,a])<=0    # −側の壁
np.savez_compressed(out,addr=V,nb=nb,near=near,xpar=x)
print('verts',len(V),'ties',ties,'near-wall',near.any(1).sum(),'deg',np.bincount((nb>=0).sum(1)))
