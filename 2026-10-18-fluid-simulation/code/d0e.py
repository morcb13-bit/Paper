import numpy as np, json
from tunnel3d import *
T=T3('box3d_g1.npz',base=4,wind=0)
NI=T.NI; V=T.V
# 縁からの深さ（内側の頂点だけを辺でたどる最短の歩数）
depth=np.full(V,-1); depth[NI:]=0; fr=np.arange(NI,V); d=0
while len(fr):
    d+=1; nb=T.nb[fr].reshape(-1); nb=nb[(nb>=0)]; nb=np.unique(nb); nb=nb[(nb<NI)]; nb=nb[depth[nb]<0]
    depth[nb]=d; fr=nb
deg=(T.nb[:NI]>=0).sum(1)
DW_=DW; DU_=DU; DP_=DP
acc=np.zeros((NI,2),np.int64); mass=np.zeros(NI,np.int64)
def pv(n):  # 3成分の運動量（P∥,T_u,T_v）をまとめた整数の組
    return np.concatenate(T.mom(n))
okb=True
for t in range(1300):
    n0=T.n[:NI].copy(); P0=pv(n0)
    pr=T.prep[T.step%5]
    # 縁の出入りの札の運動量
    flat=T.n.reshape(-1)
    T.collide()
    flat=T.n.reshape(-1); lab=np.tile(np.arange(12),V)
    def labmom(mask):
        v=flat[mask]; l=lab[mask]
        return np.concatenate([np.array([(v*DP_[l,0]).sum(),(v*DP_[l,1]).sum()]),np.array([(v*DU_[l,0]).sum(),(v*DU_[l,1]).sum()]),np.array([(v*DW_[l,0]).sum(),(v*DW_[l,1]).sum()])])
    Pin=labmom(pr['inn']); Pout=labmom(pr['out'])
    sums=np.add.reduceat(flat[pr['order']],pr['starts']); new=np.zeros(V*12,np.int64); new[pr['uq']]=sums; T.n=new.reshape(V,12)
    okb &= bool(np.array_equal(pv(T.n[:NI]),P0+Pin-Pout))
    T.n[NI:]=T.F; T.step+=1
    if t>=300:
        n=T.n[:NI]; acc[:,0]+=(n*DW_[:,0]).sum(1); acc[:,1]+=(n*DW_[:,1]).sum(1); mass+=n.sum(1)
np.savez('d0e_map.npz',acc=acc,mass=mass,depth=depth[:NI],deg=deg,band=T.band[:NI])
print('momentum ledger exact every step:',okb)
