import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching
from tunnel3d import ADJ, ADJ2, OPP, LAB, DP
from carrier3d_np import dot
Z=np.load('box3d_g1.npz'); nb=Z['nb']; NI=int(Z['inner'].sum()); V=len(nb)
def cands(k, sets):
    """A の各点から届く候補（行き先頂点）の組を返す"""
    pairs=[]
    Aall=np.where(nb[:,k]<0)[0]
    if 'self' in sets: pairs.append((Aall,Aall))
    if 'r1' in sets:
        for o in ADJ[k]: w=nb[Aall,o]; m=w>=0; pairs.append((Aall[m],w[m]))
    if 'r2' in sets:
        for o in ADJ2[k]: w=nb[Aall,o]; m=w>=0; pairs.append((Aall[m],w[m]))
    if 'h2' in sets:   # 2歩：変位の札方向成分が正
        G=dot(LAB,LAB[k][None])   # 各向きの札kとの内積 (12×2)
        pos=lambda a,b: (G[a]+G[b])   # Z[φ]
        from carrier3d_np import sgn
        for a in range(12):
            for b in range(12):
                if b==OPP[a]: continue
                if sgn(pos(a,b)[None])[0]<=0: continue
                w=nb[Aall,a]; m=w>=0; w2=np.where(m,nb[np.maximum(w,0),b],-1); m2=w2>=0
                pairs.append((Aall[m2],w2[m2]))
    s=np.concatenate([p[0] for p in pairs]); t=np.concatenate([p[1] for p in pairs])
    return Aall,s,t
for sets in (['r1','self'],['r1','self','r2'],['r1','self','r2','h2']):
    k=0
    has_out=nb[:,k]>=0; indeg=np.zeros(V,int); np.add.at(indeg,nb[has_out,k],1)
    Bin=np.where((indeg==0)&(np.arange(V)<NI))[0]
    bpos=-np.ones(V,np.int64); bpos[Bin]=np.arange(len(Bin))
    Aall,s,t=cands(k,sets)
    apos=-np.ones(V,np.int64); apos[Aall]=np.arange(len(Aall))
    m=bpos[t]>=0
    M=csr_matrix((np.ones(m.sum(),np.int8),(apos[s[m]],bpos[t[m]])),shape=(len(Aall),len(Bin)))
    mm=maximum_bipartite_matching(M,perm_type='column')
    print(sets,'中の列の始まり',len(Bin),'覆えた',int((mm>=0).sum()))
