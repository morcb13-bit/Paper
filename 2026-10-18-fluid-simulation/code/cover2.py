# 札kの「列の始まり（箱の中）」をすべて、届く範囲の「列の終わり」で一対一に埋められるか。歩数の少ない組を優先。
import numpy as np, sys
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import min_weight_full_bipartite_matching, maximum_bipartite_matching
from tunnel3d import ADJ, ADJ2, OPP, LAB
from carrier3d_np import dot, sgn
fn=sys.argv[1]; Z=np.load(fn); nb=Z['nb']; NI=int(Z['inner'].sum()); V=len(nb)
G=dot(LAB[:,None],LAB[None,:]) if False else None
def reach(Aall,k,maxh):
    """A から maxh 歩以内で、変位の札k方向成分が正（または0歩）の行き先と歩数"""
    out=[(Aall,Aall,np.zeros(len(Aall),int))]
    Dk=dot(LAB,LAB[k][None])            # 各向きと札kの内積（Z[φ]）
    cur=[(Aall,Aall,np.zeros((len(Aall),2),np.int64))]   # (起点, いまの頂点, 変位の札k成分)
    for h in range(1,maxh+1):
        nxt=[]
        for s,v,p in cur:
            for o in range(12):
                w=nb[v,o]; m=w>=0
                q=p[m]+Dk[o]
                nxt.append((s[m],w[m],q))
        s=np.concatenate([x[0] for x in nxt]); v=np.concatenate([x[1] for x in nxt]); p=np.concatenate([x[2] for x in nxt])
        key=s*V+v; _,ui=np.unique(key,return_index=True); s,v,p=s[ui],v[ui],p[ui]
        pos=sgn(p)>0
        out.append((s[pos],v[pos],np.full(pos.sum(),h)))
        cur=[(s,v,p)]
    s=np.concatenate([x[0] for x in out]); t=np.concatenate([x[1] for x in out]); h=np.concatenate([x[2] for x in out])
    key=s*V+t; o=np.lexsort((h,key)); key,s,t,h=key[o],s[o],t[o],h[o]
    first=np.r_[True,key[1:]!=key[:-1]]
    return s[first],t[first],h[first]
k=int(sys.argv[2]); maxh=int(sys.argv[3])
has_out=nb[:,k]>=0; indeg=np.zeros(V,int); np.add.at(indeg,nb[has_out,k],1)
Aall=np.where(~has_out)[0]; Bin=np.where((indeg==0)&(np.arange(V)<NI))[0]
bpos=-np.ones(V,np.int64); bpos[Bin]=np.arange(len(Bin)); apos=-np.ones(V,np.int64); apos[Aall]=np.arange(len(Aall))
s,t,h=reach(Aall,k,maxh); m=bpos[t]>=0; s,t,h=s[m],t[m],h[m]
M=csr_matrix((np.ones(len(s),np.int32),(bpos[t],apos[s])),shape=(len(Bin),len(Aall)))
mm=maximum_bipartite_matching(M,perm_type='column')
print('札',k,'歩数≤',maxh,'列の始まり(中)',len(Bin),'覆えた',int((mm>=0).sum()))
# 歩数の少ない組を優先した一対一（列の始まりをすべて埋める）
W=csr_matrix((h.astype(float)+0.0,(bpos[t],apos[s])),shape=(len(Bin),len(Aall)))
r,c=min_weight_full_bipartite_matching(W)
used=np.zeros(len(Aall),bool); used[c]=True
hh=np.zeros(len(Bin),int)
lk={(a,b):x for a,b,x in zip(apos[s].tolist(),bpos[t].tolist(),h.tolist())}
hh=np.array([lk[(ci,ri)] for ri,ci in zip(r.tolist(),c.tolist())])
print('歩数の内訳（0=その場,1,2,3）',np.bincount(hh,minlength=4))
Aint_un=Aall[(~used)&(Aall<NI)]
print('相手のいない列の終わり(中)',len(Aint_un))
# 枠まで何歩か
depth=np.full(V,-1); depth[NI:]=0; fr=np.arange(NI,V); d=0
while len(fr):
    d+=1; x=nb[fr].reshape(-1); x=np.unique(x[x>=0]); x=x[(x<NI)]; x=x[depth[x]<0]; depth[x]=d; fr=x
print('  その深さの内訳',np.bincount(depth[Aint_un]))
