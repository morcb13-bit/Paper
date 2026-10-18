# 一対一の渡りの表を作る
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching
from tunnel3d import ADJ, ADJ2, OPP
def build_bij(nb, NI):
    V=len(nb); dest=np.full((V,12),-1,np.int64); report=[]
    for k in range(12):
        has_out=nb[:,k]>=0
        d=np.where(has_out,nb[:,k],-1)
        indeg=np.zeros(V,np.int64); np.add.at(indeg,d[has_out],1)
        A=np.where(~has_out)[0]                 # 列の終わり（出口なし）
        B=np.where((indeg==0)&(np.arange(V)<NI))[0]   # 箱の中の列の始まり（入りなし）
        bpos=-np.ones(V,np.int64); bpos[B]=np.arange(len(B))
        match=-np.ones(len(A),np.int64)
        for ring in (list(ADJ[k]),['self'],list(ADJ2[k]),list(ADJ[k])+['self']+list(ADJ2[k])):
            free=np.where(match<0)[0]
            taken=np.zeros(len(B),bool); taken[match[match>=0]]=True
            rows=[];cols=[]
            for o in ring:
                w=A[free] if o=='self' else nb[A[free],o]; ok=(w>=0); ok[ok]&=(bpos[w[ok]]>=0)
                c=bpos[np.maximum(w,0)]; ok&=(c>=0)
                ok[ok]&=~taken[c[ok]]
                rows.append(np.where(ok)[0]); cols.append(c[ok])
            r=np.concatenate(rows); c=np.concatenate(cols)
            if len(r)==0: continue
            M=csr_matrix((np.ones(len(r),np.int8),(r,c)),shape=(len(free),len(B)))
            mm=maximum_bipartite_matching(M,perm_type='column')
            got=mm>=0; match[free[got]]=mm[got]
        dest[has_out,k]=nb[has_out,k]
        m=match>=0; dest[A[m],k]=B[match[m]]
        # 相手のいない列の終わり：箱の中のものは隣の枠へ、枠のものはその場（書き戻されるので無関係）
        un=A[~m]
        for v in un:
            if v>=NI: dest[v,k]=v; continue
            for o in list(ADJ[k])+list(ADJ2[k])+[OPP[k]]:
                w=nb[v,o]
                if w>=NI: dest[v,k]=w; break
        # 作りの検査：箱の中の各頂点の出どころの数
        src=np.zeros(V,np.int64); np.add.at(src,dest[dest[:,k]>=0,k],1)
        bad_in=int((src[:NI]!=1).sum()); bad_dest=int((dest[:NI,k]<0).sum())
        ring2=int(sum(1 for i in np.where(m)[0] if not np.any(nb[A[i],ADJ[k]]==B[match[i]])))
        report.append((k,len(A),len(B),int(m.sum()),bad_in,bad_dest,ring2))
    return dest, report
