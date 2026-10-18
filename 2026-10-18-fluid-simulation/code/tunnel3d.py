# 立体の担体の上の制限付きオートマトン。状態・保存・判定はすべて整数。
import numpy as np
from carrier3d_np import E, zm, sgn, dot, O, I, P
LAB=np.concatenate([E,-E])                       # 12向き ×3×2
OPP=(np.arange(12)+6)%12
G=np.stack([dot(LAB[i][None],LAB) for i in range(12)])   # 12×12×2  内積 ∈ {2+φ, φ, -φ, -(2+φ)}
def isval(z,a,b): return (z[...,0]==a)&(z[...,1]==b)
ADJ=[np.where(isval(G[k],0,1))[0] for k in range(12)]     # 内積 +φ の5本
ADJ2=[np.where(isval(G[k],0,-1))[0] for k in range(12)]   # 内積 −φ の5本
assert all(len(a)==5 for a in ADJ+ADJ2)
e0=E[0]; U=np.array([O,O,I]); W=np.array([P,[-1,0],O])
DP=dot(LAB,e0[None]); DU=dot(LAB,U[None]); DW=dot(LAB,W[None])  # 12×2

class T3:
    def __init__(s,fname,base=4,wind=4,windlab=(0,),cycle=True,obst=None):
        Z=np.load(fname); s.nb=Z['nb']; s.inner=Z['inner']; s.band=Z['band']; s.V=len(s.nb)
        s.NI=int(s.inner.sum()); V=s.V
        s.obst=np.zeros(V,bool) if obst is None else obst
        Fv=np.full(12,base,np.int64); Fv[list(windlab)]+=wind; s.F=Fv
        s.n=np.zeros((V,12),np.int64); s.n[~s.obst]=Fv; s.step=0; s.cycle=cycle
        vv=np.arange(V)[:,None]; s.tables=[]
        for c in range(5):
            dest=np.full((V,12),-1); 
            for k in range(12):
                cc=c if cycle else 0
                order=[k]+list(np.roll(ADJ[k],-cc))+list(np.roll(ADJ2[k],-cc))+[OPP[k]]
                for o in order:
                    cand=s.nb[:,o]; m=(dest[:,k]<0)&(cand>=0); dest[m,k]=cand[m]
            # 頂点ごとに巡りの位相をずらす: 位相 (c+v)%5 の表を作るため、後で頂点ごとに選ぶ
            s.tables.append(dest)
        s.prep=[]
        for c in range(5):
            ph=(np.arange(V)+c)%5
            dest=np.choose(ph[:,None].repeat(12,1),[t for t in s.tables]) if cycle else s.tables[0]
            bounce=s.obst[np.maximum(dest,0)]|(dest<0)
            kk=np.arange(12)[None,:].repeat(V,0)
            tv=np.where(bounce,vv,dest); tk=np.where(bounce,OPP[kk],kk)
            tgt=(tv*12+tk).reshape(-1)
            order=np.argsort(tgt,kind='stable'); st=tgt[order]; uq,starts=np.unique(st,return_index=True)
            src_in=np.repeat(np.arange(V)<s.NI,12); tin=(tgt//12)<s.NI
            s.prep.append(dict(order=order,uq=uq,starts=starts,out=src_in&~tin,inn=~src_in&tin,
                               hit=(bounce&s.obst[np.maximum(dest,0)]).reshape(-1)))
    def mom(s,n):
        def m(Dx): return np.array([(n*Dx[:,0]).sum(),(n*Dx[:,1]).sum()])
        return m(DP),m(DU),m(DW)
    def collide(s):
        n=s.n[:s.NI]; c=s.step%5
        Pp=np.minimum(n[:,0:6],n[:,6:12])
        n[:,0:6]-=Pp; n[:,6:12]-=Pp
        for a in range(6):
            others=[b for b in range(6) if b!=a]; b=others[(c+a)%5]
            n[:,b]+=Pp[:,a]; n[:,b+6]+=Pp[:,a]
    def run_step(s):
        n0=s.n[:s.NI]; M0=n0.sum(); P0=s.mom(n0)
        s.collide()
        n1=s.n[:s.NI]; ok_c=(n1.sum()==M0) and all(np.array_equal(a,b) for a,b in zip(P0,s.mom(n1)))
        pr=s.prep[s.step%5]; flat=s.n.reshape(-1)
        out=flat[pr['out']].sum(); inn=flat[pr['inn']].sum()
        sums=np.add.reduceat(flat[pr['order']],pr['starts'])
        new=np.zeros(s.V*12,np.int64); new[pr['uq']]=sums; s.n=new.reshape(s.V,12)
        ok=ok_c and (s.n[:s.NI].sum()==M0-out+inn)
        s.n[s.NI:]=s.F; s.n[s.obst]=0; s.step+=1
        return ok
