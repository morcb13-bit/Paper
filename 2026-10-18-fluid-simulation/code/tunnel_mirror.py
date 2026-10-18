# 鏡の壁で閉じた箱の、担体だけのオートマトン（枠の書き戻しなし）
import numpy as np
from tunnel3d import LAB, ADJ, OPP
from carrier3d_np import sgn
def labidx(vec):
    for j in range(12):
        if np.array_equal(LAB[j],vec): return j
    return -1
MIR=np.array([[labidx(LAB[k]*np.array([[-1],[1],[1]])) if a==0 else
               labidx(LAB[k]*np.array([[1],[-1],[1]])) if a==1 else
               labidx(LAB[k]*np.array([[1],[1],[-1]])) for k in range(12)] for a in range(3)])
assert (MIR>=0).all()
# 負の対照用：五回軸 e0 まわりの 72° 回転（ラベルの置換。同定だけに浮動小数を使う）
phi=(1+5**0.5)/2; Lf=LAB[...,0]+LAB[...,1]*phi
ax=Lf[0]/np.linalg.norm(Lf[0]); th=2*np.pi/5
K=np.array([[0,-ax[2],ax[1]],[ax[2],0,-ax[0]],[-ax[1],ax[0],0]]); Rm=np.eye(3)+np.sin(th)*K+(1-np.cos(th))*K@K
ROT=np.array([int(np.argmin(np.linalg.norm(Lf-(Rm@Lf[k]),axis=1))) for k in range(12)])
assert sorted(ROT.tolist())==list(range(12))
COMP=[np.array([sgn(LAB[k,a][None])[0] for k in range(12)]) for a in range(3)]
class MB:
    def __init__(s,fname,base=4,twisted=False):
        Z=np.load(fname); s.nb=Z['nb']; s.near=Z['near']; V=s.V=len(s.nb)
        s.n=np.full((V,12),base,np.int64); s.step=0; s.fallback=0
        vv=np.arange(V)[:,None]; s.prep=[]
        for c in range(5):
            tv=np.full((V,12),-1); tk=np.full((V,12),-1); wall=np.zeros((V,12),bool)
            for k in range(12):
                # 壁：外向き成分を持つなら鏡で反転して留まる
                newk=np.full(V,k)
                hit=np.zeros(V,bool)
                for a in range(3):
                    out=(s.near[:,2*a]&(COMP[a][k]>0))|(s.near[:,2*a+1]&(COMP[a][k]<0))
                    hit|=out
                    newk=np.where(out,MIR[a][newk],newk)
                if twisted: newk=np.where(hit,ROT[newk],newk)
                tv[hit,k]=np.arange(V)[hit]; tk[hit,k]=newk[hit]; wall[hit,k]=True
                order=[k]+list(np.roll(ADJ[k],-((c)%5)))
                for o in order:
                    cand=s.nb[:,o]; m=(tv[:,k]<0)&(cand>=0)
                    tv[m,k]=cand[m]; tk[m,k]=k
                m=tv[:,k]<0          # 前向きの辺が無い（端の角など）→ 反転して留まる
                s.fallback+=int(m.sum()); tv[m,k]=np.arange(V)[m]; tk[m,k]=OPP[k]; wall[m,k]=True
            tgt=(tv*12+tk).reshape(-1); order=np.argsort(tgt,kind='stable')
            uq,starts=np.unique(tgt[order],return_index=True)
            s.prep.append(dict(order=order,uq=uq,starts=starts,wall=wall.reshape(-1),dk=(tk.reshape(-1))))
    def P(s,n=None):
        n=s.n if n is None else n
        return np.array([[ (n*LAB[:,a,0]).sum(), (n*LAB[:,a,1]).sum()] for a in range(3)])
    def collide(s):
        n=s.n; c=s.step%5
        Pp=np.minimum(n[:,0:6],n[:,6:12]); n[:,0:6]-=Pp; n[:,6:12]-=Pp
        for a in range(6):
            others=[b for b in range(6) if b!=a]; b=others[(c+a)%5]
            n[:,b]+=Pp[:,a]; n[:,b+6]+=Pp[:,a]
    def run_step(s):
        M0=s.n.sum(); P0=s.P(); s.collide(); okc=(s.n.sum()==M0) and np.array_equal(P0,s.P())
        pr=s.prep[s.step%5]; flat=s.n.reshape(-1)
        w=pr['wall']; lab=np.tile(np.arange(12),s.V)
        src=flat[w]; dk=pr['dk'][w]; k0=lab[w]
        dP=np.array([[ (src*(LAB[dk,a,0]-LAB[k0,a,0])).sum(),(src*(LAB[dk,a,1]-LAB[k0,a,1])).sum()] for a in range(3)])
        P1=s.P()
        sums=np.add.reduceat(flat[pr['order']],pr['starts']); new=np.zeros(s.V*12,np.int64); new[pr['uq']]=sums; s.n=new.reshape(s.V,12)
        ok=okc and s.n.sum()==M0 and np.array_equal(s.P(),P1+dP)
        s.step+=1
        return ok, s.P()
