import numpy as np, sys, json, time
from tunnel3d import LAB, ADJ, OPP
from tunnel_mirror import MIR, COMP
from carrier3d_np import sgn, zm
fan_on=(sys.argv[1]=='fan'); body=sys.argv[2] if len(sys.argv)>2 else 'ball'
Z=np.load('mbox43.npz'); nb=Z['nb']; near=Z['near']; x=Z['xpar']; V=len(nb)
def lt(z,c): return sgn(z-np.array([c,0]))<0
def gt(z,c): return sgn(z-np.array([c,0]))>0
def sq(z): return zm(z,z)
r2=sq(x[:,1])+sq(x[:,2])
fan=gt(x[:,0],-30)&lt(x[:,0],-24)&(sgn(r2-np.array([18*18,0]))<0)
jet=gt(x[:,0],-20)&lt(x[:,0],-10)&(sgn(r2-np.array([12*12,0]))<0)
dx=x[:,0]-np.array([5,0]); s2=sq(dx)+r2
if body=='ball': ball=sgn(s2-np.array([100,0]))<0
else:
    sgnr=1 if body=='wing+' else -1
    X=dx; Y=x[:,1]; Zc=x[:,2]
    u=4*X+sgnr*3*Y; w=-sgnr*3*X+4*Y          # 5倍した回転後の座標
    # (u/5)²/14² + (w/5)²/3² + z²/20² < 1  を整数倍して比べる
    A=14*14; B=3*3; C=20*20
    val=sq(u)*(B*C)+sq(w)*(A*C)+sq(Zc)*(25*A*B)
    ball=sgn(val-np.array([25*A*B*C,0]))<0
print('fan',fan.sum(),'jet',jet.sum(),'ball',ball.sum())
n=np.full((V,12),4,np.int64); n[ball]=0
vv=np.arange(V); prep=[]
for c in range(5):
    tv=np.full((V,12),-1); tk=np.full((V,12),-1); kind=np.zeros((V,12),np.int8)  # 1壁 2球
    for k in range(12):
        newk=np.full(V,k); hit=np.zeros(V,bool)
        for a in range(3):
            out=(near[:,2*a]&(COMP[a][k]>0))|(near[:,2*a+1]&(COMP[a][k]<0)); hit|=out
            newk=np.where(out,MIR[a][newk],newk)
        tv[hit,k]=vv[hit]; tk[hit,k]=newk[hit]; kind[hit,k]=1
        for o in [k]+list(np.roll(ADJ[k],-c)):
            cand=nb[:,o]; m=(tv[:,k]<0)&(cand>=0)
            b=m&ball[np.maximum(cand,0)]
            tv[b,k]=vv[b]; tk[b,k]=OPP[k]; kind[b,k]=2
            m&=~b; tv[m,k]=cand[m]; tk[m,k]=k
        m=tv[:,k]<0; tv[m,k]=vv[m]; tk[m,k]=OPP[k]; kind[m,k]=1
    tgt=(tv*12+tk).reshape(-1); order=np.argsort(tgt,kind='stable'); uq,st=np.unique(tgt[order],return_index=True)
    prep.append((order,uq,st,kind.reshape(-1),tk.reshape(-1)))
lab=np.tile(np.arange(12),V)
def P(n): return np.array([[(n*LAB[:,a,0]).sum(),(n*LAB[:,a,1]).sum()] for a in range(3)])
def dPof(src,k0,k1): return np.array([[(src*(LAB[k1,a,0]-LAB[k0,a,0])).sum(),(src*(LAB[k1,a,1]-LAB[k0,a,1])).sum()] for a in range(3)])
xneg=np.array([sgn(LAB[k,0][None])[0]<0 for k in range(12)])
ok_all=True; Dw=[];Lyw=[];Lzw=[];Jw=[]; acc=np.zeros(4,object); frames=[]
# 表示用の粗い格子（z の薄い板、表示だけに浮動小数）
phi=(1+5**0.5)/2; xf=x[...,0]+x[...,1]*phi
slab=np.abs(xf[:,2])<4; gx=((xf[:,0]+43)//4).astype(int); gy=((xf[:,1]+43)//4).astype(int)
t0=time.time(); VACC=np.zeros((V,3,2),np.int64); SNAP=[]
for t in range(2000):
    M0=n.sum(); P0=P(n)
    # 衝突
    cc=t%5; Pp=np.minimum(n[:,0:6],n[:,6:12]); n[:,0:6]-=Pp; n[:,6:12]-=Pp
    for a in range(6):
        others=[b for b in range(6) if b!=a]; b=others[(cc+a)%5]; n[:,b]+=Pp[:,a]; n[:,b+6]+=Pp[:,a]
    dP=np.zeros((3,2),np.int64)
    # 送風機（一方通行の鏡）
    if fan_on:
        for k in np.where(xneg)[0]:
            m=fan&(n[:,k]>0); kk=MIR[0][k]
            cnt=n[m,k]; dP+=np.array([[cnt.sum()*(LAB[kk,a,0]-LAB[k,a,0]),cnt.sum()*(LAB[kk,a,1]-LAB[k,a,1])] for a in range(3)])
            n[m,kk]+=cnt; n[m,k]=0
    order,uq,st,kind,tk=prep[t%5]; flat=n.reshape(-1)
    w=kind>0; dP+=dPof(flat[w],lab[w],tk[w])
    bm=kind==2; hs=flat[bm]; hk=lab[bm]
    D=int(hs[np.isin(hk,np.where(np.array([sgn(LAB[k,0][None])[0]>0 for k in range(12)]))[0])].sum()-hs[np.isin(hk,np.where(xneg)[0])].sum())
    Ly=int(sum(hs[hk==k].sum()*int(sgn(LAB[k,1][None])[0]) for k in range(12)))
    Lz=int(sum(hs[hk==k].sum()*int(sgn(LAB[k,2][None])[0]) for k in range(12)))
    sums=np.add.reduceat(flat[order],st); new=np.zeros(V*12,np.int64); new[uq]=sums; n=new.reshape(V,12)
    ok_all&=bool(n.sum()==M0 and np.array_equal(P(n),P0+dP))
    if t>=600:
        Jp=(n[jet]*LAB[:,0,0]).sum(),(n[jet]*LAB[:,0,1]).sum()
        acc+=np.array([D,Ly,Lz,0],object); acc[3]=acc[3]+0
        if 'jacc' not in dir(): jacc=np.zeros(2,np.int64)
        jacc+=np.array(Jp,np.int64)
        if (t-600)%50==49:
            Dw.append(int(acc[0]));Lyw.append(int(acc[1]));Lzw.append(int(acc[2]));Jw.append(jacc.tolist()); acc[:]=0; jacc[:]=0
    VACC[:,:,0]+=np.stack([(n*LAB[:,a,0]).sum(1) for a in range(3)],1); VACC[:,:,1]+=np.stack([(n*LAB[:,a,1]).sum(1) for a in range(3)],1)
    if t%100==99: SNAP.append(VACC.copy()); VACC[:]=0
    if t%40==39:
        pxf=(n*(LAB[:,0,0]+LAB[:,0,1]*phi)).sum(1); pyf=(n*(LAB[:,1,0]+LAB[:,1,1]*phi)).sum(1)
        g=np.zeros((22,22)); h=np.zeros((22,22)); c_=np.zeros((22,22))
        np.add.at(g,(gx[slab],gy[slab]),pxf[slab]); np.add.at(h,(gx[slab],gy[slab]),pyf[slab]); np.add.at(c_,(gx[slab],gy[slab]),n[slab].sum(1))
        frames.append(dict(t=t+1,px=np.round(g/np.maximum(c_,1),3).tolist(),py=np.round(h/np.maximum(c_,1),3).tolist()))
cells_ball=np.zeros((22,22),int); np.add.at(cells_ball,(gx[slab&ball],gy[slab&ball]),1)
cells_fan=np.zeros((22,22),int); np.add.at(cells_fan,(gx[slab&fan],gy[slab&fan]),1)
np.savez_compressed('vfield_%s_%s.npz'%(sys.argv[1],body),snap=np.array(SNAP),mass=n.sum(1),ball=ball,fan=fan)
json.dump(dict(fan=fan_on,ok=ok_all,D=Dw,Ly=Lyw,Lz=Lzw,J=Jw,frames=frames,ball=(cells_ball>0).astype(int).tolist(),fanc=(cells_fan>0).astype(int).tolist()),open('w3_%s_%s.json'%(sys.argv[1],body),'w'))
print(sys.argv[1],'帳簿',ok_all,'sec',round(time.time()-t0))
