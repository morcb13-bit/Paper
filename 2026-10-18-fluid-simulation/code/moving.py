# 静止した空気の中を物体が進む。担体だけのオートマトン（札の入れ替えだけ）。
import numpy as np, sys, json, time
from tunnel3d import LAB, ADJ, ADJ2, OPP
from tunnel_mirror import MIR, COMP
from carrier3d_np import sgn, zm
body=sys.argv[1]; k_step=int(sys.argv[2])        # k_step=0 は止めたまま
Z=np.load('longbox.npz'); nb=Z['nb']; near=Z['near']; x=Z['xpar']; V=len(nb); vv=np.arange(V)
def sq(z): return zm(z,z)
def body_mask(cx):                                # 中心 (cx,0,0)。判定は Z[φ] の整数比較
    dx=x[:,0]-np.array([cx,0]); Y=x[:,1]; Zc=x[:,2]
    if body=='ball': return sgn(sq(dx)+sq(Y)+sq(Zc)-np.array([64,0]))<0
    s=1 if body=='wing+' else -1
    # 先に箱で絞る（値が小さいうちに整数比較）。大きな整数は Python の整数で比べる（桁あふれを避ける）
    near_=(sgn(dx-np.array([20,0]))<0)&(sgn(dx+np.array([20,0]))>0)&(sgn(Y-np.array([20,0]))<0)&(sgn(Y+np.array([20,0]))>0)&(sgn(Zc-np.array([20,0]))<0)&(sgn(Zc+np.array([20,0]))>0)
    out=np.zeros(len(dx),bool); idx=np.where(near_)[0]
    dx=dx[idx].astype(object); Y=Y[idx].astype(object); Zc=Zc[idx].astype(object)
    u=4*dx+s*3*Y; w=-s*3*dx+4*Y
    A,B,C=144,25,324                               # 12², (2.5·2)²(=25 を 4 倍した厚み²で扱う), 18²
    # (u/5)²/12² + (w/5)²/2.5² + z²/18² < 1  → 2.5²=25/4 なので w の項は 4 倍
    val=sq(u)*(4*C)*1+sq(w)*(4*A*C)+sq(Zc)*(25*4*A*B//4)
    val=sq(u)*(B*C)+sq(w)*(4*A*C)+sq(Zc)*(25*A*B)
    val=val-np.array([25*A*B*C,0],dtype=object)
    p,q=val[:,0],val[:,1]; Aa=2*p+q
    neg=np.array([ (a<=0 and b<=0 and not(a==0 and b==0)) or (a>0 and b<0 and a*a<5*b*b) or (a<0 and b>0 and a*a>5*b*b) for a,b in zip(Aa,q)],bool)
    out[idx]=neg
    return out
# 壁と担体だけの渡り表（物体は毎手番あとから当てる）
tabs=[]
for c in range(5):
    tv=np.full((V,12),-1); tk=np.full((V,12),-1); wall=np.zeros((V,12),bool)
    for k in range(12):
        newk=np.full(V,k); hit=np.zeros(V,bool)
        for a in range(3):
            out=(near[:,2*a]&(COMP[a][k]>0))|(near[:,2*a+1]&(COMP[a][k]<0)); hit|=out; newk=np.where(out,MIR[a][newk],newk)
        tv[hit,k]=vv[hit]; tk[hit,k]=newk[hit]; wall[hit,k]=True
        for o in [k]+list(np.roll(ADJ[k],-c)):
            cand=nb[:,o]; m=(tv[:,k]<0)&(cand>=0); tv[m,k]=cand[m]; tk[m,k]=k
        m=tv[:,k]<0; tv[m,k]=vv[m]; tk[m,k]=OPP[k]; wall[m,k]=True
    tabs.append((tv,tk,wall))
labs=np.tile(np.arange(12),(V,1))
def P(n): return np.array([[(n*LAB[:,a,0]).sum(),(n*LAB[:,a,1]).sum()] for a in range(3)])
def dP(cnt,k0,k1): return np.array([[(cnt*(LAB[k1,a,0]-LAB[k0,a,0])).sum(),(cnt*(LAB[k1,a,1]-LAB[k0,a,1])).sum()] for a in range(3)])
XPOS=np.array([sgn(LAB[k,0][None])[0]>0 for k in range(12)]); XNEG=np.array([sgn(LAB[k,0][None])[0]<0 for k in range(12)])
cx=-90; inb=body_mask(cx); n=np.full((V,12),4,np.int64); n[inb]=0
WARM=200; END=80; ok_all=True; rec=[]; t=0; t0=time.time()
import os
FR=os.environ.get('FRAMES')=='1'
if FR:
    phi_=(1+5**0.5)/2; xf_=x[...,0]+x[...,1]*phi_
    slab_=np.abs(xf_[:,2])<5; NX,NY=120,30
    bx_=np.clip(((xf_[:,0]+120)/2).astype(int),0,NX-1); by_=np.clip(((xf_[:,1]+30)/2).astype(int),0,NY-1)
    cnt_=np.zeros((NX,NY)); np.add.at(cnt_,(bx_[slab_],by_[slab_]),1)
    LXf=LAB[:,0,0]+LAB[:,0,1]*phi_; LYf=LAB[:,1,0]+LAB[:,1,1]*phi_
    frames=[]
def displace(newly):
    """新しく物体に入った頂点の粒子を、前へ押して外の隣へ移す（札の入れ替えだけ）"""
    global n
    d=np.zeros((3,2),np.int64)
    for v in np.where(newly)[0]:
        for k in range(12):
            c=n[v,k]
            if c==0: continue
            k2=MIR[0][k] if XNEG[k] else k
            dest=-1; cand=[o for o in [k2]+list(ADJ[k2])+list(ADJ2[k2])+[OPP[k2]] if nb[v,o]>=0]
            for o in cand:
                if not inb[nb[v,o]]: dest=nb[v,o]; break
            if dest<0: dest=nb[v,cand[0]]          # 外の隣が無ければ前の隣へ（次の段でさらに押す）
            n[dest,k2]+=c; n[v,k]=0
            d+=dP(np.array([c]),np.array([k]),np.array([k2]))
    return d
def push_all():
    d=np.zeros((3,2),np.int64)
    for _ in range(20):
        left=inb&(n.sum(1)>0)
        if not left.any(): break
        d+=displace(left)
    return d
while True:
    moving=(k_step>0 and t>=WARM)
    if moving and (t-WARM)%k_step==0 and t>WARM:
        cx+=1
        if cx>END: break
        new=body_mask(cx); newly=new&~inb; inb=new
        M0=n.sum(); P0=P(n); d=displace(newly); d+=push_all()
        ok_all&=bool(n.sum()==M0 and np.array_equal(P(n),P0+d))
        n[inb]=n[inb]  # 物体の中に残った粒子は無い（下で確認）
        dpush=d
    else: dpush=np.zeros((3,2),np.int64)
    if k_step==0 and t>=WARM+400: break
    M0=n.sum(); P0=P(n)
    cc=t%5; Pp=np.minimum(n[:,0:6],n[:,6:12]); n[:,0:6]-=Pp; n[:,6:12]-=Pp
    for a in range(6):
        others=[b for b in range(6) if b!=a]; b=others[(cc+a)%5]; n[:,b]+=Pp[:,a]; n[:,b+6]+=Pp[:,a]
    tv,tk,wall=tabs[t%5]
    hitb=inb[tv]&~wall                                   # 物体へ渡ろうとした → 跳ね返り
    tv2=np.where(hitb,vv[:,None],tv); tk2=np.where(hitb,OPP[labs],tk)
    dwall=dP(n[wall],labs[wall],tk[wall]); dbody=dP(n[hitb],labs[hitb],tk2[hitb])
    new=np.zeros(V*12,np.int64); np.add.at(new,(tv2*12+tk2).reshape(-1),n.reshape(-1)); n=new.reshape(V,12)
    ok_all&=bool(n.sum()==M0 and np.array_equal(P(n),P0+dwall+dbody) and n[inb].sum()==0)
    if (k_step==0 and t>=WARM) or (moving and -40<=cx<40):
        rec.append([int(x_) for x_ in (dbody+dpush).reshape(-1)])
    if FR and moving and (t-WARM)%3==0:
        m=n.sum(1).astype(float); vx=(n*LXf).sum(1); vy=(n*LYf).sum(1)
        G=[np.zeros((NX,NY)) for _ in range(3)]
        for g,val in zip(G,(m,vx,vy)): np.add.at(g,(bx_[slab_],by_[slab_]),val[slab_])
        frames.append(dict(cx=cx,m=(G[0]/np.maximum(cnt_,1)).round(2).tolist(),vx=(G[1]/np.maximum(cnt_,1)).round(2).tolist(),vy=(G[2]/np.maximum(cnt_,1)).round(2).tolist()))
    t+=1
if FR: json.dump(dict(frames=frames,cnt=(cnt_>0).astype(int).tolist()),open('mvfr_%s_%d.json'%(body,k_step),'w'))
json.dump(dict(body=body,k=k_step,ok=ok_all,rec=rec,steps=t),open('mv_%s_%d.json'%(body,k_step),'w'))
print(body,k_step,'帳簿',ok_all,'手番',t,'記録',len(rec),'sec',round(time.time()-t0))
