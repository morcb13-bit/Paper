# 検定DR：正5胞体を4次元で回し、固定した3次元への影で H どうしの角を測る
# 事前の基準：
#  記録の幅  CH4 109.4712（対称）/ NH3 106.759–106.779 / H2O 104.4776±0.0019 / H3O+ 111.90–112.00
#  R1  非H の頂点が張る部分空間の中で落とす向きを一様に回す → 平均の角
#  R2  面P（u,m）と面Q を同じ速さで回す二重回転（等傾）、落とす向きは実験室の第4軸に固定 → 時間平均
#  負の対照 N  落とす向きを3次元球面 S3 全体で一様に選ぶ（配置によらず同じ値になるはず → 弁別しないことの確認）
#  OK = 記録の幅に入る。NH3 と H3O+ は席の配置が同じなので、配置だけで決まる規則は片方が必ず NG になる（先に明記）
# 浮動小数は回転と角の計算だけ
import math, random
import numpy as np
d=5; V=np.array([[(d if k==i else 0)-1 for k in range(d)] for i in range(d)],float)
# 4次元の正規直交基底へ（和0の超平面）
B=np.linalg.qr(np.vstack([np.ones(5),np.eye(5)[:4]]).T)[0][:,1:5]; X=V@B   # 5頂点の4次元座標
def ang_after_drop(H,w):
    w=w/np.linalg.norm(w); P=[h-np.dot(h,w)*w for h in H]
    cs=[np.dot(P[i],P[j])/np.linalg.norm(P[i])/np.linalg.norm(P[j]) for i in range(len(P)) for j in range(i+1,len(P))]
    return np.degrees(np.arccos(np.mean(cs)))
CONF={"CH4":([0,1,2,3],[4]),"NH3/H3O+":([0,1,2],[3,4]),"H2O":([0,1],[2,3,4])}
rng=np.random.default_rng(1)
print("R1：非Hの部分空間で一様に落とす向き")
for nm,(h,o) in CONF.items():
    H=X[h]; O=X[o]; Q,_=np.linalg.qr(O.T); k=Q.shape[1]
    a=[ang_after_drop(H,Q@rng.standard_normal(k)) for _ in range(20000)]
    print(f"  {nm:9} 平均 {np.mean(a):8.3f}  幅 {np.min(a):.2f}–{np.max(a):.2f}")
print("R2：等傾の二重回転（面P・面Qを同じ速さ）、落とす向きは第4軸固定")
for nm,(h,o) in CONF.items():
    H0=X[h]; m=-H0.sum(0)
    others=np.linalg.qr(np.vstack([H0, m]).T)[0]
    # 面P：m と、H すべてに直交する向き u（あれば）。面Q：その補空間
    U=np.linalg.svd(H0)[2]; u=U[-1] if np.linalg.matrix_rank(H0)<4 else None
    mh=m/np.linalg.norm(m)
    if u is None: u=np.linalg.svd(np.vstack([H0[:3],m]))[2][-1]
    u=u-np.dot(u,mh)*mh; u/=np.linalg.norm(u)
    Qb=np.linalg.svd(np.vstack([mh,u]))[2][2:]
    lab=u.copy()            # 初めは本当の角が見える向き
    vals=[]
    for t in np.linspace(0,2*math.pi,4000,endpoint=False):
        c,s=math.cos(t),math.sin(t)
        R=np.outer(mh,c*mh+s*u)+np.outer(u,-s*mh+c*u)+np.outer(Qb[0],c*Qb[0]+s*Qb[1])+np.outer(Qb[1],-s*Qb[0]+c*Qb[1])
        vals.append(ang_after_drop(H0@R.T,lab))
    print(f"  {nm:9} 平均 {np.mean(vals):8.3f}  幅 {np.min(vals):.2f}–{np.max(vals):.2f}")
print("負の対照N：S3 全体で一様")
for nm,(h,o) in CONF.items():
    a=[ang_after_drop(X[h],rng.standard_normal(4)) for _ in range(20000)]
    print(f"  {nm:9} 平均 {np.mean(a):8.3f}")
