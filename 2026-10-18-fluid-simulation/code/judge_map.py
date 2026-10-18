import numpy as np
from judge_d0 import za, zs, sgn
Z=np.load('d0e_map.npz'); acc=Z['acc']; depth=Z['depth']; deg=Z['deg']
phi=(1+5**0.5)/2
def S(m): return (int(acc[m,0].sum()),int(acc[m,1].sum()))
def zabs(x): return x if sgn(x)>=0 else zs(x,-1)
def gt3(Sa,Na,Sb,Nb):   # |Sa|/Na > 3|Sb|/Nb
    Na,Nb=int(Na),int(Nb)
    return sgn(za(zs(zabs(Sa),Nb),zs(zabs(Sb),-3*Na)))>0
allm=np.ones(len(depth),bool); St=S(allm); Nt=len(depth)
f=lambda x:(x[0]+x[1]*phi)
print('全体 N=%d ΣT_v≈%.1f百万 頂点あたり≈%.0f'%(Nt,f(St)/1e6,f(St)/Nt))
groups=[('深さ1',depth==1),('深さ2',depth==2),('深さ3-4',(depth>=3)&(depth<=4)),('深さ5以上',depth>=5)]
for nm,m in groups: s=S(m); print('  %s N=%d ΣT_v≈%.1f百万 頂点あたり≈%.0f'%(nm,m.sum(),f(s)/1e6,f(s)/m.sum()))
m12=(depth>=1)&(depth<=2); m5=depth>=5
print('縁説（深さ1-2 が 深さ5以上 の3倍超）:',gt3(S(m12),m12.sum(),S(m5),m5.sum()))
hit=[]
for d in range(4,13):
    m=deg==d
    if m.sum()==0: continue
    s=S(m); h=gt3(s,m.sum(),St,Nt); hit.append(h)
    print('  向き%2d本 N=%6d ΣT_v≈%8.1f百万 頂点あたり≈%7.0f %s'%(d,m.sum(),f(s)/1e6,f(s)/m.sum(),'★' if h else ''))
print('並び説（どれかの組が全体の3倍超）:',any(hit))
