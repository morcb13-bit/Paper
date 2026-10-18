import numpy as np
from judge_d0 import za, zs, sgn
Z=np.load('d0e_map.npz'); acc=Z['acc']; mass=Z['mass']; deg=Z['deg']
phi=(1+5**0.5)/2
def S(m): return (int(acc[m,0].sum()),int(acc[m,1].sum())), int(mass[m].sum())
St,Mt=S(np.ones(len(deg),bool))
print('全体 ずれ≈%.4f  粒子/頂点≈%.1f'%((St[0]+St[1]*phi)/Mt, Mt/len(deg)/1000))
bad=[]
for d in range(4,13):
    m=deg==d
    if m.sum()==0: continue
    s,M=S(m)
    # s/M と St/Mt の比が 0.6〜1.4 か（符号込み）: 10·s·Mt − 6·St·M, 10·s·Mt − 14·St·M の符号
    a=zs(s,10*Mt); b=zs(St,M)
    lo=sgn(za(a,zs(b,-6))); hi=sgn(za(a,zs(b,-14)))
    ok=(lo<0 and hi>0) if sgn(St)<0 else (lo>0 and hi<0)
    if not ok: bad.append(d)
    print('  向き%2d本 粒子/頂点≈%6.1f  ずれ≈%.4f  %s'%(d,M/m.sum()/1000,(s[0]+s[1]*phi)/M,'OK' if ok else '外れ'))
print('判定:', '全体が一様にずれている（10本・12本は粒子が多いだけ）' if not bad else '並びで向きが偏る: 向き%s本'%bad)
