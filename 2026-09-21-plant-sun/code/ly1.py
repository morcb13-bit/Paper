# 検定LY1 三つの層（五角形・円環・五芒星）の網の目の大きさを整数で数える
# 事前の基準：各層の中心どうしの最も近い距離（norm2 を Z[φ] の整数組で）を出す。
#   五角形 (1,1)=φ² に対し、円環・五芒星が φ² の冪倍（norm2 で φ⁴ 倍＝(2,3)倍…）なら「相似な層」
#   そうでなければ、層は相似ではなく、別の目盛りで重なっている
import sys,pickle; sys.path.insert(0,'.')
import b13_chain_units as U
from collections import Counter,defaultdict
src=open('gen10.py').read(); exec(src[:src.index('pid={}')])
d=pickle.load(open('/home/claude/prop/g10.pkl','rb')); Q=d['Q']; S=set(Q)
RC=list({place(c,k) for (w,r,c,k) in rings})
stars=set()
for q in Q:
    for m in range(10):
        g=U.zsub(q,U.zmul(U.PHI,U.zt(m)))
        for kk in (0,1):
            if all(U.zadd(g,U.zmul(U.PHI,U.zt(kk+2*i))) in S for i in range(5)): stars.add(g)
ST=list(stars)
def val(n): return n[0]+n[1]*(1+5**.5)/2
def nearest(P,name,cell=3):
    by=defaultdict(list)
    for p in P:
        x,y=U.xy(p); by[(int(x//cell),int(y//cell))].append(p)
    c=Counter()
    for p in P:
        x,y=U.xy(p); best=None
        for dx in (-1,0,1):
            for dy in (-1,0,1):
                for r in by.get((int(x//cell)+dx,int(y//cell)+dy),()):
                    if r==p: continue
                    n=U.norm2(U.zsub(p,r))
                    if best is None or val(n)<val(best): best=n
        if best: c[best]+=1
    print(f'{name} {len(P)}個  最も近い距離の norm2（整数組: 件数）',sorted(c.items(),key=lambda t:val(t[0]))[:4])
nearest(Q,'五角形')
nearest(RC,'円環の中心',cell=6)
nearest(ST,'五芒星の中心',cell=8)
print('参考 φ² の冪の norm2: φ²=(1,1) φ⁴=(2,3) φ⁶=(5,8) φ⁸=(13,21)')
print('参考 円環の接続 2つ飛ばし',U.N_SKIP,' 連続',U.N_CONT)
