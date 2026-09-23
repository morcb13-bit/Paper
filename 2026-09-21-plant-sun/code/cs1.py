# 検定CS1 中心の五芒星＝レジスタか。PS1 の「指しを五芒星のまわりで回す」と AD3 の「扇を渡る」が同じものかを整数で見る
# 事前の基準
#   OK  中心の五芒星を囲む5枚が主の扇 0,2,4,6,8 に一枚ずつ属し、中心からの向き ζ^m の m と扇 w が一つの規則（m = w + 定数 mod 10）で対応する
#       さらに、その5枚の間（ζ^(m±1) の向き）にずらしの扇がある：中心から φ·ζ^(m+1) 方向へ伸ばした先の五角形がずらしの扇 w+1 に属する
#   NG  対応が一つの規則にならない
import sys,pickle; sys.path.insert(0,'.')
import b13_chain_units as U
from collections import defaultdict
d=pickle.load(open('/home/claude/prop/g10.pkl','rb')); Q=d['Q']; cw=d['cell_w']; S=set(Q)
g=None
for q in Q:
    for m in range(10):
        c=U.zsub(q,U.zmul(U.PHI,U.zt(m)))
        for k in (0,1):
            ring=[U.zadd(c,U.zmul(U.PHI,U.zt(k+2*i))) for i in range(5)]
            if all(r in S for r in ring) and len({w for r in ring for w in cw[r]})==5: g=c;R=ring;K=k
print('中心の五芒星',g,'= 担体の中心 z0',g==d['z0'])
pairs=[]
for i,r in enumerate(R):
    m=(K+2*i)%10; pairs.append((m,sorted(cw[r])))
print('向き m と扇',pairs)
off={(w[0]-m)%10 for m,w in pairs}
print('m と扇の差が一つ',off)
# 間の向き：中心から φ²·ζ^(m+1) の倍数方向に最初に当たる五角形の扇
res=[]
for m,w in pairs:
    for s in (1,-1):
        hit=None
        for n in range(1,40):
            for base in (U.PHI,U.zmul(U.PHI,U.PHI)):
                p=U.zadd(g,U.zmul(U.zmul(base,U.zt(m+s)),(n,0,0,0)))
                if p in S: hit=p;break
            if hit: break
        res.append((m,s,sorted(cw[hit]) if hit else None))
print('間の向き（m±1）で最初に当たる五角形の扇',res)
