# 検定NR2 同じ中心の相似なひし形を、頂点の順序に依らずに判定する
#  細ひし形 v1..v4、M=v1+v3（=2×中心）。単数 u=ζ^m φ^n（n=1..9）について 2c = M + u·(2v−M) が
#  4頂点とも「2で割り切れて、割った値が円環の中心」になれば、中心が同じで相似（比 u）の大きなひし形がある
import pickle
exec(open('nr1.py').read().split("pairsum={}")[0].replace("print(","(lambda *a,**k:None)("))
RCs=set(RC)
def half(z): 
    return tuple(c//2 for c in z) if all(c%2==0 for c in z) else None
d=pickle.load(open('/home/claude/prop/slit14.pkl','rb')); scr=set(d['screen'])
vert_of={}
for q,a in cells.items():
    for i in range(5): vert_of.setdefault(U.zadd(q,U.zt(a+2*i)),set()).add(q)
from collections import Counter
res={}; cnt=Counter()
for c in rh:
    M=U.zadd(c[0],c[2]); hits=[]
    for n in range(1,10):
        p=U.ONE
        for _ in range(n): p=U.zmul(p,U.PHI)
        for m in range(10):
            u=U.zmul(p,U.zt(m)); imgs=[half(U.zadd(M,U.zmul(u,U.zsub(U.zadd(v,v),M)))) for v in c]
            if all(x is not None and x in RCs for x in imgs): hits.append((m,n))
    res[id(c)]=hits
    for h in {n for m,n in hits}: cnt[h]+=1
print('比 φ^n ごとに、同じ中心の相似な大きなひし形（円環の中心）をもつ細ひし形の数',dict(sorted(cnt.items())),'/ 28')
print('回りの m（比 φ^3 のとき）',Counter(m for hs in res.values() for m,n in hs if n==3))
band=[c for c in rh if set().union(*[vert_of.get(v,set()) for v in c])&scr]
for c in band: print('  帯のひし形 x=%6.2f'%(sum(U.xy(v)[0] for v in c)/4),res[id(c)])
none=[c for c in rh if not res[id(c)]]
print('相似なひし形をもたない細ひし形',len(none),' その y:',sorted(round(sum(U.xy(v)[1] for v in c)/4,1) for c in none))
