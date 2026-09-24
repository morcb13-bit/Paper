# NR3 帯のひし形の相似な相手は担体の外にあるか：行を16まで積んだ円環の中心で探し直す
import os,sys,pickle
os.environ['B13_ROWS']='17'
sys.path.insert(0,'/home/claude/Paper/2026-09-21-plant-sun/code')
import b13_chain_units as U
rows,place,offs=U.build_stack(); RC=set(sum(place,[]))
d=pickle.load(open('/home/claude/prop/slit14.pkl','rb'))
os.environ['B13_ROWS']='14'
import importlib; importlib.reload(U)
r14,p14,_=U.build_stack(); print('14行の円環の中心が17行の中に含まれる',set(sum(p14,[]))<=RC,' 17行の円環',len(RC))
cells=d['cells']; scr=set(d['screen'])
F=U.gaps(cells); rh=[cyc for a,cyc in F if min(U.GAP_NAME,key=lambda x:abs(x-a))==0.8123]
vert_of={}
for q,a in cells.items():
    for i in range(5): vert_of.setdefault(U.zadd(q,U.zt(a+2*i)),set()).add(q)
half=lambda z: tuple(c//2 for c in z) if all(c%2==0 for c in z) else None
band=[c for c in rh if set().union(*[vert_of.get(v,set()) for v in c])&scr]
for c in band:
    M=U.zadd(c[0],c[2]); hits=[]
    for n in range(1,10):
        p=U.ONE
        for _ in range(n): p=U.zmul(p,U.PHI)
        for m in range(10):
            u=U.zmul(p,U.zt(m)); imgs=[half(U.zadd(M,U.zmul(u,U.zsub(U.zadd(v,v),M)))) for v in c]
            if all(x is not None and x in RC for x in imgs): hits.append((m,n))
    ext=[]
    if hits:
        m,n=hits[0]; p=U.ONE
        for _ in range(n): p=U.zmul(p,U.PHI)
        u=U.zmul(p,U.zt(m)); imgs=[half(U.zadd(M,U.zmul(u,U.zsub(U.zadd(v,v),M)))) for v in c]
        ext=[r for r in range(len(place)) if any(x in place[r] for x in imgs)]
    print('帯のひし形 x=%6.2f  相似な大きなひし形 %s  頂点の円環が乗る行 %s'%(sum(U.xy(v)[0] for v in c)/4,hits,[r+1 for r in ext]))
