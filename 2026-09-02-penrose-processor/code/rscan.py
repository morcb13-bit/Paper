import pickle,math,sys,time
from collections import defaultdict
import b13_chain_units_big as U
B=pickle.load(open("big_base.pkl","rb")); V10,XY10,CEN=B["V10"],B["XY10"],B["CEN"]
Rcar=max(math.hypot(*v) for v in XY10.values())/10.0
G=defaultdict(list)
for p in V10:
    x,y=XY10[p]; G[(int(x//50),int(y//50))].append(p)
CAP=45.0; z2=U.zt(2)
def closes(c):
    cx,cy=U.xy(c); lim=min(CAP,Rcar-math.hypot(cx,cy)/10.0)
    rr=int(lim*10)//50+2; gx,gy=int(cx//50),int(cy//50)
    g=defaultdict(list)
    for i in range(gx-rr,gx+rr+1):
        for j in range(gy-rr,gy+rr+1):
            for p in G.get((i,j),()):
                d=math.hypot(XY10[p][0]-cx,XY10[p][1]-cy)/10.0
                if d<=lim: g[round(d,9)].append(p)
    R=0.0
    for d in sorted(g):
        if all(U.zadd(c,U.zmul(U.zsub(p,c),z2)) in V10 for p in g[d]): R=d
        else: return R,lim
    return R,lim
a,b=int(sys.argv[1]),int(sys.argv[2])
t0=time.time(); out=[]
for c in CEN[a:b]:
    if Rcar-math.hypot(*U.xy(c))/10.0<=CAP: continue
    out.append((c,)+closes(c))
pickle.dump(out,open("rscan_%d.pkl"%a,"wb"))
print("%d..%d  %d 個  %.1fs"%(a,b,len(out),time.time()-t0))
