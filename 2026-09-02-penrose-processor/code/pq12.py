import json, math, sys, time
from collections import defaultdict
import b13_chain_units as U
CAP=float(sys.argv[1]); TMAX=int(sys.argv[2])
src=json.load(open("carrier_big.json"))["cells"]
cells={}
for k,a in src.items():
    z=tuple(int(x) for x in k.split(","))
    if math.hypot(*U.xy(z))<=CAP: cells[z]=a
adj=defaultdict(set)
for q,a in cells.items():
    vs=[U.zadd(q,U.zt(a+2*i)) for i in range(5)]
    for i in range(5):
        u,w=vs[i],vs[(i+1)%5]; adj[u].add(w); adj[w].add(u)
XY={p:U.xy(p) for p in adj}
V=list(adj); VI={v:i for i,v in enumerate(V)}
E=[]; EID={}
for v in V:
    for u in adj[v]: EID[(u,v)]=len(E); E.append((u,v))
n=len(E)
INC=[[EID[(u,v)] for u in adj[v]] for v in V]
OUT=[[EID[(v,u)] for u in adj[v]] for v in V]
DEG=[len(adj[v]) for v in V]
HEADV=[VI[e[1]] for e in E]
COFF={2:12,3:8,4:6}; CP={2:6,3:4,4:3}
R2=[U.norm2(v) for v in V]; RAD=[math.hypot(*XY[v]) for v in V]
print(f"半径{CAP:.0f} 五角形{len(cells)} 頂点{len(V)} 有向辺{n}")
def za(x,y): return (x[0]+y[0],x[1]+y[1])
def zm(x,y):
    a1,b1=x; a2,b2=y
    return (a1*a2+b1*b2,a1*b2+a2*b1+b1*b2)
from fractions import Fraction
_F=[1,1]
while len(_F)<200: _F.append(_F[-1]+_F[-2])
_PHI=Fraction(_F[150],_F[149])
def zq(x): return x[0]+x[1]*_PHI
def zf(x):
    v=zq(x)
    return float(v) if abs(v.numerator)<10**300 and abs(v.denominator)<10**300 else float('nan')
def ratio(L,R):
    return float(Fraction(zq(L),zq(R)))
start=min(range(len(V)),key=lambda i:RAD[i])
MEAS=sorted({5,10,20,30,40,50,60,80,100,120,140,160,200,240,280,TMAX})
def run(quantum):
    a=[0]*n
    for e in INC[start]: a[e]=1
    res={}; N0=sum(x*x for x in a); bad=0
    nv=len(V)
    for t in range(1,TMAX+1):
        b=[0]*n
        for vi in range(nv):
            inc=INC[vi]; out=OUT[vi]; d=DEG[vi]
            S=0
            for e in inc: S+=a[e]
            if quantum:
                c=COFF[d]*S
                for k in range(d): b[out[k]]=(c-12*a[inc[k]])>>1     # 毎歩 2 で割る
            else:
                c=CP[d]*S
                for k in range(d): b[out[k]]=c
        a=b
        if quantum:
            N1=sum(x*x for x in a)
            if N1!=36*N0: bad+=1
            N0=N1
        if t in MEAS:
            S=(0,0);Q=(0,0);rmax=0.0
            vw=defaultdict(int)
            for i in range(n):
                x=a[i]
                if x: vw[HEADV[i]]+= x*x if quantum else x
            for v,w in vw.items():
                S=za(S,(w,0)); Q=za(Q,zm((w,0),R2[v]))
                if RAD[v]>rmax: rmax=RAD[v]
            res[t]=(S,Q,rmax)
    return res,bad
t0=time.time(); A,bad=run(True); print("装置A %.0fs  二乗和が36倍からずれた歩: %d"%(time.time()-t0,bad))
t0=time.time(); B,_=run(False); print("古典B %.0fs"%(time.time()-t0))
print("\n  t | 装置A <r^2> | 古典B <r^2> | 前線半径")
for t in MEAS:
    SA,QA,ra=A[t]; SB,QB,_=B[t]
    print(f"{t:4d} | {float(Fraction(zq(QA),zq(SA))):11.2f} | {float(Fraction(zq(QB),zq(SB))):11.2f} | {ra:7.1f}")
print("\n倍率 t->2t")
for t in MEAS:
    if 2*t in A:
        for lab,M in (("A",A),("B",B)):
            S1,Q1,_=M[t]; S2,Q2,_=M[2*t]
            print(f"  {t:3d}->{2*t:3d} {lab}: {ratio(zm(Q2,S1),zm(Q1,S2)):.3f}")
