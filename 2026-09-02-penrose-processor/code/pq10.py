import json, math, time
from collections import defaultdict
import b13_chain_units as U

cells={tuple(int(x) for x in k.split(",")):a for k,a in json.load(open("carrier_1245.json"))["cells"].items()}
adj=defaultdict(set)
for q,a in cells.items():
    vs=[U.zadd(q,U.zt(a+2*i)) for i in range(5)]
    for i in range(5):
        u,w=vs[i],vs[(i+1)%5]; adj[u].add(w); adj[w].add(u)
XY={p:U.xy(p) for p in adj}
V=list(adj); VI={v:i for i,v in enumerate(V)}
E=[]; EID={}
for v in V:
    for u in adj[v]:
        EID[(u,v)]=len(E); E.append((u,v))     # u->v : v で散る
n=len(E)
# 各頂点 v: 入ってくる辺 (u->v) と、同じ隣 u へ出る辺 (v->u)
INC=[[EID[(u,v)] for u in adj[v]] for v in V]
OUT=[[EID[(v,u)] for u in adj[v]] for v in V]
DEG=[len(adj[v]) for v in V]
COFF={2:12,3:8,4:6}; CP={2:6,3:4,4:3}
R2=[U.norm2(v) for v in V]; RAD=[math.hypot(*XY[v]) for v in V]
print("頂点",len(V),"有向辺",n)

def za(x,y): return (x[0]+y[0],x[1]+y[1])
def zm(x,y):
    a1,b1=x; a2,b2=y
    return (a1*a2+b1*b2, a1*b2+a2*b1+b1*b2)
def zf(x): return x[0]+x[1]*1.6180339887498949

start=min(range(len(V)),key=lambda i:RAD[i])
MEAS=[5,10,20,30,40,50,60,70,80]; TMAX=80

def run(quantum):
    a=[0]*n
    for e in INC[start]: a[e]=1
    res={}; N0=sum(x*x for x in a); bad=0
    for t in range(1,TMAX+1):
        b=[0]*n
        for vi in range(len(V)):
            inc=INC[vi]; out=OUT[vi]; d=DEG[vi]
            S=0
            for e in inc: S+=a[e]
            if S==0 and all(a[e]==0 for e in inc): continue
            if quantum:
                c=COFF[d]*S
                for k in range(d): b[out[k]]=c-12*a[inc[k]]
            else:
                c=CP[d]*S
                for k in range(d): b[out[k]]=c
        a=b
        if quantum:
            N1=sum(x*x for x in a)
            if N1!=144*N0: bad+=1
            N0=N1
        if t in MEAS:
            S=(0,0); Q=(0,0); rmax=0.0
            vw=defaultdict(int)
            for i in range(n):
                if a[i]: vw[VI[E[i][1]]]+= a[i]*a[i] if quantum else a[i]
            for v,w in vw.items():
                S=za(S,(w,0)); Q=za(Q,zm((w,0),R2[v]))
                if RAD[v]>rmax: rmax=RAD[v]
            res[t]=(S,Q,rmax)
    return res,bad

t0=time.time(); A,bad=run(True);  print("装置A 済 %.0fs  二乗和が144倍からずれた歩数: %d"%(time.time()-t0,bad))
t0=time.time(); B,_  =run(False); print("古典B 済 %.0fs"%(time.time()-t0))
print()
print("  t | 装置A <r^2> | 古典B <r^2> | 前線半径")
for t in MEAS:
    SA,QA,ra=A[t]; SB,QB,rb=B[t]
    print(f"{t:3d} | {zf(QA)/zf(SA):11.2f} | {zf(QB)/zf(SB):11.2f} | {ra:7.1f}")
print()
print("倍率 t->2t （4なら弾道 t^2 / 2なら拡散 t）")
for t in (5,10,20,30,40):
    if 2*t not in A: continue
    for lab,M in (("A装置",A),("B古典",B)):
        S1,Q1,_=M[t]; S2,Q2,_=M[2*t]
        L=zm(Q2,S1); Rr=zm(Q1,S2)
        print(f"  t={t:2d}->{2*t:2d} {lab}: 倍率={zf(L)/zf(Rr):.3f}")
