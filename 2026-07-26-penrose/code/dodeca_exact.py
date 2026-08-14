#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
dodeca_exact.py ── 正12面体を Z[φ] の整数演算だけで組む
============================================================
一歩＝面の平面での鏡映。ひねりは要らない（面法線は5回対称軸なので、
鏡映一発で隣の正12面体が出て、面の5頂点がそのまま共有される）。

Z[φ] の元は (a,b) = a + bφ（φ² = φ + 1）。浮動小数は表示にしか使わない。
面法線 (0,±φ,±1) の巡回。n·n = 2+φ。面までの内積 = φ²。
鏡映 x → x − 2((x·n) − φ²)/(n·n) · n はすべて Z[φ] の中で閉じる。

事前登録（判別法6・27・28）── 数値を見る前に書いた
  EX0 独立に知れる値との照合（判別法28の系）
      外接／内接比 = 1.258409、二面角 = 116.5651°、面の頂点5個
      合わなければ、その先は走らせない
  EX1 中心1個＋12個は重なりなく組めるか
      OK なら：level 1 が成立
      NG なら：面で接する置き方では殻ができない
  EX2 2枚の鏡映の積 S_i S_j は有限位数か
      OK なら：その2つで輪ができる歩数がある
      NG なら：2種類の一歩を交互に打つ輪は存在しない（歩数によらず）
  EX3 長さ7までの全語（最初の一歩は対称性で固定）に閉じる輪はあるか
      OK なら：閉じる歩き方が実在する
      NG なら：長さ7以下には無い。長さ8以上は未探索と書く
"""
from fractions import Fraction as F
import itertools, math

# ── Z[φ] ──────────────────────────────────────────────────
def ad(x,y): return (x[0]+y[0], x[1]+y[1])
def sb(x,y): return (x[0]-y[0], x[1]-y[1])
def ml(x,y): return (x[0]*y[0]+x[1]*y[1], x[0]*y[1]+x[1]*y[0]+x[1]*y[1])
def sc(c,x): return (F(c)*x[0], F(c)*x[1])
def iv(x):
    a,b = x; d = a*a + a*b - b*b            # ノルム (a+bφ)(a+b(1−φ))
    if d == 0: raise ZeroDivisionError
    return ((a+b)/d, -b/d)
def sgn(x):
    a,b = x
    if b == 0: return (a>0)-(a<0)
    q = -a/b                                 # φ と q の大小
    t = 2*q-1
    gt = True if t < 0 else (5 > t*t)        # φ > q か
    return (1 if gt else -1) if b > 0 else (-1 if gt else 1)
def val(x): return float(x[0])+float(x[1])*(1+5**0.5)/2
Z0=(F(0),F(0)); Z1=(F(1),F(0)); PH=(F(0),F(1)); IP=sb(PH,Z1); P2=ad(PH,Z1)

def vadd(a,b): return tuple(ad(a[i],b[i]) for i in range(3))
def vsub(a,b): return tuple(sb(a[i],b[i]) for i in range(3))
def vsc(c,a): return tuple(ml(c,a[i]) for i in range(3))
def dot(a,b):
    r=Z0
    for i in range(3): r=ad(r,ml(a[i],b[i]))
    return r

# ── 正12面体 ──────────────────────────────────────────────
V=[]
for s1 in(1,-1):
    for s2 in(1,-1):
        for s3 in(1,-1): V.append((sc(s1,Z1),sc(s2,Z1),sc(s3,Z1)))
for s1 in(1,-1):
    for s2 in(1,-1):
        V+=[(Z0,sc(s1,IP),sc(s2,PH)),(sc(s1,IP),sc(s2,PH),Z0),(sc(s1,PH),Z0,sc(s2,IP))]
N=[]
for s1 in(1,-1):
    for s2 in(1,-1):
        N+=[(Z0,sc(s1,PH),sc(s2,Z1)),(sc(s1,PH),sc(s2,Z1),Z0),(sc(s1,Z1),Z0,sc(s2,PH))]
NN=dot(N[0],N[0]); OFF=P2                      # n·n = 2+φ ／ 面までの内積 = φ²
INN=iv(NN)

print("EX0 独立に知れる値との照合")
rin=val(OFF)/val(NN)**0.5; rout=max(sum(val(c)**2 for c in v)**0.5 for v in V)
ang=sorted({round(math.degrees(math.acos(max(-1,min(1,val(dot(N[i],N[j]))/val(NN))))),4)
            for i,j in itertools.combinations(range(12),2)})
nf=[sum(1 for v in V if dot(v,n)==OFF) for n in N]
print("   内接 %.6f 外接 %.6f 比 %.6f（1.258409）／ 二面角 %.4f°（116.5651）／ 面の頂点 %s"
      % (rin,rout,rout/rin,180-ang[0],sorted(set(nf))))
assert abs(rout/rin-1.258409)<1e-5 and abs((180-ang[0])-116.5651)<1e-3 and set(nf)=={5}
print("   → OK（照合が通ったので先へ進む）")

# ── 一歩＝鏡映（アフィン変換 (M,t): x → Mx + t） ──────────────
I3=((Z1,Z0,Z0),(Z0,Z1,Z0),(Z0,Z0,Z1)); T0=(Z0,Z0,Z0)
def refl_i(i):
    n=N[i]; M=[]
    for r in range(3):
        M.append(tuple(sb(I3[r][c], ml(sc(2,ml(n[r],n[c])),INN)) for c in range(3)))
    t=vsc(ml(sc(2,OFF),INN), n)
    return (tuple(M), t)
S=[refl_i(i) for i in range(12)]
def mv(M,p): return tuple(ad(ad(ml(M[i][0],p[0]),ml(M[i][1],p[1])),ml(M[i][2],p[2])) for i in range(3))
def comp(A,B):
    MA,tA=A; MB,tB=B
    M=tuple(tuple(ad(ad(ml(MA[i][0],MB[0][j]),ml(MA[i][1],MB[1][j])),ml(MA[i][2],MB[2][j]))
                  for j in range(3)) for i in range(3))
    return (M, vadd(mv(MA,tB),tA))
ID=(I3,T0)
def is_id(T): return T[0]==I3 and T[1]==T0
def verts(T): return frozenset(vadd(mv(T[0],v),T[1]) for v in V)

print("\nEX1 中心1個＋12個（各面で鏡映）")
base=verts(ID); kids=[verts(S[i]) for i in range(12)]
print("   面を共有する頂点数 %s（5なら面で接している）"
      % sorted({len(base & k) for k in kids}))
E=[(a,b) for a,b in itertools.combinations(range(20),2)
   if val(dot(vsub(V[a],V[b]),vsub(V[a],V[b])))<1.6]
ED=[]
for a,b in E:
    d=vsub(V[b],V[a])
    if not any(all(sb(ml(d[i],e[j]),ml(d[j],e[i]))==Z0 for i in range(3) for j in range(3)) for e in ED):
        ED.append(d)
def sat(TA,TB):
    PA=[vadd(mv(TA[0],v),TA[1]) for v in V]; PB=[vadd(mv(TB[0],v),TB[1]) for v in V]
    axes=[mv(TA[0],N[i]) for i in range(6)]+[mv(TB[0],N[i]) for i in range(6)]
    for u in [mv(TA[0],e) for e in ED]:
        for w in [mv(TB[0],e) for e in ED]:
            x=(sb(ml(u[1],w[2]),ml(u[2],w[1])),sb(ml(u[2],w[0]),ml(u[0],w[2])),
               sb(ml(u[0],w[1]),ml(u[1],w[0])))
            if any(z!=Z0 for z in x): axes.append(x)
    for a in axes:
        p=[dot(q,a) for q in PA]; r=[dot(q,a) for q in PB]
        mp,Mp=min(p,key=val),max(p,key=val); mr,Mr=min(r,key=val),max(r,key=val)
        if sgn(sb(mr,Mp))>=0 or sgn(sb(mp,Mr))>=0: return False
    return True
bad=sum(1 for a,b in itertools.combinations(range(12),2) if sat(S[a],S[b]))
badc=sum(1 for i in range(12) if sat(S[i],ID))
print("   独立な辺方向 %d ／ 軸 %d 本で判定" % (len(ED),12+len(ED)**2))
print("   子どうしの重なり %d / 66 ／ 子と中心の重なり %d / 12 → %s"
      % (bad,badc,"OK" if bad==0 and badc==0 else "NG"))

print("\nEX2 2枚の鏡映の積の位数")
per=set()
for j in range(1,12):
    T=comp(S[0],S[j]); k=1; U=T
    while k<=60 and not is_id(U): U=comp(U,T); k+=1
    per.add(k if k<=60 else None)
print("   S_0·S_j の位数（60まで探索）: %s → %s"
      % (sorted(x for x in per if x) or "見つからず",
         "NG（交互打ちの輪は存在しない）" if per=={None} else "OK"))
c=val(dot(N[0],N[1]))/val(NN)
print("   面のなす角 %.4f° ／ 積は %.4f° の回転 ／ 360/その値 = %.4f（整数でない）"
      % (math.degrees(math.acos(c)),2*math.degrees(math.acos(c)),
         360/(2*math.degrees(math.acos(c)))))

print("\nEX3 長さ7までの全語（最初の一歩は S_0 に固定）")
found=[]
def dfs(T,last,depth,limit,path):
    if depth>0 and is_id(T): found.append(tuple(path)); return
    if depth==limit: return
    d=math.hypot(*[val(x) for x in T[1]])
    if d > (limit-depth)*2*rin+1e-9: return
    for i in range(12):
        if i==last: continue
        dfs(comp(T,S[i]),i,depth+1,limit,path+[i])
for L in (4,5,6,7):
    found.clear(); dfs(comp(ID,S[0]),0,1,L,[0])
    print("   長さ%d まで: 閉じた語 %d 本" % (L,len(found)))
print("   → %s" % ("NG（長さ7以下に閉じる輪は無い。長さ8以上は未探索）"
                   if not found else "OK"))
