#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
dodeca_rings.py ── 正12面体の鎖が閉じるか（Z[φ] の整数演算・中間一致法）
============================================================
一歩＝面の平面での鏡映。同次座標で (2+φ) 倍して整数化すると、
長さ L の語の変換は Z[φ] 係数の 4×4 整数行列になる（分母なし）。

  H_i = [[ (2+φ)I − 2 n nᵀ , 2φ² n ],
         [        0        ,  2+φ  ]]      すべて Z[φ] の整数

  長さ L の語が閉じる ⟺ 積 = (2+φ)^L · I₄

事前登録（判別法6・27）── 数値を見る前に書いた
  R0 奇数長の語は閉じない（鏡映の個数が奇数なら向きが反転する）
     これは算術であって測定ではない。したがって偶数長だけ探す
  R1 長さ8・10・12 に閉じる語はあるか（後戻り禁止・最初の一歩は対称性で固定）
     OK なら：閉じる歩き方が実在する。監督の予想（12）が当たる
     NG なら：長さ12以下に輪は無い
  対照：H_i·H_i = (2+φ)²·I（鏡映は対合）が成り立つこと。
        成り立たなければ整数化そのものが誤り
"""
import itertools, sys

# ── Z[φ]（整数係数） a + bφ ────────────────────────────────
def ad(x,y): return (x[0]+y[0], x[1]+y[1])
def sb(x,y): return (x[0]-y[0], x[1]-y[1])
def ml(x,y): return (x[0]*y[0]+x[1]*y[1], x[0]*y[1]+x[1]*y[0]+x[1]*y[1])
def val(x): return x[0]+x[1]*(1+5**0.5)/2
Z0=(0,0); Z1=(1,0); PH=(0,1); IP=(-1,1); P2=(1,1); S2=(2,1)   # 1/φ=φ−1, φ²=1+φ, 2+φ
def divc(x):
    """(a+bφ)/(2+φ) = (a+bφ)(3−φ)/5 — 割り切れる前提"""
    y=ml(x,(3,-1))
    assert y[0]%5==0 and y[1]%5==0, ("割り切れない",x)
    return (y[0]//5, y[1]//5)

N=[]
for s1 in(1,-1):
    for s2 in(1,-1):
        N+=[(Z0,(s1*0,s1*1),(s2*1,0)),((s1*0,s1*1),(s2*1,0),Z0),((s1*1,0),Z0,(s2*0,s2*1))]
# 面法線 (0,±φ,±1) の巡回
def dot(a,b):
    r=Z0
    for i in range(3): r=ad(r,ml(a[i],b[i]))
    return r
assert all(dot(n,n)==S2 for n in N), "法線の長さ² が 2+φ でない"

def Hmat(i):
    n=N[i]; R=[]
    for r in range(3):
        R.append(tuple(sb(S2 if r==c else Z0, ml((2,0),ml(n[r],n[c]))) for c in range(3)))
    v=tuple(ml(ml((2,0),P2),n[r]) for r in range(3))
    return (tuple(R), v)
H=[Hmat(i) for i in range(12)]

# 対照：鏡映は対合か
R0,v0=H[0]
A=tuple(tuple(ad(ad(ml(R0[i][0],R0[0][j]),ml(R0[i][1],R0[1][j])),ml(R0[i][2],R0[2][j]))
              for j in range(3)) for i in range(3))
b=tuple(ad(ad(ad(ml(R0[i][0],v0[0]),ml(R0[i][1],v0[1])),ml(R0[i][2],v0[2])),ml(S2,v0[i]))
        for i in range(3))
ok = all(A[i][j]==(ml(S2,S2) if i==j else Z0) for i in range(3) for j in range(3)) and all(x==Z0 for x in b)
print("対照 H_i·H_i = (2+φ)²·I → %s" % ("OK" if ok else "NG"))
assert ok

def step(P,i):
    """P = (A|b) の左から H_i を掛ける。cs は (2+φ)^depth"""
    (Ap,bp,c)=P; (R,v)=H[i]
    A=tuple(tuple(ad(ad(ml(R[r][0],Ap[0][j]),ml(R[r][1],Ap[1][j])),ml(R[r][2],Ap[2][j]))
                  for j in range(3)) for r in range(3))
    b=tuple(ad(ad(ad(ml(R[r][0],bp[0]),ml(R[r][1],bp[1])),ml(R[r][2],bp[2])),ml(v[r],c))
            for r in range(3))
    return (A,b,ml(c,S2))
ID=(((Z1,Z0,Z0),(Z0,Z1,Z0),(Z0,Z0,Z1)),(Z0,Z0,Z0),Z1)

def words(half, first=None):
    """長さ half の後戻り無しの語を (最初の面, 最後の面, 状態) で返す"""
    out=[]
    def rec(P,last,d,f):
        if d==half: out.append((f,last,P)); return
        for i in range(12):
            if i==last: continue
            rec(step(P,i),i,d+1,f if f is not None else i)
        return
    if first is None: rec(ID,None,0,None)
    else: rec(step(ID,first),first,1,first)
    return out

def target(P):
    """P と組んで閉じる相手の (A,b)"""
    (Ap,bp,c)=P
    At=tuple(tuple(Ap[j][i] for j in range(3)) for i in range(3))
    Ab=tuple(ad(ad(ml(At[r][0],bp[0]),ml(At[r][1],bp[1])),ml(At[r][2],bp[2])) for r in range(3))
    # b_Q = −Aᵀb / c ，c = (2+φ)^half
    x=tuple((-Ab[r][0],-Ab[r][1]) for r in range(3))
    k=0; cc=c
    while cc!=Z1: cc=divc(cc); k+=1
    for _ in range(k): x=tuple(divc(y) for y in x)
    return (At,x)

def scan(half, D):
    hits=[0]
    def rec(P,last,d,f):
        if d==half:
            if last!=0:
                c=D.get(((P[0],P[1]),f))
                if c: hits[0]+=len(c)
            return
        for i in range(12):
            if i==last: continue
            rec(step(P,i),i,d+1,f if f is not None else i)
    rec(ID,None,0,None)
    return hits[0]

for L in (10,12):
    half=L//2
    P=words(half, first=0)
    D={}
    for f,last,st in P: D.setdefault((target(st),last),[]).append((f,last))
    n=scan(half,D)
    print("長さ%2d : 前半 %7d 語 ／ 閉じた語 %d" % (L,len(P),n))
    sys.stdout.flush()
