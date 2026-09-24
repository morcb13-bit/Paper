# 道具1：置き換え（ロビンソン三角形）を Z[ζ10] の整数で持ち、毎段 φ 倍に育てて切る
# 判定は整数と加算だけ（norm2 は p+qφ の整数対、面積は φ の整数係数で数える）
import sys, b13_chain_units as U
Z=U.zadd; S=U.zsub; M=U.zmul; PHI=U.PHI
def sc(z): return M(z,PHI)
def step(T, swap=False):
    out=[]
    for col,A,B,C in T:
        a,b,c=sc(A),sc(B),sc(C)
        if (col==0)!=swap:
            P=Z(a,S(B,A))                      # φA+(B−A)
            out+=[(0,c,P,b),(1,P,c,a)]
        else:
            Q=Z(b,S(A,B)); R=Z(b,S(C,B))
            out+=[(1,R,c,a),(1,Q,R,b),(0,R,Q,a)]
    return out
def sides(t):
    _,A,B,C=t; return tuple(U.norm2(S(x,y)) for x,y in ((A,B),(A,C),(B,C)))
def F(n):
    a,b=0,1
    if n<0: return (-1)**(n+1)*F(-n)
    for _ in range(n): a,b=b,a+b
    return a
def fadd(x,y): return (x[0]+y[0],x[1]+y[1])
def fmul(x,y): return (x[0]*y[0]+x[1]*y[1], x[0]*y[1]+x[1]*y[0]+x[1]*y[1])   # φ²=φ+1
def area(T):                         # 赤＝1、青＝φ（同じ段の脚の長さで）
    s=(0,0)
    for t in T: s=fadd(s,(1,0) if t[0]==0 else (0,1))
    return s
def check(seed, N, swap=False):
    T=seed; rows=[]
    for n in range(N+1):
        r=sum(1 for t in T if t[0]==0); b=len(T)-r
        shapes={t[0]:set() for t in T}
        for t in T: shapes[t[0]].add(tuple(sorted(sides(t))))
        rows.append((n,r,b,{k:len(v) for k,v in shapes.items()},area(T)))
        if n<N: T=step(T,swap)
    return rows,T
O=U.ZERO; e=U.zt
red=[(0,O,e(0),e(1))]; blue=[(1,O,e(0),U.zmul(e(3),(1,0,0,0)))]
N=12
print("種：赤一枚（36°の三角）"); ok=True
rows,_=check(red,N)
for n,r,b,sh,a in rows:
    exp=(F(2*n-1),F(2*n)); g=(r,b)==exp and all(v==1 for v in sh.values())
    ok&=g; print(n,r,b,"期待",exp,"形の種類",sh,"面積",a,"OK" if g else "NG")
# 面積：一段で φ² 倍（整数対で）
ar=[x[4] for x in rows]; ga=all(ar[i+1]==fmul(ar[i],(1,1)) for i in range(N))
print("面積が一段ごとに φ² 倍:", ga); ok&=ga
print("種：赤10枚の輪（中心のまわり）")
wheel=[(0,O,e(i),e(i+1)) if i%2==0 else (0,O,e(i+1),e(i)) for i in range(10)]
rows,_=check(wheel,8)
for n,r,b,sh,a in rows: print(n,r,b,"= 10×",(F(2*n-1),F(2*n)), "OK" if (r,b)==(10*F(2*n-1),10*F(2*n)) else "NG", sh)
print("負の対照：赤と青の規則を入れ替える")
rows,_=check(red,6,swap=True)
for n,r,b,sh,a in rows: print(n,r,b,"期待",(F(2*n-1),F(2*n)),"形の種類",sh,"面積",a)
print("総合", "OK" if ok else "NG")
