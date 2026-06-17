# -*- coding: utf-8 -*-
"""
宿題② の整数不変量 = 欠陥位置の「黄金体 Q(√5) 所属」。
Schmidt の定理（φ は Pisot）: ξ の φ進展開が（最終的に）周期的 ⟺ ξ ∈ Q(√5)。
 - pentagon（正20面体・∈Q(√5)）-> 周期的/有限 = 署名層（閉じる整数模様）
 - umbilic （Fib楕円体・∉Q(√5)）-> 非周期 = 計量層（無限・収束だが閉じない）
細分写像 p:2->6 は計量側の umbilic を署名側の pentagon 3個へ運ぶ（電荷 12=6χ 保存）。
"""
from mpmath import mp, mpf, sqrt, nint, fabs, isint
from fractions import Fraction
mp.dps = 60
PHI=(1+sqrt(5))/2

def in_Q5(r):
    """√r ∈ Q(√5)? r∈Q. (p+q√5)^2=p^2+5q^2+2pq√5=r が有理解を持つ ⟺ pq=0 かつ
       r/1 または r/5 が有理数の平方。"""
    r=Fraction(r).limit_denominator(10**12)
    def is_sq(fr):
        n,d=fr.numerator,fr.denominator
        if n<0: return False
        import math
        sn=math.isqrt(n); sd=math.isqrt(d)
        return sn*sn==n and sd*sd==d
    return is_sq(r) or is_sq(r/5)

def phi_digits(x,kmax=6,kmin=-60):
    x=mpf(x); out=[]
    for k in range(kmax,kmin-1,-1):
        w=PHI**k; d=max(-1,min(1,int(nint(x/w)))); out.append((k,d)); x=x-d*w
    return out, x

def classify(x, label):
    digs, res = phi_digits(x)
    frac=[d for k,d in digs if k<0]
    # 有限(末尾全0)?
    tail_zero = all(d==0 for d in frac[-30:])
    # 短周期検出（末尾40桁で最小周期）
    seq=frac[-40:]; period=None
    for p in range(1,13):
        if all(seq[i]==seq[i-p] for i in range(p,len(seq))):
            period=p; break
    state = "有限終端" if tail_zero else (f"周期{period}" if period else "非周期(窓内)")
    print(f"  {label:<26} 残差@φ^-60={mp.nstr(fabs(res),3):>10}  状態={state}")
    return state

print("=== 欠陥位置の黄金体所属（整数不変量）===")
A,B,C=21,13,8
ux2=Fraction(A*(A-B),(A-C))      # 168/13
uz2=Fraction(C*(B-C),(A-C))      # 40/13
print(f"  umbilic x^2={ux2} ,  ∈Q(√5)? {in_Q5(ux2)}   (168/13 も 168/65 も有理平方でない)")
print(f"  umbilic z^2={uz2} ,  ∈Q(√5)? {in_Q5(uz2)}")
print(f"  pentagon(icosa) 頂点座標 φ, 1 は自明に ∈Q(√5)\n")

print("=== Schmidt 周期性で metric/signature を判別 ===")
classify(sqrt(mpf(ux2.numerator)/ux2.denominator),"umbilic x (∉Q5,計量層)")
classify(sqrt(mpf(uz2.numerator)/uz2.denominator),"umbilic z (∉Q5,計量層)")
classify(PHI,                                      "pentagon φ (∈Q5,署名層)")
classify(mpf(1),                                   "pentagon 1 (∈Q5,署名層)")
classify(mpf(1)/3,                                 "1/3 (∈Q5: 周期の参照)")

print("\n=== 電荷整数（Sturm 署名側・保存量）===")
print(f"  umbilic: 1/6量子で 3, 個数4 -> 12 ;  pentagon: 1, 個数12 -> 12 ;  総和=6χ=12")
print(f"  写像: umbilic(電荷3,非周期φ位置) --細分3--> pentagon×3(各電荷1,周期φ位置)")
