# 平衡構造（re・θe）の記録の幅に入る、最小分母の組 (p,q)：cos = p/q
# 角の換算だけ浮動小数。幅は記録の誤差、または記録どうしのばらつき
import math
REC=[("H2O",104.4776-0.0019,104.4776+0.0019),
     ("NH3 (CCSD(T) 106.77 と Herzberg 106.67 の間)",106.67,106.77),
     ("NH3 (CCSD(T) 106.769 のみ ±0.01)",106.759,106.779),
     ("PH3 (実験 93.4 〜 ab initio 93.56)",93.4,93.56),
     ("H2S (92.11 ±0.01)",92.10,92.12)]
def pairs(lo,hi,qmax=60,k=4):
    c1,c2=sorted([math.cos(math.radians(lo)),math.cos(math.radians(hi))]); out=[]
    for q in range(1,qmax+1):
        for p in range(-q,1):
            if math.gcd(p,q)==1 and c1<=p/q<=c2: out.append((p,q))
        if len(out)>=k: break
    return out
for nm,lo,hi in REC:
    print(nm, "→", pairs(lo,hi))
print("模型 (-2,7)=",round(math.degrees(math.acos(-2/7)),4),"  (-1,4)=",round(math.degrees(math.acos(-1/4)),4))
