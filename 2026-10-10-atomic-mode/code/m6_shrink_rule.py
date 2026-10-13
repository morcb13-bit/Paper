# 検定3：電子を H に置き換えると腕が縮む規則を、未調査の元素（N殻・O殻）まで当てる
# 規則（事前に書いたもの）：手どうしの内積 -2、手の長さ² = 6 + d·(4-h)  （h = H の数、L殻は d=1 → 10-h）
#   ※ L殻 B側（CH4・NH3・H2O）は規則を作った材料なので検定ではなく「当てはめ」
#   M殻以降で外れたら、d を殻の大きさ n² に比例させる（L:M:N:O = 4:9:16:25）一つだけを足す
#   A側（h<4）は「A の席が空くたびに縮む」：長さ² = 6 - (4-h)
# 負の対照：縮まない（長さ² 10 固定）、縮み方を逆（長さ² = 6 - (4-h) を B側にも）
# 判定は整数の交差掛け算だけ。浮動小数は実測角を組に換算する一回と、必要な d の窓を出す一回だけ
import math
same=lambda a,b:a[0]*b[1]==b[0]*a[1]
cross=lambda a,b:abs(a[0]*b[1]-b[0]*a[1])
def to_pair(deg,tol=0.1,qmax=200):
    lo,hi=sorted([math.cos(math.radians(deg-tol)),math.cos(math.radians(min(deg+tol,180)))])
    for q in range(1,qmax+1):
        for p in range(-q,q+1):
            if lo<=p/q<=hi: return (p,q)
def d_window(deg,h,tol=0.1):          # 内積-2・長さ² 6+d(4-h) で角が ±tol に入る整数 d
    if h>=4: return None
    out=[]
    for d in range(0,400):
        c=-2/(6+d*(4-h)); a=math.degrees(math.acos(c))
        if abs(a-deg)<=tol: out.append(d)
    return out
# 実測（[知識] CCCBDB 気相の実験値。GeH4・SnH4 は正四面体、GaH3・AlH3・BH3 は平面で 120、BeH2・MgH2 は直線）
DATA=[("L",2,"BeH2",2,180.0),("L",2,"BH3",3,120.0),("L",2,"CH4",4,109.47),("L",2,"NH3",3,106.67),("L",2,"H2O",2,104.48),
      ("M",3,"MgH2",2,180.0),("M",3,"AlH3",3,120.0),("M",3,"SiH4",4,109.47),("M",3,"PH3",3,93.35),("M",3,"H2S",2,92.11),
      ("N",4,"GaH3",3,120.0),("N",4,"GeH4",4,109.47),("N",4,"AsH3",3,92.08),("N",4,"H2Se",2,90.90),
      ("O",5,"SnH4",4,109.47),("O",5,"SbH3",3,91.56),("O",5,"H2Te",2,90.26)]
Aside={"BeH2","BH3","MgH2","AlH3","GaH3"}
def rule(nm,n,h):
    if nm in Aside: return (-2,6-(4-h))
    return (-2,6+1*(4-h))
def rule_n2(nm,n,h):                   # d ∝ n²（L を 4 として 4倍の整数で）：内積 -8、長さ² 24 + n²(4-h)
    if nm in Aside: return (-8,24-n*n*(4-h))
    return (-8,24+n*n*(4-h))
ctrl_fix=lambda nm,n,h:(-2,10)
ctrl_rev=lambda nm,n,h:(-2,6-(4-h))
print("殻 分子  H  実測角   実測の組   規則      判定 差 | n²版      判定 差 | 固定 反転 | 合うd")
tot={"r":0,"n2":0,"fix":0,"rev":0}
for sh,n,nm,h,deg in DATA:
    m=to_pair(deg); r=rule(nm,n,h); r2=rule_n2(nm,n,h)
    j=lambda x:"OK" if same(x,m) else "NG"
    for k,x in (("r",r),("n2",r2),("fix",ctrl_fix(nm,n,h)),("rev",ctrl_rev(nm,n,h))):
        tot[k]+=same(x,m)
    dw=None if nm in Aside else d_window(deg,h)
    print(f"{sh}  {nm:5} {h}  {deg:7}  {str(m):10} {str(r):9} {j(r)} {cross(r,m):3} | {str(r2):9} {j(r2)} {cross(r2,m):3} | {j(ctrl_fix(nm,n,h))}   {j(ctrl_rev(nm,n,h))}  | {dw if dw is not None else '―'}")
print("OK の数（17 分子中）：", tot)
