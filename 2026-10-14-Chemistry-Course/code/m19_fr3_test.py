# m19: 検定 FR3（事前登録 FR3_prereg.md のとおり）。判定は整数の比較だけ。
import random
V={}
for e in "Li Na K Rb Cs".split(): V[e]=1
for e in "Be Mg Ca Sr Ba".split(): V[e]=2
for e in "B Al Ga In Tl".split(): V[e]=3
for e in "C Si Ge Sn Pb".split(): V[e]=4
for e in "N P As Sb Bi".split(): V[e]=5
for e in "O S Se Te Po".split(): V[e]=6
for e in "F Cl Br I At".split(): V[e]=7
XSET=set("Ga Ge As Se Br In Sn Sb Te I Tl Pb Bi Po At".split())

def side_fr(e,X): return "D" if V[e]<=3 else "S"
def side_fr3(e,X):
    if e in X: return "D" if V[e]+V[e]<=8 else "S"
    return "D" if V[e]<=3 else "S"
def elem(e,X,fr3):
    return "金属" if (side_fr3 if fr3 else side_fr)(e,X)=="D" else "共有"
def comp(a,b,X,fr3):
    sd=side_fr3 if fr3 else side_fr
    p,q=sd(a,X),sd(b,X)
    if p==q=="D": return "金属"
    if p==q=="S": return "共有"
    if fr3:
        d=a if p=="D" else b
        if d in X: return "共有"
    return "イオン"

# 照合 1（昨日の 52）
ELEM={}
for e in "Li Na K Rb Cs Be Mg Ca Sr Ba Al Ga In Tl Sn Pb Bi".split(): ELEM[e]="金属"
for e in "C N P O S Se F Cl Br I".split(): ELEM[e]="共有"
C1=[("NaCl","Na","Cl","イオン"),("KBr","K","Br","イオン"),("LiF","Li","F","イオン"),("MgO","Mg","O","イオン"),
 ("CaF2","Ca","F","イオン"),("CaO","Ca","O","イオン"),("Na2O","Na","O","イオン"),("K2S","K","S","イオン"),
 ("MgCl2","Mg","Cl","イオン"),("BaCl2","Ba","Cl","イオン"),("CO2","C","O","共有"),("CS2","C","S","共有"),
 ("CCl4","C","Cl","共有"),("SiO2","Si","O","共有"),("SiC","Si","C","共有"),("NCl3","N","Cl","共有"),
 ("PCl3","P","Cl","共有"),("SO2","S","O","共有"),("SF6","S","F","共有"),("OF2","O","F","共有"),
 ("BN","B","N","共有"),("GaAs","Ga","As","共有"),("NaK","Na","K","金属"),("CaMg2","Ca","Mg","金属"),("Mg17Al12","Mg","Al","金属")]
def score1(fr3):
    s=0;miss=[]
    for e,r in ELEM.items():
        p=elem(e,XSET,fr3); s+=p==r
        if p!=r: miss.append(e)
    for n,a,b,r in C1:
        p=comp(a,b,XSET,fr3); s+=p==r
        if p!=r: miss.append(n)
    return s,miss
for fr3 in (False,True):
    s,m=score1(fr3); print("照合1",("FR3" if fr3 else "FR "),f"{s}/52 外れ",m)

# 照合 2（新しい 25）
C2=[("GaP","Ga","P"),("GaN","Ga","N"),("InP","In","P"),("InAs","In","As"),("InSb","In","Sb"),("AlP","Al","P"),
 ("AlAs","Al","As"),("AlSb","Al","Sb"),("BP","B","P"),("GaCl3","Ga","Cl"),("SnCl4","Sn","Cl"),("GeCl4","Ge","Cl"),
 ("SbCl3","Sb","Cl"),("AsCl3","As","Cl"),("SiCl4","Si","Cl"),("PBr3","P","Br")]
C2=[c+("共有",) for c in C2]+[(n,a,b,"イオン") for n,a,b in [("KI","K","I"),("CsI","Cs","I"),("RbBr","Rb","Br"),
 ("NaI","Na","I"),("NaBr","Na","Br"),("LiI","Li","I"),("SrBr2","Sr","Br"),("BaI2","Ba","I"),("CsBr","Cs","Br")]]
def score2(X,fr3,show=False):
    s=0
    for n,a,b,r in C2:
        p=comp(a,b,X,fr3); s+=p==r
        if show and p!=r: print(f"   外れ {n}: 予言 {p} / 記録 {r}")
    return s
print("== 照合2 FR"); F=score2(XSET,False,True); print(f"  FR {F}/25")
print("== 照合2 FR3"); F3=score2(XSET,True,True); print(f"  FR3 {F3}/25")
els=sorted({x for c in C2 for x in c[1:3]}); k=sum(e in XSET for e in els)
rng=random.Random(13); T=10000; tot=0; ge=0; lo=99; hi=-1; hist={}
for _ in range(T):
    X=set(rng.sample(els,k)); s=score2(X,True); tot+=s; ge+=s>=F3; lo=min(lo,s); hi=max(hi,s); hist[s]=hist.get(s,0)+1
print(f"== 負の対照（X {k} 個を {len(els)} 元素にくじ引き）総和 {tot} 平均≈{tot//T}.{(tot%T)*100//T:02d} 最小 {lo} 最大 {hi}")
print("  分布:"," ".join(f"{a}:{hist[a]}" for a in sorted(hist)),f" FR3以上 {ge}/{T}")
c1=F3>=F+3; c2=T*F3>=tot+3*T
print(f"== 判定  {F3} ≥ {F}+3 → {c1}   {T*F3} ≥ {tot}+{3*T}={tot+3*T} → {c2}   FR3: {'OK' if c1 and c2 else 'NG'}")
