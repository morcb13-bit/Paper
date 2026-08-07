"""
「全部1ではないか」を検定する。担体も窓も使わない。Z[zeta10] の整数演算のみ。

#  検定ON1 四つの一歩の絶対ノルム
#      OK なら：飛・連・奇数の連2本は絶対ノルムが等しい＝同じ一つの元の同伴でありうる
#      NG なら：別々の元であり、「全部1」は成り立たない
#
#  検定ON2 互いの商が単数か
#      OK なら：四つは同伴。実際に一つの元の単数倍で、見え方の違いは単数の作用
#      NG なら：ノルムが同じでも別の元。同伴ではない
#
#  検定ON3 商の中身（どの単数か）
#      OK なら：商は phi と zeta だけで書ける。phi が長さの単数、zeta が向きの単数
#      NG なら：他の単数が要る。話が phi と zeta で閉じない
#
#  検定ON4 五芒星から出る一歩のノルム
#      OK なら：ノルム1＝文字どおり単数。五芒星の一歩は「1」そのもの
#      NG なら：単数でない
#
#  検定ON5 1になりきらない残りは何か
#      OK なら：円環をつなぐ一歩は (1-zeta^2) の同伴として書ける
#      NG なら：残りは (1-zeta^2) では書けない
#      ※ 初版は 1-zeta を取った。10 が素数冪でないので 1-zeta は可逆量（N=1）で、
#         これは誤り。分岐する側は 1-zeta^2（= 1 - zeta5）。検定が NG を返して判明した。
#
#  検定ON6 残渣で偶数の鎖が往復する理由
#      OK なら：phi の共役 psi が負であること「だけ」で説明が尽きる
#      NG なら：符号以外の要素が要る
"""
from fractions import Fraction as F
import math, cmath
from qphi import Qp, zmul, zconj, zsigma, zre, PHI, ONE

def zt(k):
    k %= 10
    base = [(1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1),(-1,1,-1,1)]
    if k < 5: return tuple(F(x) for x in base[k])
    return tuple(-F(x) for x in base[k-5])

def zadd(a,b): return tuple(a[i]+b[i] for i in range(4))
def zsub(a,b): return tuple(a[i]-b[i] for i in range(4))
def zrot(a,k): return zmul(a, zt(k))
def zabs2(a): return zre(zmul(a, zconj(a)))
def res2(a):  return zabs2(zsigma(a))

def znorm(a):
    """絶対ノルム N(a) = |a|^2 * sigma(|a|^2)。有理数で返る。"""
    n = zabs2(a); s = Qp(n.p+n.q, -n.q)
    r = n*s
    assert r.q == 0, "ノルムが有理数にならない"
    return r.p

def zinv_times_norm(a):
    """conj(a)*sigma(a)*sigma(conj(a))。a を掛けると N(a) になる。"""
    return zmul(zmul(zconj(a), zsigma(a)), zsigma(zconj(a)))

def zdiv(a, b):
    """a/b。整数係数にならなければ None。"""
    n = znorm(b)
    if n == 0: return None
    num = zmul(a, zinv_times_norm(b))
    q = tuple(x/n for x in num)
    if any(x.denominator != 1 for x in q): return None
    return q

PHI_Z = zadd(zt(1), zt(9))          # zeta + zeta^-1 = 2cos36 = phi
PHI2_Z = zmul(PHI_Z, PHI_Z)
CONT = [zrot(zmul(PHI2_Z, zsub(zt(4), zt(0))), k) for k in range(10)]
SKIP = [zrot(zmul(PHI2_Z, zsub(zt(2), zt(0))), k) for k in range(10)]
C_A, C_B = CONT[5], CONT[7]
S_AX, C_AX = SKIP[7], CONT[6]

def num(a):
    z = cmath.exp(1j*math.pi/5)
    return sum(float(a[k])*z**k for k in range(4))
def ang(a): return math.degrees(cmath.phase(num(a))) % 360
def show(a): return "(" + ",".join(str(int(x)) for x in a) + ")"

STEPS = [("飛 S_AX", S_AX), ("連 C_AX", C_AX), ("奇数の連 C_A", C_A), ("奇数の連 C_B", C_B)]

# ---------------- ON1 ----------------
print("検定ON1 四つの一歩の絶対ノルム")
norms = []
for nm, v in STEPS:
    N = znorm(v); norms.append(N)
    print(f"  {nm:14s} {show(v):16s} |v|^2={zabs2(v)}  残渣^2={res2(v)}  N={N}")
on1 = len(set(norms)) == 1
print(f"検定ON1: {'OK' if on1 else 'NG'}  すべて等しいか = {on1}（値 {sorted(set(norms))}）")

# ---------------- ON2 ----------------
print("\n検定ON2 互いの商が単数か")
on2 = True
for i in range(4):
    for j in range(4):
        if i == j: continue
        q = zdiv(STEPS[i][1], STEPS[j][1])
        if q is None or abs(znorm(q)) != 1:
            on2 = False
            print(f"  {STEPS[i][0]} / {STEPS[j][0]} → 単数でない")
print(f"  商 連/飛 = {show(zdiv(C_AX, S_AX))}  N={znorm(zdiv(C_AX,S_AX))}")
print(f"  商 C_A/C_AX = {show(zdiv(C_A, C_AX))}  N={znorm(zdiv(C_A,C_AX))}")
print(f"  商 C_B/C_AX = {show(zdiv(C_B, C_AX))}  N={znorm(zdiv(C_B,C_AX))}")
print(f"検定ON2: {'OK' if on2 else 'NG'}  12組すべて単数か = {on2}")

# ---------------- ON3 ----------------
print("\n検定ON3 商の中身")
on3 = True
for nm, v in STEPS:
    q = zdiv(v, S_AX)
    found = None
    for e in range(-3, 4):
        pw = (F(1),F(0),F(0),F(0))
        for _ in range(abs(e)): pw = zmul(pw, PHI_Z)
        if e < 0:
            inv = zdiv((F(1),F(0),F(0),F(0)), pw)
            pw = inv if inv is not None else None
        if pw is None: continue
        for k in range(10):
            if zmul(pw, zt(k)) == q: found = (e, k); break
        if found: break
    print(f"  {nm:14s} / 飛 = {show(q):14s}  = phi^{found[0]} * zeta^{found[1]}" if found
          else f"  {nm:14s} / 飛 = {show(q):14s}  = phi^e * zeta^k では書けない")
    if not found: on3 = False
print(f"検定ON3: {'OK' if on3 else 'NG'}  四つとも phi^e * zeta^k で書けるか = {on3}")

# ---------------- ON4 ----------------
print("\n検定ON4 五芒星から出る一歩")
p3 = zmul(PHI_Z, PHI2_Z); p5 = zmul(p3, PHI2_Z)
on4 = True
for nm, v in (("phi^3 五芒星→円環", p3), ("phi^5 外側の殻", p5), ("phi", PHI_Z)):
    N = znorm(v)
    print(f"  {nm:18s} {show(v):16s} 実長={math.sqrt(zabs2(v).val()):9.6f}"
          f" 残渣長={math.sqrt(res2(v).val()):9.6f}  N={N}")
    if abs(N) != 1: on4 = False
print(f"検定ON4: {'OK' if on4 else 'NG'}  すべて単数（N=±1）か = {on4}")

# ---------------- ON5 ----------------
print("\n検定ON5 1になりきらない残り")
for nm, x in (("1-zeta   （zeta10 側）", zsub(zt(0), zt(1))),
              ("1-zeta^2 （zeta5 側）", zsub(zt(0), zt(2)))):
    print(f"  {nm} = {show(x):14s} N={znorm(x)}  実長={math.sqrt(zabs2(x).val()):.6f}"
          f"  残渣長={math.sqrt(res2(x).val()):.6f}")
lam = zsub(zt(0), zt(2))    # 1 - zeta^2
q = zdiv(S_AX, lam)
five = zmul(zmul(lam,lam), zmul(lam,lam))
q5 = zdiv((F(5),F(0),F(0),F(0)), five)
on5 = (q is not None and abs(znorm(q)) == 1 and q5 is not None and abs(znorm(q5)) == 1)
print(f"  飛 / (1-zeta^2)  = {show(q) if q else 'None'}  N={znorm(q) if q else '-'}")
print(f"  5  / (1-zeta^2)^4 = {show(q5) if q5 else 'None'}  N={znorm(q5) if q5 else '-'}")
print(f"検定ON5: {'OK' if on5 else 'NG'}  円環の一歩は (1-zeta^2) の同伴か = {on5}")

# ---------------- ON6 ----------------
print("\n検定ON6 残渣で往復する理由")
psi = zsigma(PHI_Z)
print(f"  phi = {show(PHI_Z)}  値={zre(PHI_Z).val():.6f}")
print(f"  sigma(phi) = psi = {show(psi)}  値={zre(psi).val():.6f}")
print(f"  連 = 飛 x phi なので、残渣では 連* = 飛* x psi")
print(f"  飛*  向き={ang(zsigma(S_AX)):7.3f}deg 長さ={math.sqrt(res2(S_AX).val()):.6f}")
print(f"  連*  向き={ang(zsigma(C_AX)):7.3f}deg 長さ={math.sqrt(res2(C_AX).val()):.6f}")
on6 = (zdiv(C_AX, S_AX) == PHI_Z) and (zre(psi).val() < 0)
print(f"  連/飛 が phi ちょうどか = {zdiv(C_AX,S_AX)==PHI_Z}／psi が負か = {zre(psi).val()<0}")
print(f"検定ON6: {'OK' if on6 else 'NG'}")

# ---------------- 鎖を一つの元で書き直す ----------------
print("\n鎖を『飛 x 単数』だけで書く")
for n in (9, 10):
    v = [C_A, C_B] if n%2 else [S_AX, C_AX]
    ws = []
    for s in v:
        q = zdiv(s, S_AX)
        for e in range(0,3):
            pw = (F(1),F(0),F(0),F(0))
            for _ in range(e): pw = zmul(pw, PHI_Z)
            hit = [k for k in range(10) if zmul(pw, zt(k)) == q]
            if hit: ws.append(f"phi^{e}*zeta^{hit[0]}"); break
    print(f"  n={n}（{'奇数' if n%2 else '偶数'}）: 飛 x [{ws[0]}] と 飛 x [{ws[1]}] の交互")
