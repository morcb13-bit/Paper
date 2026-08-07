"""
奇数の鎖（ジグザグ）と偶数の鎖（直線）が同じ直線上に並ぶ理由を、
実空間と残渣の両方で見る。担体も窓も使わない。すべて Z[zeta10] の整数演算。

#  検定CD1 二歩ぶんの和
#      OK なら：奇数の一往復と偶数の飛+連は、格子ベクトルとして同一。
#               同じ直線に並ぶのは作図の偶然ではなく、二歩の和が同じ一歩だから
#      NG なら：一致は近似にすぎず、共線の理由は別にある
#
#  検定CD2 中間の環の位置（実空間と残渣）
#      OK なら：実空間で軸に近いほうが、残渣では原点から遠い（順序が逆転する）。
#               監督の「見かけの視角は同じだが実は遠方」に対応する量が存在する
#      NG なら：どちらの座標でも同じ側が近い。遠近の読み替えは立たない
#
#  検定CD3 実長 × 残渣長 は一歩の種類によらず一定か
#      OK なら：その一定値が「実物の大きさ」にあたり、遠近の読みが一意に立つ
#      NG なら：一定値が複数ある。遠近の読みは類ごとに分かれる
#
#  検定CD4 残渣での折れ角
#      OK なら：偶数の鎖は残渣で 180°（一本の線の上を往復する）＝真横から見たジグザグ。
#               奇数の鎖は残渣でも折れる
#      NG なら：偶数は残渣でも折れずに前進する。潰れたジグザグという読みは立たない
"""
from fractions import Fraction as F
import math, cmath
from qphi import Qp, zmul, zconj, zsigma, zre, PHI, ONE, ZERO

# ---- Z[zeta10] の道具（b13_chain_units.py と同じ構成）----
def zt(k):
    k %= 10
    base = [(1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1),(-1,1,-1,1)]
    if k < 5: return tuple(F(x) for x in base[k])
    return tuple(-F(x) for x in base[k-5])

def zadd(a,b): return tuple(a[i]+b[i] for i in range(4))
def zsub(a,b): return tuple(a[i]-b[i] for i in range(4))
def zrot(a,k): return zmul(a, zt(k))
def zabs2(a): return zre(zmul(a, zconj(a)))          # |v|^2 を Q(phi) で
def res2(a):  return zabs2(zsigma(a))                # 残渣^2
ZE = (F(0),)*4

PHI2 = (F(-1),F(1),F(-1),F(1))   # zeta^4 = phi^2 * ... 確認は下で
# phi^2 を Z[zeta10] の元として：zeta+zeta^-1 = 2cos36 = phi
PHI_Z = zadd(zt(1), zt(9))
PHI2_Z = zmul(PHI_Z, PHI_Z)

CONT = [zrot(zmul(PHI2_Z, zsub(zt(4), zt(0))), k) for k in range(10)]
SKIP = [zrot(zmul(PHI2_Z, zsub(zt(2), zt(0))), k) for k in range(10)]

C_A, C_B = CONT[5], CONT[7]        # 奇数の鎖：342 度と 54 度
S_AX, C_AX = SKIP[7], CONT[6]      # 偶数の鎖：ともに 18 度

def num(a): 
    z = cmath.exp(1j*math.pi/5)
    return sum(float(a[k])*z**k for k in range(4))
def ang(a): return math.degrees(cmath.phase(num(a))) % 360

print("一歩の確認（観察者側の計算）")
for nm, v in (("連 C_A", C_A), ("連 C_B", C_B), ("飛 S_AX", S_AX), ("連 C_AX", C_AX)):
    print(f"  {nm:9s} v={tuple(int(x) for x in v)}  |v|^2={zabs2(v)}={zabs2(v).val():.6f}"
          f"  実長={math.sqrt(zabs2(v).val()):.6f}  向き={ang(v):7.3f}deg"
          f"  残渣長={math.sqrt(res2(v).val()):.6f}  残渣向き={ang(zsigma(v)):7.3f}deg")

# ---------------- 検定CD1 ----------------
odd_pair  = zadd(C_A, C_B)
even_pair = zadd(S_AX, C_AX)
cd1 = (odd_pair == even_pair)
print(f"\n検定CD1 二歩ぶんの和")
print(f"  奇数 C_A+C_B   = {tuple(int(x) for x in odd_pair)}  実長={math.sqrt(zabs2(odd_pair).val()):.6f} 向き={ang(odd_pair):.3f}deg")
print(f"  偶数 S_AX+C_AX = {tuple(int(x) for x in even_pair)}  実長={math.sqrt(zabs2(even_pair).val()):.6f} 向き={ang(even_pair):.3f}deg")
print(f"  残渣長 {math.sqrt(res2(odd_pair).val()):.6f} / {math.sqrt(res2(even_pair).val()):.6f}")
print(f"検定CD1: {'OK' if cd1 else 'NG'}  格子ベクトルとして同一か = {cd1}")

# ---------------- 検定CD2 ----------------
# 中間の環（1歩目の着地点）を比べる
print(f"\n検定CD2 中間の環（1歩目）")
ax = C_AX   # 軸方向（18度）
def axis_comp(v):   # 軸方向成分（正規化前。両者で同じ分母なので比較に足りる）
    return zre(zmul(v, zconj(ax)))
def perp2(v):       # 軸からの外れの二乗（|v|^2 - 軸成分^2/|ax|^2）
    a = axis_comp(v); n = zabs2(ax)
    return zabs2(v) - (a*a)/n

for nm, v in (("奇数の1歩 C_A", C_A), ("偶数の1歩 S_AX", S_AX)):
    print(f"  {nm:15s} 実長={math.sqrt(zabs2(v).val()):.6f}"
          f"  軸からの外れ={math.sqrt(abs(perp2(v).val())):.6f}"
          f"  残渣の原点からの距離={math.sqrt(res2(v).val()):.6f}")
real_ratio = math.sqrt(zabs2(S_AX).val()) / math.sqrt(zabs2(C_A).val())
res_ratio  = math.sqrt(res2(S_AX).val())  / math.sqrt(res2(C_A).val())
# 厳密判定：実で偶数が短い かつ 残渣で偶数が長い
cd2 = (zabs2(S_AX) < zabs2(C_A)) and (res2(C_A) < res2(S_AX))
print(f"  実長の比（偶数/奇数）  = {real_ratio:.6f}")
print(f"  残渣長の比（偶数/奇数）= {res_ratio:.6f}")
print(f"  積 = {real_ratio*res_ratio:.6f}")
print(f"検定CD2: {'OK' if cd2 else 'NG'}  実で近い側が残渣で遠い = {cd2}")

# ---------------- 検定CD3 ----------------
print(f"\n検定CD3 実長 x 残渣長")
targets = [("飛 2つ飛ばし", S_AX), ("連 連続", C_AX),
           ("クラスタ間", even_pair)]
# 五芒星から出る一歩（phi^3）: 五芒星中心から円環中心まで
star_step = zmul(PHI_Z, zmul(PHI_Z, PHI_Z))   # phi^3（向きは問わない、長さだけ）
targets.append(("五芒星→円環 phi^3", star_step))
targets.append(("外側の殻 phi^5", zmul(star_step, zmul(PHI_Z,PHI_Z))))
prods = []
for nm, v in targets:
    pr = zmul(zabs2(v).__class__ and v, v)  # placeholder
    p2 = zabs2(v) * res2(v)                 # (実長 x 残渣長)^2 を Q(phi) で
    prods.append(p2)
    print(f"  {nm:20s} 実長={math.sqrt(zabs2(v).val()):10.6f}"
          f"  残渣長={math.sqrt(res2(v).val()):10.6f}"
          f"  積={math.sqrt(p2.val()):.6f}   積^2={p2}")
distinct = []
for p in prods:
    if not any(p == d for d in distinct): distinct.append(p)
cd3 = (len(distinct) == 1)
print(f"  相異なる積は {len(distinct)} 種: {distinct}")
print(f"検定CD3: {'OK' if cd3 else 'NG'}  一定か = {cd3}（NG なら遠近の読みは類ごとに分かれる）")

# ---------------- 検定CD4 ----------------
def bend(u, v):
    """u で来て v へ行くときの折れ角（度）。厳密には cos を Q(phi) で持てるが表示は度で。"""
    return (180 - (ang(v) - ang(u))) % 360

print(f"\n検定CD4 折れ角")
print(f"  実空間  奇数 C_A→C_B = {bend(C_A,C_B):.4f}deg   C_B→C_A = {bend(C_B,C_A):.4f}deg")
print(f"  実空間  偶数 S_AX→C_AX = {bend(S_AX,C_AX):.4f}deg   C_AX→S_AX = {bend(C_AX,S_AX):.4f}deg")
oa, ob = zsigma(C_A), zsigma(C_B)
sa, ca = zsigma(S_AX), zsigma(C_AX)
print(f"  残渣    奇数 = {bend(oa,ob):.4f}deg / {bend(ob,oa):.4f}deg")
print(f"  残渣    偶数 = {bend(sa,ca):.4f}deg / {bend(ca,sa):.4f}deg")
# 厳密判定：偶数の2歩が残渣で平行（実係数倍）かどうか。sa と ca が Q(phi) 上で同一直線か
def parallel(u, v):
    """u, v が残渣平面で同一直線上か。外積 Im(u * conj(v)) = 0 で判定（厳密）。"""
    cr = zmul(u, zconj(v))
    # Im(x) = 0 <=> x = conj(x)
    return cr == zconj(cr)
cd4 = parallel(sa, ca) and not parallel(oa, ob)
print(f"  残渣で偶数の2歩は平行か = {parallel(sa,ca)} ／ 奇数の2歩は平行か = {parallel(oa,ob)}")
print(f"検定CD4: {'OK' if cd4 else 'NG'}")

# ---------------- 鎖を実際に敷いて残渣の軌跡を出す ----------------
print(f"\n鎖の軌跡（n=9 と n=10）")
def chain(n):
    cs=[ZE]; v=[C_A,C_B] if n%2 else [S_AX,C_AX]
    for i in range(n-1): cs.append(zadd(cs[-1], v[i%2]))
    return cs
for n in (9,10):
    cs = chain(n)
    print(f"  n={n}（{'奇数ジグザグ' if n%2 else '偶数直線'}）")
    for i,c in enumerate(cs):
        a = axis_comp(c); nn = zabs2(ax)
        along = a.val()/math.sqrt(nn.val())
        off = math.sqrt(abs(perp2(c).val()))
        sg = zsigma(c)
        print(f"    環{i}: 実 軸方向={along:9.4f} 外れ={off:7.4f}"
              f"   残渣 |.|={math.sqrt(res2(c).val()):8.5f} 向き={ang(sg) if any(sg) else 0:8.3f}deg")
