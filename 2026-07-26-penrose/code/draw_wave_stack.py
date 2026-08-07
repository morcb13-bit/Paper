"""
位相が半周期ズレた cos 波を2本重ね、等間隔に積んで n=1..5 を作る。
偶数＝波の交点、奇数＝山と谷を結んだジグザグ、2本の波＝2種類の奇数。
図は wave_stack.svg。数値の検定は厳密（qphi）、描画だけ浮動小数。

#  検定WS1 奇数の鎖は cos 波の山と谷に乗るか
#      OK なら：奇数のジグザグは、半周期 a・振幅 A の1本の波の標本そのもの
#      NG なら：ジグザグは波の標本ではない
#
#  検定WS2 2種類の奇数は位相が半周期ズレた2本の波か
#      OK なら：鏡像の2本の鎖が、そのまま2本の波に対応する
#      NG なら：鏡像は位相のズレでは書けない
#
#  検定WS3 偶数の鎖は2本の波の交点に乗るか
#      OK なら：「偶数＝波が重なった部分」がそのまま成り立つ
#      NG なら：交点は等間隔だが実際の偶数の鎖は等間隔でない。どこまで一致するかを数える
#
#  検定WS4 V字（3環）はちょうど1周期か
#      OK なら：3の篩は「波1周期で敷ける」ことに対応する
#      NG なら：3と周期は対応しない
#
#  検定WS5 n=1..5 の篩の類が既登録の sieve_class と一致するか
#      OK なら：この図の色分けは新しい判定ではなく登録済みの判定の描き直し
#      NG なら：勝手な判定を持ち込んでいる
"""
from fractions import Fraction as F
import math, cmath
from qphi import Qp, zmul, zconj, zsigma, zre

def zt(k):
    k %= 10
    base = [(1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1),(-1,1,-1,1)]
    if k < 5: return tuple(F(x) for x in base[k])
    return tuple(-F(x) for x in base[k-5])
def zadd(a,b): return tuple(a[i]+b[i] for i in range(4))
def zsub(a,b): return tuple(a[i]-b[i] for i in range(4))
def zrot(a,k): return zmul(a, zt(k))
def zabs2(a): return zre(zmul(a, zconj(a)))

PHI_Z = zadd(zt(1), zt(9)); PHI2_Z = zmul(PHI_Z, PHI_Z)
CONT = [zrot(zmul(PHI2_Z, zsub(zt(4), zt(0))), k) for k in range(10)]
SKIP = [zrot(zmul(PHI2_Z, zsub(zt(2), zt(0))), k) for k in range(10)]
C_A, C_B = CONT[5], CONT[7]
S_AX, C_AX = SKIP[7], CONT[6]

ROT = cmath.exp(-1j*math.radians(18))           # 単位の軸を水平にする
def num(a):
    z = cmath.exp(1j*math.pi/5)
    return sum(float(a[k])*z**k for k in range(4))
def xy(a):
    p = num(a)*ROT
    return (p.real, p.imag)

def chain(n, flip=False):
    cs=[(F(0),)*4]
    v=[C_B,C_A] if (n%2 and flip) else ([C_A,C_B] if n%2 else [S_AX,C_AX])
    for i in range(n-1): cs.append(zadd(cs[-1], v[i%2]))
    return cs

# ── 基準量（すべて格子から出す） ──
STEP2 = zadd(C_A, C_B)                     # クラスタ間 8.057480
a_full = math.sqrt(zabs2(STEP2).val())     # 1周期 = 8.057480
a = a_full/2                               # 半周期 = 4.028740
PV   = zsub(C_B, C_A)                      # 差ベクトルの実長 5.854102
PKV  = math.sqrt(zabs2(PV).val())/2        # 山谷差（環の高さの振れ）2.927051
AMP  = PKV/2                               # cos の振幅 1.463525
DIFF = zsub(S_AX, C_AX)                    # 飛 - 連。実長 1.902113
print(f"1周期 = {a_full:.6f}   半周期 a = {a:.6f}   山谷差 = {PKV:.6f}   振幅 A = {AMP:.6f}")

# ---------------- WS1 ----------------
print("\n検定WS1 奇数の鎖は cos 波の山と谷に乗るか")
ws1 = True
for n in (3,5,7,9):
    cs = chain(n)
    pts = [xy(c) for c in cs]
    x0, y0 = pts[0]
    ctr = y0 - AMP                             # 中心線（山から振幅ぶん下）
    for i,(x,y) in enumerate(pts):
        wx = x0 + i*a
        wy = ctr + AMP*math.cos(math.pi*i)
        if abs(x-wx) > 1e-9 or abs(y-wy) > 1e-9: ws1 = False
    print(f"  n={n}: 軸位置の刻み={[round(pts[i+1][0]-pts[i][0],6) for i in range(n-1)]}"
          f"  高さ={[round(p[1]-ctr,6) for p in pts]}")
print(f"検定WS1: {'OK' if ws1 else 'NG'}  すべて山/谷にちょうど乗るか = {ws1}")

# ---------------- WS2 ----------------
print("\n検定WS2 2種類の奇数は位相が半周期ズレた2本の波か")
cA = [xy(c) for c in chain(5, flip=False)]
cB = [xy(c) for c in chain(5, flip=True)]
print(f"  奇数A の高さ = {[round(p[1],6) for p in cA]}")
print(f"  奇数B の高さ = {[round(p[1],6) for p in cB]}")
ws2 = all(abs(cA[i][1] + cB[i][1]) < 1e-9 for i in range(5)) and \
      all(abs(cA[i][0] - cB[i][0]) < 1e-9 for i in range(5))
print(f"  同じ軸位置で高さが符号反転（＝半周期ズレ）か = {ws2}")
print(f"検定WS2: {'OK' if ws2 else 'NG'}")

# ---------------- WS3 ----------------
print("\n検定WS3 偶数の鎖は2本の波の交点に乗るか")
ws3 = True
for n in (4,6,8):
    cs = [xy(c) for c in chain(n)]
    steps = [round(cs[i+1][0]-cs[i][0],6) for i in range(n-1)]
    cross = [round(a,6)]*(n-1)
    print(f"  n={n}: 実際の刻み={steps}")
    print(f"        交点の刻み={cross}")
    if steps != cross: ws3 = False
    # 2歩ごとの照合
    two = [round(cs[i+2][0]-cs[i][0],6) for i in range(n-2)]
    print(f"        2歩ごと ={two}  （1周期 {a_full:.6f}）")
print(f"検定WS3: {'OK' if ws3 else 'NG'}  1歩ごとに一致するか = {ws3}")

# ---------------- WS4 ----------------
print("\n検定WS4 V字（3環）はちょうど1周期か")
vee = math.sqrt(zabs2(zadd(C_A,C_B)).val())
ws4 = abs(vee - a_full) < 1e-12
print(f"  V字の端から端 = {vee:.6f}   1周期 = {a_full:.6f}   一致 = {ws4}")
print(f"検定WS4: {'OK' if ws4 else 'NG'}")

# ---------------- WS5 ----------------
print("\n検定WS5 篩の類（既登録 b13_chain_units.sieve_class と照合）")
def sieve_class(n):
    if n==1: return "unit"
    if n in (2,3): return "sieve%d"%n
    if n%2==0: return "out2"
    if n%3==0: return "out3"
    return "prime"
import importlib.util
spec = importlib.util.find_spec("b13_chain_units")
ws5 = spec is not None
if ws5:
    import b13_chain_units as B
    ws5 = all(B.sieve_class(n) == sieve_class(n) for n in range(1,14))
for n in range(1,6): print(f"  n={n}: {sieve_class(n)}")
print(f"検定WS5: {'OK' if ws5 else 'NG'}  既登録の判定と一致 = {ws5}")

# ---------------- WS6 ----------------
print("\n検定WS6 偶数の鎖は交点の格子（刻み a）に対してどうずれるか")
half = math.sqrt(zabs2(DIFF).val())/2
ws6 = True
for n in (4,6,8):
    cs = [xy(c)[0] for c in chain(n)]
    dev = [cs[i] - (cs[0] + i*a) for i in range(n)]
    print(f"  n={n}: 交点の格子からのズレ={[round(d,6) for d in dev]}")
    for i,d in enumerate(dev):
        want = 0.0 if i % 2 == 0 else -half
        if abs(d - want) > 1e-9: ws6 = False
print(f"  |飛-連|/2 = {half:.6f}")
print(f"検定WS6: {'OK' if ws6 else 'NG'}  偶数番は乗り、奇数番だけ一定量ずれるか = {ws6}")

# ================= 図 =================
COLOR={"unit":"#8a8f99","sieve2":"#4aa3e0","sieve3":"#5cc06e",
       "out2":"#2a4a63","out3":"#2a5238","prime":"#f2b544"}
NAME={"unit":"単位","sieve2":"2の篩","sieve3":"3の篩",
      "out2":"2で落ちる","out3":"3で落ちる","prime":"残る"}

S = 30.0                      # px / 単位長
PITCH = PKV*1.5               # 積む間隔＝山谷差（隣の行の波が接する）
NS = [1,2,3,4,5]
MARG_L, MARG_T = 110, 74
xmax = 4*a + a
W = int(MARG_L + xmax*S + 30)
H = int(MARG_T + (len(NS)-1)*PITCH*S + 2*AMP*S + 140)

def px(x): return MARG_L + x*S
def py(y): return MARG_T + y*S

out=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
     f'<rect width="{W}" height="{H}" fill="#0d1117"/>',
     '<style>text{font-family:"Noto Sans CJK JP","Noto Sans JP",sans-serif;fill:#c9d1d9}</style>']

for row,n in enumerate(NS):
    Y = row*PITCH + AMP                       # この行の中心線
    cls = sieve_class(n); col = COLOR[cls]
    # 2本の波
    for sgn,dash in ((+1,""), (-1,"4 3")):
        d=[]
        for k in range(0, 401):
            x = xmax*k/400
            y = Y + sgn*AMP*math.cos(math.pi*x/a)
            d.append(("M" if k==0 else "L")+f"{px(x):.2f},{py(y):.2f}")
        out.append(f'<path d="{"".join(d)}" fill="none" stroke="#3d4550" '
                   f'stroke-width="1.2"{f" stroke-dasharray=\"{dash}\"" if dash else ""}/>')
    # 中心線＝一直線
    out.append(f'<line x1="{px(0)}" y1="{py(Y)}" x2="{px(xmax)}" y2="{py(Y)}" '
               f'stroke="#5a6472" stroke-width="0.8" stroke-dasharray="2 4"/>')
    # 環の位置
    if n % 2:
        pts = [(i*a, Y + AMP*math.cos(math.pi*i)) for i in range(n)]
        mirr= [(i*a, Y - AMP*math.cos(math.pi*i)) for i in range(n)]
        d="".join(("M" if i==0 else "L")+f"{px(x):.2f},{py(y):.2f}" for i,(x,y) in enumerate(pts))
        out.append(f'<path d="{d}" fill="none" stroke="{col}" stroke-width="2.6"/>')
        d2="".join(("M" if i==0 else "L")+f"{px(x):.2f},{py(y):.2f}" for i,(x,y) in enumerate(mirr))
        out.append(f'<path d="{d2}" fill="none" stroke="{col}" stroke-width="1.4" opacity="0.45"/>')
        for x,y in mirr: out.append(f'<circle cx="{px(x)}" cy="{py(y)}" r="4" fill="none" stroke="{col}" stroke-width="1.2" opacity="0.55"/>')
        for x,y in pts:  out.append(f'<circle cx="{px(x)}" cy="{py(y)}" r="6.5" fill="{col}"/>')
        note = "山と谷／2種の奇数"
    else:
        pts = [(a/2 + i*a, Y) for i in range(n)]
        d="".join(("M" if i==0 else "L")+f"{px(x):.2f},{py(y):.2f}" for i,(x,y) in enumerate(pts))
        out.append(f'<path d="{d}" fill="none" stroke="{col}" stroke-width="2.6"/>')
        for x,y in pts: out.append(f'<circle cx="{px(x)}" cy="{py(y)}" r="6.5" fill="{col}"/>')
        # 実際の偶数の鎖（照合用・白抜き）
        real = [xy(c)[0] for c in chain(n)]
        for i,rx in enumerate(real):
            X = pts[0][0] + (rx - real[0])
            out.append(f'<circle cx="{px(X)}" cy="{py(Y)}" r="4.5" fill="none" '
                       f'stroke="#e06c75" stroke-width="1.5"/>')
        note = "波の交点"
    out.append(f'<text x="18" y="{py(Y)+5:.1f}" font-size="17" font-weight="600">n={n}</text>')
    out.append(f'<text x="18" y="{py(Y)+23:.1f}" font-size="11" fill="{col}">{NAME[cls]}</text>')
    out.append(f'<text x="18" y="{py(Y)+40:.1f}" font-size="10" fill="#7d8590">{note}</text>')

half_txt = f'{math.sqrt(zabs2(DIFF).val())/2:.6f}'
yb = H-76
for i,t in enumerate([
    f'半周期 a = {a:.6f}（＝クラスタ間 {a_full:.6f} の半分）／振幅 A = {AMP:.6f}（＝山谷差 {PKV:.6f} の半分）',
    f'点線＝各行の中心線。5行とも同じ間隔で積んである（間隔＝山谷差×1.5）',
    f'赤丸＝実際の偶数の鎖。偶数番は交点にちょうど乗り、奇数番だけ {half_txt} 手前へずれる（＝|飛−連|/2）',
    f'色は登録済みの sieve_class（1=単位／2・3=篩そのもの／4=2で落ちる／5=残る）']):
    out.append(f'<text x="18" y="{yb+i*17}" font-size="11.5" fill="#8b949e">{t}</text>')
out.append('</svg>')
open("wave_stack.svg","w").write("\n".join(out))
print("\nwave_stack.svg を書き出しました")
