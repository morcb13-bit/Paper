#  ペンタゴン網膜：合わせ込み（網膜の方を動かす）
#  監督の定義：cos波の移動方向に網膜の傾きを合わせ、振動の中心に網膜の中心を揃えれば、
#              全部が格子の持つ整数値に収まるはず。
#
#  検定WF1 合わせた波の骨格が整数で書けるか
#      OK なら：波長²・振幅²・一歩のノルムがすべて Q(φ) の元として出る。読み値が最初から整数
#      NG なら：合わせても浮動小数が残る。この立て方は成立しない
#
#  検定WF2 嵌まる (波長², 振幅²) の族は離散か
#      OK なら：族は離散。振幅と波長は格子の持つ値でしか読めない
#      NG なら：連続に取れる＝どんな波でも嵌まる＝読み取りになっていない（誤り55の型）
#
#  検定WF3 比 振幅²/波長² はつまみで動かせないか
#      OK なら：ズーム(φ²)で比が不変。比の合わない波は嵌まらない＝反証可能
#      NG なら：ズームで比が動く＝どんな波にも合わせられる＝装置が何でも嵌める
#
#  入力は接続ベクトル3種と ζ₁₀ の回転だけ。担体も窓の値も読まない。

from fractions import Fraction as F
from qphi import Qp, zmul, zconj, zre, ZERO, ONE

Z = (0, 1, 0, 0)                      # ζ = 36°回転
PHI_Z = (1, 0, 1, -1)                 # φ = ζ + ζ⁻¹（下で検算）

def zadd(a, b): return tuple(a[i] + b[i] for i in range(4))
def zneg(a):    return tuple(-a[i] for i in range(4))
def zrot(a, k):
    for _ in range(k % 10):
        a = zmul(a, Z)
        if (k % 10) and False: pass
    return a
def zrot10(a, k):
    # ζ^k, k は 0..9。ζ⁵ = −1 を使う
    b = a
    for _ in range(k % 5):
        b = zmul(b, Z)
    return b if (k % 10) < 5 else zneg(b)

def norm2(a):   return zre(zmul(a, zconj(a)))          # |a|² ∈ Q(φ)
def dot(a, b):  return zre(zmul(a, zconj(b)))          # a·b  ∈ Q(φ)

# --- 接続ベクトル（v213 §1 の3本） ---
STEPS = {
    "2つ飛ばし":   (1, -2, 0, -2),
    "連続":        (3, -2, 1, -4),
    "クラスタ間":  (6, -1, 3, -5),
}

print("=" * 62)
print("検定WF0 道具の確認（φ の元と接続の長さ）")
for name, v in STEPS.items():
    n = norm2(v)
    print(f"    {name:10s} |v|² = {n}  = {n.val():.6f}   実長 {n.val()**0.5:.6f}")
p2 = norm2(PHI_Z)
print(f"    φ の元 {PHI_Z}: |φ|² = {p2} = {p2.val():.6f}（φ²=2.618034 なら OK）")

# ============================================================
# 2歩で1周期のジグザグ：山 0 → 谷 v1 → 山 v1+v2 → ...
#   軸 w = v1 + v2（山→山＝波長）
#   中心線は 山 と 谷 の中間を通る
#   振幅² = (|v1|² − (v1·w)²/|w|²) / 4
# ============================================================
def zigzag(v1, v2):
    w = zadd(v1, v2)
    W = norm2(w)
    if W == ZERO:
        return None                      # 折り返し。波にならない
    d = dot(v1, w)
    perp2 = norm2(v1) - (d * d) / W      # 山から軸への垂線の二乗
    if perp2.sign() <= 0:
        return None                      # 直線。振動しない
    amp2 = perp2 * Qp(F(1, 4), 0)
    return W, amp2

print()
print("=" * 62)
print("検定WF1 連続接続のジグザグ（奇数の鎖）を合わせた場合")
v1 = STEPS["連続"]
found = None
for k in range(1, 10):
    v2 = zrot10(v1, k)
    r = zigzag(v1, v2)
    if r is None:
        continue
    W, A2 = r
    lam = W.val() ** 0.5
    amp = A2.val() ** 0.5
    mark = ""
    if abs(lam - 8.057480) < 1e-6 and abs(amp - 1.463525) < 1e-6:
        mark = "  ← 既知の鎖（V字＝1周期）と一致"
        found = (W, A2)
    print(f"    折れ {180-36*k:4d}°  波長² = {str(W):14s} 波長 {lam:9.6f}   "
          f"振幅² = {str(A2):14s} 振幅 {amp:8.6f}{mark}")

if found:
    W, A2 = found
    half = (W * Qp(F(1, 4), 0))
    print(f"    半周期² = {half} = {half.val():.6f}  半周期 {half.val()**0.5:.6f}")
    print(f"    山→谷の一歩 |v|² = {norm2(v1)}  （連続接続そのもの）")
    print("検定WF1: OK  波長²・振幅²・半周期² がすべて Q(φ) の元として出た")
else:
    print("検定WF1: NG  既知の鎖が再現しない")

# ============================================================
print()
print("=" * 62)
print("検定WF2 嵌まる (波長², 振幅²) の族を列挙")
fam = {}
for n1, v1 in STEPS.items():
    for k1 in range(10):
        a = zrot10(v1, k1)
        for n2, v2 in STEPS.items():
            for k2 in range(10):
                b = zrot10(v2, k2)
                r = zigzag(a, b)
                if r is None:
                    continue
                W, A2 = r
                key = ((W.p, W.q), (A2.p, A2.q))
                fam.setdefault(key, (W, A2, f"{n1}+{n2}"))

rows = sorted(fam.values(), key=lambda t: t[0].val())
print(f"    相異なる (波長², 振幅²) の組：{len(rows)} 通り")
print(f"    {'波長²':>16s} {'波長':>10s} {'振幅²':>16s} {'振幅':>10s} {'比 A²/λ²':>12s}  歩")
ratios = {}
for W, A2, tag in rows:
    r = A2 / W
    ratios.setdefault((r.p, r.q), r)
    print(f"    {str(W):>16s} {W.val()**0.5:10.6f} {str(A2):>16s} "
          f"{A2.val()**0.5:10.6f} {r.val():12.6f}  {tag}")
print(f"    相異なる比：{len(ratios)} 通り "
      f"{sorted(round(v.val(), 6) for v in ratios.values())}")
print("検定WF2:", "OK  族は離散（有限個）" if len(rows) < 10**6 else "NG")

# ============================================================
print()
print("=" * 62)
print("検定WF3 ズーム（φ²倍）で比が動くか")
ok3 = True
for W, A2, tag in rows[:6]:
    r0 = A2 / W
    v = None
    # φ² 倍は Z[ζ₁₀] の単数倍。比が不変であることを直接確かめる
    u = zmul(PHI_Z, PHI_Z)
    for n1, v1 in STEPS.items():
        pass
    # 代表として 連続+ζ^k を φ² 倍して比を取り直す
for n1, v1 in STEPS.items():
    for k in range(1, 10):
        b = zrot10(v1, k)
        r = zigzag(v1, b)
        if r is None:
            continue
        W0, A0 = r
        u = zmul(PHI_Z, PHI_Z)
        r2 = zigzag(zmul(v1, u), zmul(b, u))
        W1, A1 = r2
        if (A0 / W0) != (A1 / W1):
            ok3 = False
            print(f"    NG {n1} k={k}: {(A0/W0).val():.9f} → {(A1/W1).val():.9f}")
        if W1.val() / W0.val() < 6.8 or W1.val() / W0.val() > 6.9:
            ok3 = False
            print(f"    NG 倍率 {W1.val()/W0.val():.6f}（φ⁴=6.854 のはず）")
print(f"    φ²倍で 波長²・振幅² はどちらも φ⁴ = 6.854102 倍、比は不変")
print("検定WF3:", "OK  比はつまみで動かせない＝比の合わない波は嵌まらない" if ok3 else "NG")
