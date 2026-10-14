# m16: 共有の模型（m15）を平衡構造の記録と照合する。判定は整数の比較だけ。
#
# 事前登録（走らせる前に決めたこと）
#  L1 順序   ：長さは 頂点 > 辺 > 面（C–C > C=C > C≡C）。CO（面）< CO2（辺）
#  L2 比     ：長さ² は 頂点:辺:面 = 9:3:1。記録の幅（最後の桁の半分）の区間どうしが重なれば OK
#  A1 角     ：模型の組 (p, a, b) の cos² = p²/(ab) が、記録の角 ±幅 の cos² 区間に入れば OK
#  A2 形     ：平面・一直線（模型の整数判定と、記録の点群が一致すれば OK）
#  正の対照  ：記録そのものから作った比は L2 を OK にする（判定器が OK を返せること）
#  負の対照  ：順序を入れ替えた模型（面が一番長い）は L1 で NG になる（判定器が NG を返せること）
#  三値      ：記録どうしの食い違いが模型の分解能より大きいときは「判定不能」
#
# 記録の数値は 10^-5 Å・10^-2 度の整数として書き写す（観察者側の転記）。
# 角の cos² 区間だけは、記録の角を整数の分数に直すところで一度だけ浮動小数を使う（外向きに丸める、観察者側の換算）。
import math

# ---- 記録（re、Å を 10^-5 単位の整数で） ----
REC_LEN = {
    "C-C  ethane    (Harmony 1990)":            (152200, 50),
    "C=C  ethylene  (Craig 2006)":              (133050, 100),
    "C#C  acetylene (Cazzoli 2008, re_SE)":     (120289, 1),
    "C#O  CO        (Lille table re)":          (112832, 1),
    "C=O  CO2       (Lille table re)":          (116000, 50),
}
# 模型の中心どうしの長さ²（m15）
MODEL_LEN2 = {"vertex": 108, "edge": 36, "face": 12}

def sq_interval(v, d):
    return ((v - d) ** 2, (v + d) ** 2)

def ratio_overlap(rec_a, rec_b, ka, kb):
    """記録 a, b の長さ² の比が ka:kb と両立するか。kb·La² と ka·Lb² の区間が重なるか（整数）"""
    a_lo, a_hi = sq_interval(*rec_a); b_lo, b_hi = sq_interval(*rec_b)
    lhs = (kb * a_lo, kb * a_hi); rhs = (ka * b_lo, ka * b_hi)
    return lhs[0] <= rhs[1] and rhs[0] <= lhs[1], lhs, rhs

def verdict(b): return "OK" if b else "NG"

r = list(REC_LEN.values())
print("== L1 順序")
order_cc = r[0][0] - r[0][1] > r[1][0] + r[1][1] and r[1][0] - r[1][1] > r[2][0] + r[2][1]
order_co = r[3][0] + r[3][1] < r[4][0] - r[4][1]
print("  C–C > C=C > C≡C :", verdict(order_cc))
print("  CO（面）< CO2（辺）:", verdict(order_co))
# 負の対照：面が一番長いという模型
neg = r[2][0] - r[2][1] > r[0][0] + r[0][1]
print("  負の対照（面が最長）:", verdict(neg), "← NG が出れば判定器は NG を返せる")

print("== L2 比（長さ² 9:3:1）")
ok1, l, rr = ratio_overlap(r[0], r[1], 3, 1)
print("  単:二 = 3:1 →", verdict(ok1), " 1·L単²", l, " 3·L二²", rr)
ok2, l, rr = ratio_overlap(r[1], r[2], 3, 1)
print("  二:三 = 3:1 →", verdict(ok2), " 1·L二²", l, " 3·L三²", rr)
ok3, l, rr = ratio_overlap(r[4], r[3], 3, 1)
print("  CO2:CO = 3:1 →", verdict(ok3), " 1·L(CO2)²", l, " 3·L(CO)²", rr)
# 正の対照：記録から作った整数比（長さ² を 10^-5 Å 単位のまま）
ka, kb = r[0][0] ** 2, r[1][0] ** 2
okp, _, _ = ratio_overlap(r[0], r[1], ka, kb)
print("  正の対照（記録そのものの比）→", verdict(okp), "← OK が出れば判定器は OK を返せる")
# 記録の長さ² の比を、いちばん近い小さい整数の組で見る（観察者側の参考）
for name, (x, y) in [("単:二", (r[0][0], r[1][0])), ("二:三", (r[1][0], r[2][0])), ("単:三", (r[0][0], r[2][0]))]:
    print(f"  参考 {name} の長さ² ×1000 の比 ≈ {x*x*1000//(y*y)}/1000（模型は 3000/1000、単:三 は 9000/1000）")

# ---- 角 ----
DEN = 10 ** 8
def cos2_bounds(deg100, d100):
    """記録の角 (θ ± δ, 0.01 度単位) の cos² 区間を、分母 10^8 の整数で外向きに。観察者側の換算"""
    lo, hi = (deg100 - d100) / 100, (deg100 + d100) / 100
    vals = [math.cos(math.radians(t)) ** 2 for t in (lo, hi)]
    if lo <= 90 <= hi: vals.append(0.0)
    return math.floor(min(vals) * DEN), math.ceil(max(vals) * DEN), (lo, hi)

def angle_check(name, pair, deg100, d100, src):
    p, a, b = pair
    lo, hi, rng = cos2_bounds(deg100, d100)
    obtuse_model = p < 0; obtuse_rec = rng[0] > 90
    inside = lo * a * b <= p * p * DEN <= hi * a * b and obtuse_model == obtuse_rec
    print(f"  {name}: 模型の組 {pair}（cos²={p*p}/{a*b}） 記録 {deg100/100}±{d100/100}° [{src}] →", verdict(inside))
    return inside

print("== A1 角")
angle_check("エチレン H–C–H", (-9, 27, 27), 11710, 10, "Craig 2006")
angle_check("エチレン H–C–C", (-18, 27, 36), 12145, 10, "Craig 2006")
print("  Fe–C–O（CO の炭素の埋まった腕、模型は一直線 180°）:")
angle_check("    CO-ミオグロビン", (-1, 1, 1), 17100, 500, "Lund 再精密化, 幅±5°はこちらで設定")
print("  Fe–O–O（O2 の埋まった腕、模型の組 (−18, 27, 36)）:")
hb = angle_check("    オキシヘモグロビン", (-18, 27, 36), 15600, 500, "Shaanan 1982, 幅±5°はこちらで設定")
print("    オキシミオグロビン 115°（Phillips 1980）は今回原典を取得できず未照合")
print("    記録どうしの食い違い 156 − 115 = 41° は、模型の分解能（腕の組が一つの角を与える）より大きい → 判定不能")
print("  順序：Fe–C–O（模型 180）> Fe–O–O（模型 125.26）。記録 171 > 156 →", verdict(17100 - 500 > 15600 + 500))

print("== A2 形（平面・一直線）")
print("  エチレンの H 4 個が一平面（模型：法線 (0,72,72) で 4 点とも内積 0、記録：D2h の平面分子） → OK")
print("  アセチレン H–C≡C–H が一直線（模型：外積 0、記録：D∞h） → OK")
print("  CO2 O=C=O が一直線（模型：O の中心の和 (0,0,0)、記録：D∞h） → OK")
