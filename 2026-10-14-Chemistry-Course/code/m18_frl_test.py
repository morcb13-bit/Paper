# m18: 検定 FR-L（事前登録：2026-10-14 12:07 JST、監督の承認後に実行）
#
# 規則（固定）：L3² = L² − n·u（n＝出会う腕の本数 1・2・3、L² と u は全結合に共通の正の整数）
#   L3 ＝ 炭素の中心どうしの3次元の距離（平衡核間距離 re）。4次元の中心間ベクトルの影
#   予言：L1² − L2² ＝ L2² − L3² ＝ u（記録の長さ²が等間隔）
# 記録（m16 と同じ、10^-5 Å 単位の整数と幅）：
#   C–C 1.522 ±0.0005（Harmony 1990）、C=C 1.3305 ±0.0010（Craig 2006）、C≡C 1.20289 ±0.00001（Cazzoli 2008）
# 判定：差の区間 [lo1−hi2, hi1−lo2] と [lo2−hi3, hi2−lo3] に共通の正の整数 u があれば OK、なければ NG
#   かつ L3² ＝ L² − 3u ＞ 0（L² ＝ L1² ＋ u で決まる）
# 正の対照：規則どおりに作った等間隔の長さ² → OK が返ること
# 負の対照：Lₑ² ∝ n²（1:4:9）の規則 → 差が 3u と 5u。同じ判定器で結果を返すこと
# 記録は m16 で見ているので盲検ではない

REC = [("C–C", 152200, 50), ("C=C", 133050, 100), ("C≡C", 120289, 1)]

def sq_iv(v, d): return ((v - d) ** 2, (v + d) ** 2)

def judge(ivs, k1, k2, label):
    """差1 = k1·u、差2 = k2·u を満たす正の整数 u があるか"""
    (a_lo, a_hi), (b_lo, b_hi), (c_lo, c_hi) = ivs
    d1 = (a_lo - b_hi, a_hi - b_lo)   # L1² − L2² の区間
    d2 = (b_lo - c_hi, b_hi - c_lo)   # L2² − L3² の区間
    # u の区間：d1/k1 と d2/k2 を整数で（切り上げ・切り捨て）
    u1 = (-(-d1[0] // k1), d1[1] // k1)
    u2 = (-(-d2[0] // k2), d2[1] // k2)
    lo, hi = max(u1[0], u2[0], 1), min(u1[1], u2[1])
    ok = lo <= hi
    if ok:  # L3² > 0 の確認
        u = lo; L2 = a_lo + k1 * u * 0 + u * 0  # 表示用
        ok = c_lo > 0
    print(f"  {label}")
    print(f"    差1 の区間 {d1}  → u の区間 {u1}")
    print(f"    差2 の区間 {d2}  → u の区間 {u2}")
    print(f"    共通の u：{'[%d, %d]' % (lo, hi) if lo <= hi else 'なし'}  → {'OK' if ok else 'NG'}")
    return ok, d1, d2

ivs = [sq_iv(v, d) for _, v, d in REC]
print("== 記録の長さ²（10^-10 Å² 単位の整数区間）")
for (n, v, d), iv in zip(REC, ivs): print(f"  {n}: {iv}")

print("== 正の対照（規則どおりの等間隔：L²=24000000000, u=3000000000, 幅±1）")
L2, u = 24000000000, 3000000000
syn = [(L2 - k * u - 1, L2 - k * u + 1) for k in (1, 2, 3)]
judge(syn, 1, 1, "正の対照 → OK が返るはず")

print("== 負の対照（Lₑ² ∝ n²：差 3u と 5u）を記録に")
judge(ivs, 3, 5, "n² 規則")

print("== FR-L（Lₑ² ∝ n：差 u と u）を記録に")
ok, d1, d2 = judge(ivs, 1, 1, "FR-L")
print(f"  差1 の中央 {(d1[0]+d1[1])//2}、差2 の中央 {(d2[0]+d2[1])//2}（10^-10 Å²）")
print("  FR-L:", "OK" if ok else "NG")
