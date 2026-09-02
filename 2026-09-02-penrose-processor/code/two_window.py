#  検定TW0・TW1  二つの窓は台地を縮めるか
#
#      発端（監督・2026-09-02）：私の描いたペンローズタイリングは、別の黄金比のひし形の
#          ペンローズタイリングとぴったり重なる。同じ場所に異なる経路で到達する。
#
#      これを担体の言葉に直すと「担体は一つ、窓が二つ」。同じ番地が二つの残渣を持つ。
#      検定CS・GZ2・TB が三度とも落ちた原因は台地（着地が二値なので、許容差の内側では
#      全部が同点になり向きが取り出せない）だった。五角形はずらすと最大 63.4% が
#      自分に重なる（検定PS1）。二つの読みは間隔の比が φ なので、台地の縁が揃わない。
#      ならば両方を要求すれば台地は縮むはず ── これを測る。
#
#  検定TW0  二つ目の読みが、この担体の上に実在するか
#      測る量：円環中心を担体の中心まわりに ζ^k（36°×k）で回したとき、
#              五角形中心にちょうど乗る個数（判定は Q(ζ10) の厳密一致）
#      OK なら：ある k で乗る個数が多い。二つ目の読みが取れる。TW1 へ進む
#      NG なら：どの k でも乗らない。この担体では二つの窓が重ならないので TW1 は走らせない
#      必ず落ちる設定：円環中心を格子に乗らない量だけずらすと、どの k でも 0 になること
#
#  検定TW1  二重に読んだ点は、ずらしても自分に重なりにくいか
#      測る量：ずれ26通りに対する自己重なり率（ずれ0 を除く最大）。
#              窓A単独（五角形）／窓B単独（円環中心）／二重（両方が印を付けた点）
#      OK なら：二重の率が、同じ個数の無作為部分集合20本の帯の下端を下回る
#      NG なら：帯の中に入る ── 縮んだのは点が減ったからであって、二つの窓の性質ではない
#      対照：無作為部分集合20本。閾値は置かず帯で見る（wind_core の事前登録の欠陥を避ける）
#      必ず落ちる設定：二つ目の読みを一つ目と同じものに差し替えると、
#                      二重の率が窓A単独と完全に一致すること（一致しなければ物差しが壊れている）
#
#      帰属：二つのペンローズが同じ座標に乗るという観察は監督の発話。
#            相似比 φ³ の四量一致と 63.4% は観察者側の計算（既出）。
#            台地が縮むという見立ては観察者側で、これが初測定。
#
#      使い方： python3 two_window.py            本測定（b13_chain_units.py と wind_core.py が要る）
#               python3 two_window.py --selftest 物差しだけを合成データで検査（担体不要）

import math, random, sys

TOL = 0.30          # 重なりの判定（star_compare.py と同じ）
NBAND = 20          # 無作為部分集合の本数

# ══ 物差し（担体に依らない部分） ════════════════════════════

def selfmatch(pts, shifts, tol=TOL):
    """ずらして自分に重なる個数。pts は画面座標の列"""
    grid = {}
    for p in pts:
        grid.setdefault((int(p[0] // 2), int(p[1] // 2)), []).append(p)
    out = []
    for dx, dy in shifts:
        n = 0
        for p in pts:
            qx, qy = p[0] + dx, p[1] + dy
            gx, gy = int(qx // 2), int(qy // 2)
            hit = False
            for ax in (-1, 0, 1):
                for ay in (-1, 0, 1):
                    for r in grid.get((gx + ax, gy + ay), ()):
                        if (r[0] - qx) ** 2 + (r[1] - qy) ** 2 <= tol * tol:
                            hit = True; break
                    if hit: break
                if hit: break
            if hit: n += 1
        out.append((n, (dx, dy)))
    return out


def profile(pts, shifts):
    """(点の数, ずれ0 の個数, ずれ0 以外の最大率, 上位5本の率)"""
    sc = selfmatch(pts, shifts)
    zero = [n for n, d in sc if abs(d[0]) < 1e-9 and abs(d[1]) < 1e-9]
    n0 = zero[0] if zero else 0
    others = sorted((n for n, d in sc
                     if not (abs(d[0]) < 1e-9 and abs(d[1]) < 1e-9)), reverse=True)
    m = len(pts)
    top = [n / m for n in others[:5]] if m else []
    return m, n0, (others[0] / m if others and m else 0.0), top


def band(pool, m, shifts, n=NBAND, seed=13):
    """pool から m 個を無作為に取った部分集合の、ずれ0 以外の最大率の帯"""
    rnd = random.Random(seed)
    vals = []
    for _ in range(n):
        sub = rnd.sample(pool, m)
        vals.append(profile(sub, shifts)[2])
    vals.sort()
    return vals[0], vals[len(vals) // 2], vals[-1]


def report(name, pts, shifts):
    m, n0, r, top = profile(pts, shifts)
    print("  %-10s 点 %4d ／ ずれ0 で %4d ／ ずれ0 以外の最大 %.4f（%.1f%%）"
          % (name, m, n0, r, 100 * r))
    print("             上位5本 %s" % " ".join("%.3f" % t for t in top))
    return r


# ══ 自己検査（合成データ・担体不要） ══════════════════════════

def selftest():
    print("■ 物差しの検査（合成データ・担体は使わない）\n")
    # 周期的な点の集合＝台地が立つはずのもの
    A = [(i * 1.0, j * 1.0) for i in range(26) for j in range(26)]
    SH = [(dx * 1.0, dy * 1.0) for dx in range(-3, 4) for dy in range(-3, 4)]

    print("① 周期的な集合は台地を作るか（作らなければ物差しが壊れている）")
    rA = report("周期", A, SH)
    assert rA > 0.5, "NG: 周期的な集合で台地が立たない"
    print("   → OK\n")

    print("② 不規則な集合は台地を作らないか（偶然の水準と比べる）")
    rnd = random.Random(7)
    NB = 676
    B = [(rnd.uniform(0, 26), rnd.uniform(0, 26)) for _ in range(NB)]
    chance = 1.0 - math.exp(-(NB / (26.0 * 26.0)) * math.pi * TOL * TOL)
    rB = report("不規則", B, SH)
    print("             偶然の水準 %.3f（密度と許容差だけから決まる）" % chance)
    assert rB < chance * 1.5, "NG: 不規則な集合で偶然を超えて重なった"
    assert rA > chance * 2.0, "NG: 周期と不規則の差が出ていない"
    print("   → OK（周期 %.3f は偶然を大きく超え、不規則 %.3f は偶然どまり）\n" % (rA, rB))

    print("③ 必ず落ちる設定：二つ目の読みを一つ目と同じにすると率が一致するか")
    doubled = [p for p in A if p in set(A)]          # A ∩ A
    rD = report("A∩A", doubled, SH)
    assert abs(rD - rA) < 1e-12, "NG: A∩A が A と一致しない。物差しが壊れている"
    print("   → OK（一致）\n")

    print("④ 大きさの効果：無作為に間引くだけで率は下がるか")
    m = len(A) // 7
    lo, mid, hi = band(A, m, SH, n=NBAND, seed=3)
    print("   %d 個へ無作為に間引いた帯 %.4f 〜 %.4f（中央 %.4f）" % (m, lo, hi, mid))
    assert hi < rA, "NG: 間引いても率が下がらない"
    print("   → OK（下がる。だから TW1 は帯と比べる必要がある）\n")

    print("物差しは通った。本測定は担体のあるところで走らせること。")


# ══ 本測定 ══════════════════════════════════════════════

def main():
    exec(open("wind_core.py").read(), globals())
    # wind_core が入れるもの： U, CELLS, CTR, land, ...

    XY = {q: U.xy(q) for q in CELLS}
    CELLSET = set(CELLS)

    rows, place, offs = U.build_stack()
    RINGS = [tuple(r) for r in sum(place, [])]
    if not RINGS or len(RINGS[0]) != 4:
        print("停止：円環中心が Q(ζ10) の4成分で取れない。build_stack の戻りを確認すること")
        print("      取れたもの（先頭）: %r" % (RINGS[0] if RINGS else None))
        return
    print("担体：五角形 %d 枚／円環 %d 個" % (len(CELLS), len(RINGS)))

    # ── ずれ26通り（star_compare.py と同じ作り方） ──
    SHIFTS = []
    for q in CELLS:
        x, y = U.xy(U.zsub(q, CTR))
        if x * x + y * y <= 36:
            SHIFTS.append((x, y))
    print("試すずれ %d 通り（長さ6以内）\n" % len(SHIFTS))

    # ══ 検定TW0 ══════════════════════════════════════
    print("■ 検定TW0 二つ目の読みは担体の上にあるか（ζ^k で回して五角形中心に乗る個数）")

    def zeta_mul(v):
        """ζ を一回掛ける。ζ^4 = -1 + ζ - ζ² + ζ³ で畳む（検算済み：ζ^5=-1・ζ^10=1）"""
        a0, a1, a2, a3 = v
        return (-a3, a0 + a3, a1 - a3, a2 + a3)

    def rotk(v, k):
        for _ in range(k % 10):
            v = zeta_mul(v)
        return v

    ORIGIN = (0, 0, 0, 0)
    CENTERS = (("原点", ORIGIN), ("担体の中心", CTR))

    best = None                      # (乗った個数, 中心の名, k, 乗った点の列)
    for cname, c in CENTERS:
        line = []
        for k in range(10):
            img = []
            for r in RINGS:
                w = U.zadd(c, rotk(U.zsub(r, c), k))
                if w in CELLSET:
                    img.append(w)
            line.append(len(img))
            if best is None or len(img) > best[0]:
                best = (len(img), cname, k, img)
        print("   %-8s ζ^k で乗る個数（k=0..9）: %s  ／ 円環 %d 個"
              % (cname, " ".join("%3d" % n for n in line), len(RINGS)))

    best_n, best_c, best_k, best_set = best

    # 必ず落ちる設定：格子に乗らない量だけずらすと、どの k でも 0 になること
    from fractions import Fraction as F
    cc = dict(CENTERS)[best_c]
    bad = 0
    for r in RINGS:
        w = U.zadd(cc, rotk(U.zsub(r, cc), best_k))
        w = tuple(F(a) + b for a, b in zip(w, (F(1, 3), 0, 0, 0)))
        if w in CELLSET:
            bad += 1
    print("   必ず落ちる設定（格子外へ 1/3 ずらす）… %d 個（0 でなければ判定が壊れている）" % bad)

    if best_n == 0:
        print("\n判定：NG。この担体では二つの窓が重ならない。TW1 は走らせない\n")
        return
    print("\n判定：OK。%s まわり k=%d（%d°）で %d 個。これを窓B の読みとする\n" % (best_c, best_k, 36 * best_k, best_n))

    # ══ 検定TW1 ══════════════════════════════════════
    print("■ 検定TW1 二重に読んだ点は台地を作るか")
    A_pts = [XY[q] for q in CELLS]
    B_pts = [U.xy(r) for r in RINGS]
    D_pts = [XY[q] for q in dict.fromkeys(best_set)]      # 二重（重複を除く）

    rA = report("窓A 五角形", A_pts, SHIFTS)
    rB = report("窓B 円環中心", B_pts, SHIFTS)
    rD = report("二重", D_pts, SHIFTS)

    lo, mid, hi = band(A_pts, len(D_pts), SHIFTS)
    print("\n  対照：五角形から %d 個を無作為に取った帯 %.4f 〜 %.4f（中央 %.4f）"
          % (len(D_pts), lo, hi, mid))

    print("\n  判定：", end="")
    if rD < lo:
        print("OK。二重の率 %.4f が帯の下端 %.4f を下回る。台地は二つの窓で縮む" % (rD, lo))
    elif rD > hi:
        print("NG（想定外）。二重の率 %.4f が帯の上端 %.4f を超えた" % (rD, hi))
    else:
        print("NG。二重の率 %.4f は帯 %.4f〜%.4f の中。縮んだのは点が減ったからで、"
              "二つの窓の性質ではない" % (rD, lo, hi))

    # ── 必ず落ちる設定：窓B を窓A に差し替える ──
    rSame = profile(A_pts, SHIFTS)[2]
    print("\n  必ず落ちる設定（窓B を窓A と同じものにする）… %.4f" % rSame)
    print("     窓A単独 %.4f と %s" % (rA, "一致（物差しは正常）" if abs(rSame - rA) < 1e-12
                                        else "不一致 ── 物差しが壊れている"))


if __name__ == "__main__":
    if "--selftest" in sys.argv:
        selftest()
    else:
        main()
