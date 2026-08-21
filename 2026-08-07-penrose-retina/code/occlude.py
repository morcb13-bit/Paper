#  検定M4  隠れを跨ぐ ── 条件を緩めるのではなく、一歩を外挿する
#
#      発端（監督）：赤ん坊は、動くアヒルの軌道を一部隠すとその先を追えない。
#
#      前のターンで「隠れを扱うには『部品が全部乗る』を『m個中k個乗る』に緩める
#      ことになり、単調に狭まる性質が失われる」と書いた。だがそれは隠れの扱い方を
#      一つしか考えていなかった。**緩めずに、一歩を外挿して跨ぐ**という手がある。
#      条件は全称のまま、跨いだ先の一枚で当て直す。
#      赤ん坊が跨げないのは、跨ぐ前に一歩が決まっていないからだと読める。
#
#      対象：三角形が一定の一歩で動く列。途中の g 枚で対象が隠れる（セルが消える）。
#            隠れる前の n 枚で一歩を決め、隠れの先の一枚で当て直す
#      測る量：(a) 隠れる前の n 枚で一歩が正しく決まるか（山の1位が実際の一歩か）
#              (b) 跨いだあとに残った支持セルの適合率・再現率
#              (c) 一枚ずつ続くことを要求する追い方（前の検定の規則）が跨げるか
#
#      OK なら：一歩が決まっていれば g が伸びても跨げる。決まっていなければ跨げない
#      NG なら：一歩が決まっていても跨げない／決まっていなくても跨げてしまう
#      必ず落ちる設定：外挿に間違った一歩を使うと、跨いだ先で支持セルが消えること

exec(open('mg_lib.py').read())

import random

def seq(n_before, gap, v, nz=80, seed=40):
    """隠れのある列を作る。隠れている枚は対象のセルを消す"""
    cA, cB = spot(-7, 6), spot(8, -7)
    Es = []; truth = None
    total = n_before + gap + 1
    for t in range(total):
        a = put(TRI, U.zadd(cA, tuple(t * x for x in v)))
        if t == 0:
            truth = a
        hidden = n_before <= t < n_before + gap
        b = put(SQR, cB)
        E = (set() if hidden else a) | b | noise(nz, seed + 7 * t)
        Es.append(E)
    return Es, truth

def step_from(Es, n, moving=True):
    """最初の n 枚だけで一歩を決める。moving なら「止まっている」(一歩ゼロ)は外す
       ── 動いているものを探しているのだから、止まっている群は別の答えである"""
    cand = [(len(track(Es[:n], d)), i) for i, d in enumerate(STEPS)
            if not (moving and d == ZERO)]
    sc = sorted(cand, reverse=True)
    return STEPS[sc[0][1]], sc[0][0], sc[1][0]

def cross(Es, n, gap, d):
    """最初の n 枚で追い、隠れを跨いで最後の一枚で当て直す（条件は全称のまま）"""
    keep = track(Es[:n], d)
    t = n + gap
    out = set()
    for q in keep:
        r = fast_land(tuple(F(a) + b for a, b in zip(q, tuple(t * x for x in d))))
        if r is not None and r in Es[t]:
            out.add(q)
    return out

V = min(STEPS[1:], key=lambda v: abs(math.hypot(*U.xy(v)) - 3.1))
WRONG = max(STEPS[1:], key=lambda v: math.hypot(*U.xy(v)) if math.hypot(*U.xy(v)) < 5.2 else 0)
print("追う一歩 %.3f（画面座標の長さ）／雑音80セル\n" % math.hypot(*U.xy(V)))

print("(a) 隠れる前の枚数で一歩が決まるか（止まっている群＝一歩ゼロは別の答えとして外す）")
print("  枚数   山の1位   2位   1位は実際の一歩か   止まっている群の大きさ")
for n in (2, 3, 4, 5, 6, 8):
    Es, _ = seq(n, 0, V)
    d, top, snd = step_from(Es, n)
    still = len(track(Es[:n], ZERO))
    print("  %3d   %6d %6d   %-8s %14d" % (n, top, snd, "はい" if d == V else "いいえ", still))

print("\n(b) 隠れを跨ぐ（一歩を外挿・条件は全称のまま）")
print("  隠れ枚数   前2枚          前4枚          前6枚")
for gap in (0, 1, 2, 3, 5):
    cells = []
    for n in (2, 4, 6):
        Es, truth = seq(n, gap, V)
        d, _, _ = step_from(Es, n)
        g = cross(Es, n, gap, d)
        p, r = pr(g, truth)
        cells.append("%2d個 %.2f/%.2f" % (len(g), p, r))
    print("  %6d     %s" % (gap, "   ".join(cells)))
print("  各欄は 支持セル数 適合率/再現率")

print("\n(c) 一枚ずつ続くことを要求する追い方（前の検定の規則）")
for gap in (0, 1, 2):
    Es, truth = seq(6, gap, V)
    g = track(Es, V)
    print("  隠れ %d 枚：支持セル %d 個" % (gap, len(g)))

print("\n必ず落ちる設定：間違った一歩で外挿する")
for gap in (0, 2):
    Es, truth = seq(6, gap, V)
    g = cross(Es, 6, gap, WRONG)
    p, r = pr(g, truth)
    print("  隠れ %d 枚：支持セル %d 個・適合率 %.2f" % (gap, len(g), p))
