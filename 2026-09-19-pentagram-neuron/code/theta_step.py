#  五芒星ニューロン模型  検定NE5a〜NE5c  θ の周期と一歩の大きさ
#
#      発端（引継書 v256 §4-1）：一歩あたり θ を 72°/60 = 1.2° 進める割り当ては
#                    別の記事からの借用で、五芒星に伏せた状態で成り立つかは未検算。
#
#  θ の定義（借用元の記事を持っていないので、装置の中で立て直したもの）
#      面12枚それぞれを、自分の面の中心を通る法線まわりに θ だけ回す。
#      頂点は3枚の面に属するので θ≠0 では像が3つに分かれ 20×3 = 60 点。
#      θ=0 では重なって20点。これを「分裂頂点60」と読む。
#      ※ 借用元の定義がこれと違えば、別のものを測ったことになる。
#
#  事前登録（走らせる前に書いた）
#
#  検定NE5a  縮退する θ
#      OK なら：像が20点に落ちるのは θ が 72° の倍数のときだけ
#      NG なら：他の θ でも落ちる＝周期が72°でない
#
#  検定NE5b  一歩の大きさ
#      OK なら：δ = 72/60 = 1.2°/歩 で、60歩ちょうどで縮退に戻る
#      負の対照：δ = 1.0°/歩 では60歩後に戻らない。戻るならこの検査は NG を返せない
#
#  検定NE5c  伏せた面と担体
#      OK なら：伏せた面の5頂点が担体の五芒星に嵌まる θ と、NE5a で縮退する θ が一致
#      NG なら：ずれる＝書き込みと読み出しの交替が構造から出ていない
#
#  前提：wind_core.py を exec し、pentagram_neuron.py の carrier()・dodeca() を使う。

import math, itertools

PHI = (1 + 5 ** 0.5) / 2
TOL = 1e-7


def face_frame(V, f):
    """面 f の中心と法線。"""
    P = [V[i] for i in f]
    c = tuple(sum(p[k] for p in P) / 5 for k in range(3))
    n = tuple(c)                      # 正十二面体は中心が原点なので面の中心が法線向き
    L = math.dist((0, 0, 0), n)
    return c, tuple(x / L for x in n)


def rot(p, c, n, th):
    """点 p を、c を通る軸 n のまわりに th 回す（ロドリゲス）。"""
    v = [p[k] - c[k] for k in range(3)]
    ct, st = math.cos(th), math.sin(th)
    d = sum(v[k] * n[k] for k in range(3))
    cr = (n[1] * v[2] - n[2] * v[1], n[2] * v[0] - n[0] * v[2], n[0] * v[1] - n[1] * v[0])
    return tuple(c[k] + v[k] * ct + cr[k] * st + n[k] * d * (1 - ct) for k in range(3))


def split_points(V, FA, deg):
    """面ごとに deg 度回したときの、頂点の像60点。重複を潰した個数も返す。"""
    th = math.radians(deg)
    pts = []
    for f in FA:
        c, n = face_frame(V, f)
        for i in f:
            pts.append(rot(V[i], c, n, th))
    uniq = []
    for p in pts:
        if not any(math.dist(p, q) < 1e-9 for q in uniq):
            uniq.append(p)
    return pts, len(uniq)


def run_theta():
    NG = 0
    V, E, FA, FD, D = dodeca()

    # ── 検定NE5a ─────────────────────────────────────────
    print("検定NE5a  縮退する θ（0.05°刻みで 0〜144° を掃く）")
    hits = []
    deg = 0.0
    while deg <= 144.0 + 1e-9:
        _, m = split_points(V, FA, deg)
        if m == 20:
            hits.append(round(deg, 2))
        deg += 0.05
    print(f"  20点に落ちた θ：{hits}")
    ok = hits == [0.0, 72.0, 144.0]
    print(f"  72°の倍数だけか  {'OK' if ok else 'NG'}")
    NG += 0 if ok else 1

    print("  θ ごとの点の数（抜き出し）")
    for d in (0, 1.2, 12, 36, 60, 71.9, 72, 108, 144):
        _, m = split_points(V, FA, d)
        print(f"    θ = {d:6.2f}°　点 {m:2d}")

    # ── 検定NE5b ─────────────────────────────────────────
    print("\n検定NE5b  一歩の大きさ（60歩で戻るか）")
    print("  一歩      60歩後の θ    点の数    縮退")
    for delta in (1.2, 1.0, 0.6, 2.4, 1.25):
        th = delta * 60
        _, m = split_points(V, FA, th % 360)
        print(f"    {delta:4.2f}°　　{th:7.2f}°　　　{m:2d}　　  {'はい' if m == 20 else 'いいえ'}")
    _, m12 = split_points(V, FA, 1.2 * 60)
    _, m10 = split_points(V, FA, 1.0 * 60)
    ok = (m12 == 20) and (m10 != 20)
    print(f"  1.2°で戻り 1.0°で戻らない  {'OK' if ok else 'NG'}")
    NG += 0 if ok else 1

    print("  途中の歩で縮退しないか（1.2°／歩・1〜59歩）")
    mid = [k for k in range(1, 60) if split_points(V, FA, 1.2 * k)[1] == 20]
    print(f"    縮退した歩：{mid if mid else 'なし'}  {'OK' if not mid else 'NG'}")
    NG += 0 if not mid else 1

    # ── 検定NE5c ─────────────────────────────────────────
    print("\n検定NE5c  伏せた面と担体の五芒星")
    F, faces, SC = carrier()
    st = None
    for a, c in faces:
        if abs(a - 2.9389) < 0.01:
            st = [tuple(float(t) for t in U.xy(p)) for p in c]
            break
    cx = sum(p[0] for p in st) / 10
    cy = sum(p[1] for p in st) / 10
    rad = sorted(set(round(math.hypot(p[0] - cx, p[1] - cy), 6) for p in st))
    tips = [p for p in st if abs(math.hypot(p[0] - cx, p[1] - cy) - rad[-1]) < 1e-6]
    print(f"  五芒星の頂点10個／中心からの距離 {rad}（尖り {len(tips)} 個）")

    tipang = sorted(math.degrees(math.atan2(p[1] - cy, p[0] - cx)) % 360 for p in tips)
    base = tipang[0]
    print("  伏せた面の5頂点を θ だけ回したときの、尖りとの最大ずれ")
    print("    θ        最大ずれ    嵌まる")
    fits = []
    deg = 0.0
    while deg <= 144.0 + 1e-9:
        ang = sorted((base + deg + 72 * k) % 360 for k in range(5))
        err = max(min(abs(a - b), 360 - abs(a - b)) for a, b in zip(ang, tipang))
        if err < 1e-6:
            fits.append(round(deg, 2))
        deg += 0.05
    for d in (0, 1.2, 36, 71.9, 72, 108, 144):
        ang = sorted((base + d + 72 * k) % 360 for k in range(5))
        err = max(min(abs(a - b), 360 - abs(a - b)) for a, b in zip(ang, tipang))
        print(f"    {d:6.2f}°　{err:9.4f}°　 {'はい' if err < 1e-6 else 'いいえ'}")
    print(f"  嵌まる θ：{fits}")
    ok = fits == hits
    print(f"  縮退する θ と一致するか  {'OK' if ok else 'NG'}")
    NG += 0 if ok else 1

    print(f"\nNG {NG} / 4")


if __name__ == "__main__":
    exec(open("wind_core.py").read())
    exec(open("pentagram_neuron.py").read().split("def run_tests")[0])
    run_theta()
