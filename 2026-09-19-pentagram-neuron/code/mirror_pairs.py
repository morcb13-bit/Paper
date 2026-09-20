#  五芒星ニューロン模型  検定MR1〜MR3  左右の対・軸の上の残渣・三つめの焦点
#
#      発端（監督）：整数化した図は2を除き必ず左右対称になる。2はアンバランスなので
#                    残渣を持った鏡像。──検定IF2 で星2と星3、星4と星5 が同じ (p,q)
#                    を持っていたことに対応する。
#
#  事前登録（走らせる前に書いた）
#
#  検定MR1  二つの焦点からの (p,q) が一致する星の組を全部拾う
#      OK なら：一致する組はすべて、焦点を結ぶ直線を軸とする鏡像である
#      NG なら：一致しても鏡像でない組がある＝(p,q) の一致は鏡の言い換えではない
#
#  検定MR2  対を持たない星を拾い、軸からの外れ（残渣）を出す
#      測定。ゼロなら軸の上、ゼロでないならその値が残渣
#
#  検定MR3  軸から外した三つめの焦点で、左右の対が破れるか
#      OK なら：MR1 の組が (p,q) の三つ組では別になる
#      NG なら：三つめを足しても区別できない
#
#  負の対照  鏡を軸でない線に取ると、対が成立しない
#      NG なら：MR1 は鏡を測っていない
#
#  判定は golden() が返す整数対 (p,q) の一致で行う。鏡像かどうかの照合は座標で別に取る。
#
#  前提：wind_core.py を exec し、pentagram_neuron.py の carrier() と golden() を使う。

import math, itertools

PHI = (1 + 5 ** 0.5) / 2
TOL = 1e-6


def reflect(P, A, B):
    """直線 AB についての鏡像。"""
    ax, ay = A; bx, by = B; px, py = P
    dx, dy = bx - ax, by - ay
    L2 = dx * dx + dy * dy
    t = ((px - ax) * dx + (py - ay) * dy) / L2
    qx, qy = ax + t * dx, ay + t * dy
    return (2 * qx - px, 2 * qy - py)


def off_axis(P, A, B):
    """直線 AB からの垂線の長さ（符号つき）。"""
    ax, ay = A; bx, by = B; px, py = P
    dx, dy = bx - ax, by - ay
    return ((px - ax) * dy - (py - ay) * dx) / math.hypot(dx, dy)


def tag(S, foci):
    """焦点それぞれからの r² を黄金整数 (p,q) に直した組。None があれば None。"""
    out = []
    for Fp in foci:
        g = golden(math.dist(S, Fp) ** 2)
        if g is None:
            return None
        out.append(g)
    return tuple(out)


def run_mirror():
    F, faces, SC = carrier()
    NG = 0

    # 焦点は距離 φ⁴ の隣り合う二つ（検定IF1 と同じ取り方）
    AB = None
    for i, j in itertools.combinations(range(30), 2):
        if abs(math.dist(SC[i], SC[j]) - PHI ** 4) < 1e-6:
            AB = (i, j); break
    A, B = SC[AB[0]], SC[AB[1]]
    print(f"焦点　星{AB[0]} と 星{AB[1]}　距離 {math.dist(A,B):.6f} = φ⁴")

    # ── 検定MR1 ───────────────────────────────────────────
    print("\n検定MR1  (p,q) が一致する組は鏡像か")
    groups = {}
    for i, S in enumerate(SC):
        t = tag(S, (A, B))
        groups.setdefault(t, []).append(i)
    pairs = [g for g in groups.values() if len(g) > 1]
    singles = [g[0] for g in groups.values() if len(g) == 1]
    print(f"  一致する組 {len(pairs)} 組／対を持たない星 {len(singles)} 個"
          f"／組の大きさ {sorted(len(g) for g in pairs)}")

    bad = []
    for g in pairs:
        for i, j in itertools.combinations(g, 2):
            R = reflect(SC[i], A, B)
            if math.dist(R, SC[j]) > TOL:
                bad.append((i, j, math.dist(R, SC[j])))
    ok = not bad
    print(f"  鏡像でない組 {len(bad)} 件  {'OK' if ok else 'NG'}")
    for i, j, d in bad[:5]:
        print(f"    星{i} と 星{j}　鏡像からの外れ {d:.6f}")
    NG += 0 if ok else 1

    for g in pairs:
        i, j = g[0], g[1]
        (pa, qa), (pb, qb) = groups_key(groups, g)
        print(f"    星{i} ／ 星{j}　rA² = {pa}+{qa}φ　rB² = {pb}+{qb}φ")

    # ── 検定MR2 ───────────────────────────────────────────
    print("\n検定MR2  対を持たない星の、軸からの外れ（残渣）")
    res = sorted((abs(off_axis(SC[i], A, B)), i) for i in singles)
    on_axis = [i for d, i in res if d < TOL]
    print(f"  対を持たない星 {len(singles)} 個／うち軸の上（外れ<{TOL}）{len(on_axis)} 個")
    for d, i in res:
        g = golden(d ** 2)
        gs = f"{g[0]}+{g[1]}φ" if g else "黄金整数でない"
        print(f"    星{i:2d}　外れ {d:9.6f}　外れ² = {gs}")

    # ── 検定MR3 ───────────────────────────────────────────
    print("\n検定MR3  軸から外した三つめの焦点で対が破れるか")
    C = None
    for i in singles:
        if abs(off_axis(SC[i], A, B)) > TOL:
            C = i; break
    if C is None:
        C = max(range(30), key=lambda i: abs(off_axis(SC[i], A, B)))
    print(f"  三つめの焦点　星{C}　軸からの外れ {abs(off_axis(SC[C],A,B)):.6f}")

    allg = all(golden(math.dist(S, SC[C]) ** 2) for S in SC)
    print(f"  三つめからの r² が30個すべて黄金整数  {'はい' if allg else 'いいえ'}")

    broke = 0
    for g in pairs:
        for i, j in itertools.combinations(g, 2):
            ti = tag(SC[i], (A, B, SC[C]))
            tj = tag(SC[j], (A, B, SC[C]))
            if ti != tj:
                broke += 1
    total = sum(len(g) * (len(g) - 1) // 2 for g in pairs)
    ok = broke == total and total > 0
    print(f"  破れた組 {broke}/{total}  {'OK' if ok else 'NG'}")
    NG += 0 if ok else 1

    # ── 負の対照 ─────────────────────────────────────────
    print("\n負の対照  鏡を軸でない線に取る")
    print("  線の取り方　　　　　　　　鏡像が星に重なる個数（30中）")
    hit_axis = sum(1 for S in SC
                   if any(math.dist(reflect(S, A, B), T) < TOL for T in SC))
    print(f"    軸 AB　　　　　　　　　　　{hit_axis}")
    th_axis = math.degrees(math.atan2(B[1] - A[1], B[0] - A[0])) % 180
    print(f"  （軸 AB の方位は {th_axis:.6f}°。同じ方位の線は軸そのものなので対照から除く）")
    worst = 0
    for k in range(1, 20):
        deg = 180 * k / 20
        if min(abs(deg - th_axis), 180 - abs(deg - th_axis)) < 1e-9:
            continue
        th = math.radians(deg)
        B2 = (A[0] + math.cos(th), A[1] + math.sin(th))
        h = sum(1 for S in SC
                if any(math.dist(reflect(S, A, B2), T) < TOL for T in SC))
        worst = max(worst, h)
        print(f"    星{AB[0]} を通る {deg:5.1f}°の線　　　{h}")
    ok = worst < hit_axis
    print(f"  軸以外の最大 {worst} ＜ 軸 {hit_axis}  {'OK' if ok else 'NG'}")
    NG += 0 if ok else 1

    # ── 検定MR1b ─────────────────────────────────────────
    print("\n検定MR1b  担体の対称軸に乗っていない焦点対で同じことが起きるか")
    print("  焦点対（φ⁴）　方位　　一致する組　鏡像でない組　対を持たない星")
    rows = 0
    for i, j in itertools.combinations(range(30), 2):
        if abs(math.dist(SC[i], SC[j]) - PHI ** 4) > 1e-6:
            continue
        P, Q = SC[i], SC[j]
        deg = math.degrees(math.atan2(Q[1] - P[1], Q[0] - P[0])) % 180
        gr = {}
        skip = False
        for k, S in enumerate(SC):
            t = tag(S, (P, Q))
            if t is None:
                skip = True; break
            gr.setdefault(t, []).append(k)
        if skip:
            print(f"    星{i:2d}／星{j:2d}　{deg:6.2f}°　r² が黄金整数にならない星がある")
            rows += 1
            continue
        pr = [g for g in gr.values() if len(g) > 1]
        nb = sum(1 for g in pr for a, b in itertools.combinations(g, 2)
                 if math.dist(reflect(SC[a], P, Q), SC[b]) > TOL)
        sg = sum(1 for g in gr.values() if len(g) == 1)
        print(f"    星{i:2d}／星{j:2d}　{deg:6.2f}°　　{len(pr):2d} 組　　　{nb} 件　　　{sg} 個")
        rows += 1
        if rows >= 12:
            print("    （以下略）")
            break

    print(f"\nNG {NG} / 3")


def groups_key(groups, g):
    for k, v in groups.items():
        if v is g:
            return k
    raise KeyError


if __name__ == "__main__":
    exec(open("wind_core.py").read())
    src = open("pentagram_neuron.py").read().split("def run_tests")[0]
    exec(src)
    run_mirror()
