#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
stack_rings.py ── kabai_ring.py で閉じた輪を、軸方向に積んで管にする。

積む一歩も6次元の整数ベクトル。輪と同じく、通るか通らないかは整数と
支持関数だけで決まる。

  事前登録

  検定1 積む一歩が接触しているか
      OK なら：上下の層のあいだに、面接触ちょうどの対が1組以上ある
      NG なら：層が浮いている（ただ離して置いただけ）

  検定2 重なりゼロ
      OK なら：層内・層間のすべての対で、中心差が 2·RT の内部に入らない
      NG なら：食い込んでいる

  検定3 らせんか
      OK なら：積む一歩に回転が伴い、何層かで純平行移動に戻る
      NG なら：回転を伴わない（ただの平行移動の積み重ね）

  対照   接触しない一歩・重なる一歩が実在すること（＝検定1・2が落ちうる）
"""

import itertools
import json
import math
import kabai_ring as K


def load_ring(n):
    d = json.load(open("kabai_ring.json"))["rings"]
    return [tuple(x) for x in d[str(n)]["centers6"]]


def overlap(a, b):
    """2つの RT が食い込むか（中心差が 2·RT の内部にあるか）。"""
    d = (a[0]-b[0], a[1]-b[1], a[2]-b[2])
    return K.support_ok(d, 2.0)


def main():
    N = 15
    C = load_ring(N)
    P = [K.phys(c) for c in C]
    g = K.signed_perm(K.GEN[2], 2*math.pi/5)
    axis = K.GEN[2]
    L = K.norm(axis)
    axis = K.mul3(axis, 1.0/L)

    stepv = [s[0] for s in K.STEP]
    contact = set(stepv)

    # 積む一歩の候補：面接触1本、または2本の和
    cand = list(stepv)
    for a, b in itertools.product(stepv, repeat=2):
        cand.append(tuple(a[i]+b[i] for i in range(6)))
    cand = list(dict.fromkeys(cand))

    print("== 積む一歩の探索 ==")
    print("   候補 %d 本 × 回転 5 通り" % len(cand))

    ok = []
    no_contact = overlapped = 0
    for r in range(5):
        M = [[1 if i == j else 0 for j in range(6)] for i in range(6)]
        for _ in range(r):
            M = K.matmul(g, M)
        RC = [K.matvec(M, c) for c in C]
        for t in cand:
            up = [tuple(c[i]+t[i] for i in range(6)) for c in RC]
            pu = [K.phys(c) for c in up]
            # 輪の重心がどう動いたかで判定する。回転の中心は原点であって
            # 輪の軸ではないので、t だけを見ると横ずれを見落とす
            m0 = (sum(p[0] for p in P)/len(P), sum(p[1] for p in P)/len(P),
                  sum(p[2] for p in P)/len(P))
            m1 = (sum(p[0] for p in pu)/len(pu), sum(p[1] for p in pu)/len(pu),
                  sum(p[2] for p in pu)/len(pu))
            dv = (m1[0]-m0[0], m1[1]-m0[1], m1[2]-m0[2])
            h = K.dot(dv, axis)
            lat = math.sqrt(max(0.0, K.dot(dv, dv) - h*h))
            if h <= 1e-9 or lat > 3.0:
                continue
            bad = False
            touch = 0
            for a in range(len(P)):
                for b in range(len(pu)):
                    if overlap(P[a], pu[b]):
                        bad = True
                        break
                    d = tuple(up[b][i] - C[a][i] for i in range(6))
                    if d in contact:
                        touch += 1
                if bad:
                    break
            if bad:
                overlapped += 1
                continue
            if touch == 0:
                no_contact += 1
                continue
            ok.append((r, t, h, lat, touch))
    ok_narrow = [o for o in ok if o[3] <= 1e-6]

    print("   重なって落ちた %d 本／接触せず落ちた %d 本／成立 %d 本"
          % (overlapped, no_contact, len(ok)))
    print("== 検定1・2 の対照 ==")
    print("   重なりで落ちる例 %s／接触なしで落ちる例 %s"
          % ("有り" if overlapped else "**無し**",
             "有り" if no_contact else "**無し**"))
    if not ok:
        print("   成立 0 → NG")
        return

    # 一番低く積める一歩（＝一番よく詰まる）を採る
    # らせんとして閉じるものだけ採る：S(x)=g^r x + t を5回合成した並進部分が
    # 軸にぴったり乗る（横成分ちょうど0）ことを整数のまま確かめる
    def screw_pitch(r, t):
        M = [[1 if i == j else 0 for j in range(6)] for i in range(6)]
        Mr = [[1 if i == j else 0 for j in range(6)] for i in range(6)]
        for _ in range(r):
            Mr = K.matmul(g, Mr)
        T = [0]*6
        for _ in range(5):
            Tv = K.matvec(M, t)
            T = [T[k]+Tv[k] for k in range(6)]
            M = K.matmul(Mr, M)
        p = K.phys(tuple(T))
        hh = K.dot(p, axis)
        ll = math.sqrt(max(0.0, K.dot(p, p) - hh*hh))
        return tuple(T), hh, ll

    screw = [o for o in ok if o[0] != 0 and screw_pitch(o[0], o[1])[2] < 1e-9]
    pick = screw if screw else ok
    pick.sort(key=lambda r: (r[2], -r[4]))
    r, t, h, lat, touch = pick[0]
    print("\n== 採った一歩 ==")
    print("   回転 g^%d／t = %s" % (r, t))
    print("   輪ごと軸方向に %.6f 上がる／横ずれ %.9f／接触対 %d 組"
          % (h, lat, touch))

    print("== 検定3 らせんか ==")
    period = 1 if r == 0 else 5 // math.gcd(5, r)
    Tv, pitch, plat = screw_pitch(r, t)
    print("   5層ぶんの並進 = %s" % (Tv,))
    print("   軸成分 %.9f／**横成分 %.12f**／1層 %.9f"
          % (pitch, plat, pitch/5))
    print("   4+12φ = %.9f との差 %.3e"
          % (4+12*K.PHI, abs(pitch-(4+12*K.PHI))))
    print("   → %s（%d層で純平行移動に戻るらせん軸）"
          % ("OK" if plat < 1e-9 and r else "NG", period))

    # 積む
    LAYERS = 10
    M = [[1 if i == j else 0 for j in range(6)] for i in range(6)]
    for _ in range(r):
        M = K.matmul(g, M)

    layers = []
    cur = list(C)
    for k in range(LAYERS):
        layers.append(list(cur))
        cur = [tuple(K.matvec(M, c)[i] + t[i] for i in range(6)) for c in cur]

    allc = [c for lay in layers for c in lay]
    allp = [K.phys(c) for c in allc]
    print("\n== 全層の検査 ==")
    print("   層 %d ／ RT %d 個" % (LAYERS, len(allc)))
    bad = 0
    for a in range(len(allp)):
        for b in range(a+1, len(allp)):
            if overlap(allp[a], allp[b]):
                bad += 1
    print("   食い込んでいる対 %d → %s" % (bad, "OK" if bad == 0 else "NG"))
    print("   重複した中心 %d" % (len(allc) - len(set(allc))))

    lo = min(K.dot(p, axis) for p in allp)
    hi = max(K.dot(p, axis) for p in allp)
    print("   管の高さ %.6f（1層あたり %.6f）" % (hi - lo, h))

    # 描画
    mid = K.FRAME((sum(p[0] for p in allp)/len(allp),
                   sum(p[1] for p in allp)/len(allp),
                   sum(p[2] for p in allp)/len(allp)))
    shift = K.mul3(mid, -1.0)
    for tag, az, el in (("iso", 0.55, 1.02), ("top", 0.0, 0.0),
                        ("side", 0.20, 1.5708)):
        name = "kabai_tube_%s.svg" % tag
        n, W, H = K.draw(allp, name, az, el, shift)
        print("  %-22s 多角形 %d 枚" % (name, n))

    json.dump({"ring": N, "layers": LAYERS, "rot": r, "step": list(t),
               "rise": h, "lateral": lat, "contacts_per_layer": touch,
               "screw_period": period, "pitch": pitch,
               "pitch_exact": "4+12phi", "rt_total": len(allc),
               "overlaps": bad},
              open("kabai_tube.json", "w"), indent=1)


if __name__ == "__main__":
    main()
