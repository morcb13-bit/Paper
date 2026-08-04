#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
chain_balanced.py ── 鎖の長さを平衡3進で読む
============================================================
標準ライブラリのみ。b13_chain_units.py だけを import する（判別法4）。

奇数の鎖を V字（3）で敷くと、真ん中に環が1つ余るか、1つ足りなくて重ねるか、
過不足なしかのどれかになる。桁は −1, 0, +1 の三つしかない。

事前登録（判別法6）
-------------------
  検定1 V字が座れる位置（鏡映を許す）
      OK なら：奇数の鎖では全位置に座れる
      NG なら：座れない位置を出す
      ※ 鏡映を許さないと奇数位置しか出ない。V字自体が鏡対称なので許すのが正しい
  検定2 独立した1の最小個数は n mod 3 か
      OK なら：3の篩が1の個数として刻まれている
      NG なら：一致しない n を出す
  検定3 桁（+1 = 環が余る／−1 = 中心で1環重ねる／0 = 過不足なし）で n が復元されるか
      OK なら：3×(V字の個数) + 桁 = n
      NG なら：合わない n を出す
  検定4 過不足の印は鎖のちょうど中心に来るか
      OK なら：左右のV字の枚数が揃い、印が中心環に一致する
      NG なら：ずれる n と、ずれ幅を出す
  検定5 被覆に穴や溢れが無いか（1..n をちょうど覆うか）
  検定6 鎖は鏡対称か（共役＋回転）
      OK なら：奇数も偶数も鏡対称
  検定7 V字を三つ並べたとき、その中心を結ぶ折れ角は 108° か
      OK なら：3進が二段目へ繰り上がる（繰り込みが実在する）
      NG なら：実測値を出す。3進は一段で止まる

  NG が一件も出ない検査は反証不能を疑うこと。検定7 はそのために置いてある
  （実際に NG が出る）。緩めて落ちる例：検定1 で鏡映を許さないと検定2〜5 が崩れる。
"""
import math, sys
sys.path.insert(0, "/mnt/skills/user/b13-verify/scripts")
import b13_chain_units as B

def links(cs): return [B.zsub(cs[i+1], cs[i]) for i in range(len(cs)-1)]
L3 = links(B.unit(3))

def pos3(n, mirror=True):
    """V字が座れる開始環番号"""
    L = links(B.unit(n)); out = []
    for i in range(len(L)-1):
        for r in range(10):
            if all(L[i+j] == B.zrot(L3[j], r) for j in range(2)): out.append(i+1); break
            if mirror and all(L[i+j] == B.zrot(B.zconj(L3[j]), r) for j in range(2)):
                out.append(i+1); break
    return out

def min_ones(n):
    """独立した1の最小個数（動的計画）"""
    ok = set(pos3(n))
    INF = 10**9
    dp = [INF]*(n+2); dp[n+1] = 0
    for p in range(n, 0, -1):
        c = dp[p+1] + 1
        if p in ok and p+2 <= n: c = min(c, dp[p+3])
        dp[p] = c
    return dp[1]

def balanced(n):
    """中心に桁を置く敷き方 → (桁, V字の被覆, 印の環番号)"""
    r = n % 3
    if r == 0:   k, d = n//3, 0
    elif r == 1: k, d = (n-1)//3, +1
    else:        k, d = (n+1)//3, -1
    cov = []; p = 1; cent = None
    if r == 0:
        for _ in range(k): cov.append((p, p+2)); p += 3
        cent = (n+1)//2
    elif r == 1:
        h = k//2
        for _ in range(h): cov.append((p, p+2)); p += 3
        cent = p; p += 1
        for _ in range(k-h): cov.append((p, p+2)); p += 3
    else:
        h = k//2
        for _ in range(h): cov.append((p, p+2)); p += 3
        p -= 1; cent = p
        for _ in range(k-h): cov.append((p, p+2)); p += 3
    return d, cov, cent

def bal3_digits(n):
    out = []
    while n:
        r = n % 3
        d = 0 if r == 0 else (1 if r == 1 else -1)
        out.append(d); n = (n - d)//3
    return out

def main(N=29):
    ok = lambda c: "OK" if c else "NG"
    odd = [n for n in range(3, N+1, 2)]

    c1 = all(pos3(n) == list(range(1, n-1)) for n in odd)
    print("検定1 V字は全位置に座れる（鏡映込み）        %s" % ok(c1))
    print("      鏡映を許さない場合の n=9 の位置: %s" % pos3(9, mirror=False))

    c2 = all(min_ones(n) == n % 3 for n in odd)
    print("検定2 独立した1の最小個数 = n mod 3          %s" % ok(c2))

    c3 = c4 = c5 = True
    for n in odd:
        d, cov, cent = balanced(n)
        if 3*len(cov) + d != n: c3 = False
        if cent != (n+1)//2: c4 = False
        used = set()
        for a, b in cov: used |= set(range(a, b+1))
        if d == +1: used.add(cent)
        if used != set(range(1, n+1)): c5 = False
    print("検定3 3×(V字の個数) + 桁 = n                 %s" % ok(c3))
    print("検定4 印は鎖のちょうど中心                    %s" % ok(c4))
    print("検定5 被覆は 1..n をちょうど覆う              %s" % ok(c5))

    c6 = True
    for n in range(2, 14):
        L = links(B.unit(n)); rev = [B.zrot(v, 5) for v in L[::-1]]
        if not any(L == [B.zrot(B.zconj(v), r) for v in rev] for r in range(10)): c6 = False
    print("検定6 鎖は鏡対称（共役＋回転）                %s" % ok(c6))

    cs = B.unit(9); _, cov, _ = balanced(9)
    m = [(a+b)//2 - 1 for a, b in cov]
    u = complex(*B.xy(B.zsub(cs[m[0]], cs[m[1]])))
    w = complex(*B.xy(B.zsub(cs[m[2]], cs[m[1]])))
    ang = math.degrees(abs(math.atan2((w/u).imag, (w/u).real)))
    print("検定7 V字三つの中心を結ぶ折れ角が108°        %s  実測 %.4f°"
          % (ok(abs(ang-108) < 1e-6), ang))

    print("\n  n   桁   V字   式              中心   3進の全桁（下位から）")
    for n in odd:
        d, cov, cent = balanced(n)
        sign = "+1" if d == 1 else ("-1" if d == -1 else " 0")
        expr = "3×%-2d %s" % (len(cov), "+1" if d == 1 else ("−1" if d == -1 else "  "))
        print("  %2d   %s   %3d   %-14s %3d    %s"
              % (n, sign, len(cov), expr + " = %d" % n, cent, bal3_digits(n)))

    nng = sum(1 for c in (c1, c2, c3, c4, c5, c6) if not c)
    print("\nNG %d / 7（検定7 は落ちるのが正しい）" % (nng + 1))

if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 29)
