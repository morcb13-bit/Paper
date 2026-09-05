#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
2n+3 の歩数の側で篩が閉じることの検査。標準ライブラリのみ。

測る前に決めた検査項目（何が OK で何が NG か）:

  検定 SR1  n <= N で「2n+3 が合成数」と「n = 2ab+3a+3b+3 と書ける」が完全一致するか。
            不一致 1 件で NG。一致率では判定しない。
  検定 SR2  （落ちるのが正しい検査）係数を一つ壊した n = 2ab+3a+3b+2 で同じ判定をする。
            NG を返すのが正常。OK が返ったら検査が反証不能。
  検定 SR3  行 a が落とす数の集合が、2a+3 の奇数倍の集合と一致するか。
            a = 0..A の全行で一致しなければ NG。
  検定 SR4  行 a の対角 b=a が (2a+3)^2 か。1 行でも外れれば NG。
  検定 SR5  行 a が素数を篩う行であることと、a 自身が生き残る歩数であることが同値か。
            n <= N の全域で不一致 0 を要求。
  検定 SR6  歩数の側の 13 倍の一手 n -> 13n + 18 が、実際に 13 倍になっているか。
            -1 から 52 回まわして 13^52 に一致しなければ NG。
  検定 SR7  （落ちるのが正しい検査）一手を n -> 13n + 17 に変えたとき 13^52 になるか。
            NG を返すのが正常。OK が返ったら検査が壊れている。
  検定 SR8  記事前半の桁（真ん中の環の過不足）が n mod 3 と対応するか。
            +1 型 <-> n = 3k+2、-1 型 <-> n = 3k+1、桁 0 <-> n = 3k。不一致 0 を要求。
"""

N = 200000
A = 200

# --- 素数表（判定の外側の道具。歩数の側の主張には使わない） ---
lim = 2 * N + 3
comp = bytearray(lim + 1)
comp[0] = comp[1] = 1
i = 2
while i * i <= lim:
    if not comp[i]:
        for j in range(i * i, lim + 1, i):
            comp[j] = 1
    i += 1
def isprime(x):
    return x >= 2 and not comp[x]

def mark(N, c):
    """歩数の側だけで塗る。c=3 が本番、c=2 が負の対照。"""
    hit = bytearray(N + 1)
    a = 0
    while 3 * a + c <= N:
        b = 0
        while True:
            n = 2 * a * b + 3 * a + 3 * b + c
            if n > N:
                break
            hit[n] = 1
            b += 1
        a += 1
    return hit

ok = []

# SR1
hit = mark(N, 3)
bad = [n for n in range(N + 1) if isprime(2 * n + 3) == bool(hit[n])]
ok.append(("SR1", len(bad) == 0, "不一致 %d 件" % len(bad)))

# SR2（負の対照。落ちることを確認する）
hit2 = mark(N, 2)
bad2 = [n for n in range(N + 1) if isprime(2 * n + 3) == bool(hit2[n])]
ok.append(("SR2", len(bad2) == 0, "壊した形での不一致 %d 件。NG が正常" % len(bad2)))

# SR3
allrow = True
for a in range(A + 1):
    row, b = [], 0
    while True:
        n = 2 * a * b + 3 * a + 3 * b + 3
        if n > N:
            break
        row.append(2 * n + 3)
        b += 1
    q = 2 * a + 3
    if row != [q * (2 * b + 3) for b in range(len(row))]:
        allrow = False
        break
ok.append(("SR3", allrow, "行 0..%d" % A))

# SR4
diag = all(2 * (2 * a * a + 6 * a + 3) + 3 == (2 * a + 3) ** 2 for a in range(A + 1))
ok.append(("SR4", diag, "対角 b=a"))

# SR5
live = [not hit[n] for n in range(N + 1)]          # 生き残る歩数
prow = [isprime(2 * a + 3) for a in range(N + 1)]  # 素数を篩う行
ok.append(("SR5", live == prow, "行番号と生き残りの同値"))

# SR6
n = -1
for _ in range(52):
    n = 13 * n + 18
ok.append(("SR6", 2 * n + 3 == 13 ** 52, "13^52"))

# SR7（落ちるのが正しい）
m = -1
for _ in range(52):
    m = 13 * m + 17
ok.append(("SR7", 2 * m + 3 == 13 ** 52, "壊した一手。NG が正常"))

# SR8
d = {0: 0, 1: -1, 2: +1}
bad8 = 0
for n in range(0, N + 1):
    p = 2 * n + 3
    if not isprime(p) or p == 3:
        continue
    digit = +1 if p % 3 == 1 else -1
    if d[n % 3] != digit:
        bad8 += 1
ok.append(("SR8", bad8 == 0, "桁と n mod 3、不一致 %d 件" % bad8))

for name, good, note in ok:
    print("検定 %s : %s  (%s)" % (name, "OK" if good else "NG", note))
