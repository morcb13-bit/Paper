#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
行番号を辿って何段まで潜れるかを数える。標準ライブラリのみ。

  検定 SD1  p -> (p-3)/2 の鎖が、各段で生き残る歩数の上にだけ乗っているか。
            1 段でも外れれば NG。
  検定 SD2  段が上がるごとに本数が単調に減るか。増える段があれば NG。
  検定 SD3  （落ちるのが正しい検査）一手を p -> (p-1)/2 に変えたとき、
            同じ段数分布になるか。NG を返すのが正常。
"""

M = 10 ** 6

lim = 2 * M + 3
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

def chain(p, c=3):
    out = [p]
    while True:
        n = (out[-1] - c) // 2
        if n >= 3 and (out[-1] - c) % 2 == 0 and isprime(n):
            out.append(n)
        else:
            break
    return out

from collections import Counter
cnt = Counter()
longest = []
for p in range(3, M):
    if isprime(p):
        ch = chain(p)
        cnt[len(ch)] += 1
        if len(ch) >= 7:
            longest.append(ch)

print("鎖の長さの分布（p < %d）" % M)
for k in sorted(cnt):
    print("  %d 段 : %d 本" % (k, cnt[k]))

# SD1
sd1 = all(all(isprime(x) for x in chain(p)) for p in range(3, 20000) if isprime(p))
# SD2
ks = sorted(cnt)
sd2 = all(cnt[ks[i]] > cnt[ks[i + 1]] for i in range(len(ks) - 1))
# SD3（落ちるのが正しい）
cnt3 = Counter(len(chain(p, c=1)) for p in range(3, 200000) if isprime(p))
cnt0 = Counter(len(chain(p, c=3)) for p in range(3, 200000) if isprime(p))
sd3 = (cnt3 == cnt0)

print()
print("検定 SD1 : %s" % ("OK" if sd1 else "NG"))
print("検定 SD2 : %s" % ("OK" if sd2 else "NG"))
print("検定 SD3 : %s  （NG が正常）" % ("OK" if sd3 else "NG"))
if longest:
    print()
    print("7 段以上の鎖:")
    for ch in longest:
        print("  ", " -> ".join(str(x) for x in ch))
