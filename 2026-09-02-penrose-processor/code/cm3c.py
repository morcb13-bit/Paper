# 検定CM3 最終
#   前の走りで分かったこと：基線ごとに「平行移動して担体に残る星」の数が大きく違う。
#   10.5784 は 411個中 1個しか残らないので、|W∩(W+b)|=1 は幾何ではなく縁の産物。
#   よって基線は「残る星」で正規化し、必ず同大の無作為Wと並べて読む。
#
#   OK   縁を揃えた上で、三眼の中央が二眼より小さくなり、かつ対照より小さい
#   NG   三眼が二眼と変わらない、または全滅する、または対照と同じ
#   必ず落ちる設定  c=0 で三眼は二眼と一致
import pickle, math, random
from collections import defaultdict
import b13_chain_units as U

N = 12
B = pickle.load(open("cm3_base.pkl", "rb"))
O, ST, CEN = B["O"], B["ST"], B["CEN"]
n = len(CEN); key = {c: i for i, c in enumerate(CEN)}

def flows(w, starts):
    for a, b in starts:
        prev, cur = a, b; ok = True
        for x in w:
            o = O[(prev, cur)]
            if x not in o: ok = False; break
            prev, cur = cur, o[x]
        if ok: return True
    return False

rng = random.Random(20260903)
words = [tuple(rng.choice((-1, 0, 1)) for _ in range(N)) for _ in range(600)]
one = [S for S in (set(i for i in range(n) if flows(w, ST[i])) for w in words) if S]
print("検定CM3 最終  五芒星 %d 個  語長 %d  受理された語 %d" % (n, N, len(one)))
rng2 = random.Random(11); allst = list(range(n))
rnd = [set(rng2.sample(allst, len(S))) for S in one]

c0 = CEN[0]; seen = {}
for j in range(n):
    d = U.zsub(CEN[j], c0); L = math.hypot(*U.xy(d)) / 10.0
    if 0 < L < 26: seen[d] = round(L, 4)
BASE = sorted(seen.items(), key=lambda kv: (kv[1], kv[0]))

def shift(S, d):
    out = set()
    for i in S:
        j = key.get(U.zadd(CEN[i], d))
        if j is not None: out.add(j)
    return out
ALL = set(range(n))
def med(sets):
    ne = [len(x) for x in sets if x]
    if not ne: return (0, 0, 0)
    ne.sort(); return (len(ne), ne[len(ne) // 2], sum(1 for x in ne if x == 1))

print()
print("二眼（本番と対照を並べる）")
print("  基線 b     残る星  本番:語/中央/1点   対照:語/中央/1点")
rows = []
for d, L in BASE:
    keep = len(shift(ALL, d))
    a = med([S & shift(S, d) for S in one])
    b_ = med([S & shift(S, d) for S in rnd])
    rows.append((L, d, keep, a, b_))
for L, d, keep, a, b_ in rows:
    print("  %8.4f   %4d    %4d /%4d /%3d      %4d /%4d /%3d" % ((L, keep) + a + b_))

# 縁が効きすぎない基線だけを残す
ok = [r for r in rows if r[2] >= 80]
ok.sort(key=lambda r: (r[3][1], -r[3][0]))
L, bd, keep, a, b_ = ok[0]
print()
print("採る基線 b = %.4f（残る星 %d、本番の中央 %d に対し対照 %d）" % (L, keep, a[1], b_[1]))
two = [S & shift(S, bd) for S in one]
twoR = [S & shift(S, bd) for S in rnd]

print()
print("三眼（b 固定、c を全基線で走査）")
ZERO = (0, 0, 0, 0)
print("  必ず落ちる設定 c=0 :", "OK" if [x & shift(x, ZERO) for x in two] == two else "NG")
print("  基線 c     残る星  本番:語/中央/1点   対照:語/中央/1点")
best = None
for d, L2 in BASE:
    keep2 = len(shift(ALL, d))
    if keep2 < 80: continue
    a2 = med([x & shift(x, d) for x in two])
    b2 = med([x & shift(x, d) for x in twoR])
    print("  %8.4f   %4d    %4d /%4d /%3d      %4d /%4d /%3d" % ((L2, keep2) + a2 + b2))
    if best is None or (a2[0], -a2[1]) > (best[1][0], -best[1][1]): best = (L2, a2, d)
print()
if best:
    L2, a2, d = best
    three = [x & shift(x, d) for x in two]
    surv = [(len(p), len(q)) for p, q in zip(two, three) if q]
    print("c = %.4f で残った語 %d 本：二眼|W| 中央 %d → 三眼|W| 中央 %d"
          % (L2, len(surv),
             sorted(s[0] for s in surv)[len(surv)//2] if surv else 0,
             sorted(s[1] for s in surv)[len(surv)//2] if surv else 0))
    print("三眼で消えた語 %d / %d" % (sum(1 for p, q in zip(two, three) if p and not q),
                                      sum(1 for p in two if p)))
