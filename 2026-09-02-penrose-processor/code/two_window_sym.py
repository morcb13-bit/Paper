#  検定TW1-S  5回対称の担体で測り直す
#
#      監督の図（2026-09-02 提示）は、鎖に積んだ担体ではなく5回対称の担体である。
#      大きなタイリングの頂点＝円環中心。窓B はその点集合。
#      事前登録・対照・必ず落ちる設定は two_window.py（v2）と同一。担体だけを差し替える。

import json, math, random, sys
import b13_chain_units as U

src = open("two_window.py").read()
exec(src[src.index("TOL = 0.30"):src.index("def selftest")])   # 物差しだけを取り込む

d = json.load(open("rings_integer.json"))
RINGS = [tuple(r) for r in d["rings"]]
CELLS = [tuple(json.loads(k)) if k.startswith("[") else tuple(eval(k)) for k in d["cells"]]

A = [U.xy(q) for q in CELLS]
B = [U.xy(r) for r in RINGS]
print("担体（5回対称）：五角形 %d 枚（窓A）／円環中心 %d 個（窓B）" % (len(A), len(B)))

CTR = min(CELLS, key=lambda q: U.xy(q)[0] ** 2 + U.xy(q)[1] ** 2)
SHIFTS = []
for q in CELLS:
    x, y = U.xy(U.zsub(q, CTR))
    if x * x + y * y <= 36:
        SHIFTS.append((x, y))
print("試すずれ %d 通り（長さ6以内）\n" % len(SHIFTS))

d_rc = [min(math.dist(r, c) for c in A) for r in B]
print("■ 検定TW0 円環中心から最も近い五角形中心まで %.6f 〜 %.6f ／ 一致 %d 個\n"
      % (min(d_rc), max(d_rc), sum(1 for t in d_rc if t < 1e-9)))

print("■ 検定TW1-S")
rA = report("窓A 五角形", A, SHIFTS)
rB = report("窓B 円環中心", B, SHIFTS)
rU = report("A∪B 両方", A + B, SHIFTS)

rnd = random.Random(13)
v1 = []
for _ in range(NBAND):
    th = rnd.uniform(0, 2 * math.pi); rr = rnd.uniform(0.6, 1.4)
    o = (rr * math.cos(th), rr * math.sin(th))
    v1.append(profile(A + [(p[0] + o[0], p[1] + o[1]) for p in B], SHIFTS)[2])
v1.sort()
xs = [p[0] for p in A]; ys = [p[1] for p in A]
v2 = []
for _ in range(NBAND):
    R = [(rnd.uniform(min(xs), max(xs)), rnd.uniform(min(ys), max(ys))) for _ in B]
    v2.append(profile(A + R, SHIFTS)[2])
v2.sort()
print("\n  対照1 円環中心を格子外へずらす  帯 %.4f 〜 %.4f" % (v1[0], v1[-1]))
print("  対照2 一様乱数に置き換える      帯 %.4f 〜 %.4f" % (v2[0], v2[-1]))

print("\n  判定：", end="")
if rU >= rA - 1e-12:
    print("NG。A∪B %.4f が A 単独 %.4f を下回らない" % (rU, rA))
elif rU < v1[0] and rU < v2[0]:
    print("OK。A∪B %.4f は A 単独 %.4f より低く、帯の下端 %.4f / %.4f も下回る"
          % (rU, rA, v1[0], v2[0]))
else:
    print("NG。A∪B %.4f は帯の中（%.4f〜%.4f / %.4f〜%.4f）"
          % (rU, v1[0], v1[-1], v2[0], v2[-1]))

rS = profile(A + A, SHIFTS)[2]
print("\n  必ず落ちる設定（B を A 自身に）… %.4f ／ 窓A単独 %.4f → %s"
      % (rS, rA, "一致" if abs(rS - rA) < 1e-12 else "不一致・物差しが壊れている"))

ra = {d_: n / len(A) for n, d_ in selfmatch(A, SHIFTS)}
rb = {d_: n / len(B) for n, d_ in selfmatch(B, SHIFTS)}
print("\n  台地のずれベクトル（窓A の上位5本）")
for n, d_ in sorted(selfmatch(A, SHIFTS), reverse=True)[1:6]:
    print("    (%7.3f,%7.3f) 長さ%.3f  窓A %.3f / 窓B %.3f"
          % (d_[0], d_[1], math.hypot(*d_), ra[d_], rb[d_]))
