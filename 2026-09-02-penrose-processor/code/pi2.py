# 検定PI2（符号つきの和の対照）
#
#   §2-15 の |Σω^S| = 7.5498（最短経路60本）が装置の量か、本数のばらつきかを分ける。
#
#   出すもの
#     (0) 本番の再現（60本・S分布・|Σω^S|・Σ(−1)^S）
#     (1) 無位相の目安 √n（位相が無相関なら |Σ| は √n のあたりに来る）
#     (2) 対照A 桁の付け替え：盤面も経路も60本もそのまま、頂点への桁の割り当てだけを
#         個数を保って入れ替える。位置と桁の結び付きだけが壊れる
#     (3) 対照B 独立引き：各経路の S を、39歩ぶんの桁を個数比で独立に引いて作る
#     (4) §6 が指定した配置換えグラフが d=39 を作れるかを確かめる
#
#   OK   本番の |Σ| が対照AとBの分布から外れる（下側 5% 未満など）＝装置の量
#   NG   対照の分布の真ん中に来る＝本数のばらつき。符号の筋を切る
#   必ず落ちる設定  桁を全部 0 にすると |Σ| = n（打ち消しゼロ）
import json, math, cmath, random, statistics
from collections import defaultdict, Counter, deque
import b13_chain_units as U

cells = {tuple(int(x) for x in k.split(",")): a
         for k, a in json.load(open("carrier_1245.json"))["cells"].items()}
adj = defaultdict(set)
for q, a in cells.items():
    vs = [U.zadd(q, U.zt(a + 2 * i)) for i in range(5)]
    for i in range(5):
        u, w = vs[i], vs[(i + 1) % 5]
        adj[u].add(w); adj[w].add(u)
XY = {p: U.xy(p) for p in adj}
DIG = {p: {2: -1, 3: 0, 4: 1}[len(adj[p])] for p in adj}
O = {}
for cur in adj:
    for prev in adj[cur]:
        outs = [w for w in adj[cur] if w != prev]
        if not outs: O[(prev, cur)] = {}; continue
        back = math.atan2(XY[prev][1] - XY[cur][1], XY[prev][0] - XY[cur][0])
        ws = [w for _, w in sorted(((math.atan2(XY[w][1] - XY[cur][1],
                                                XY[w][0] - XY[cur][0]) - back) % (2 * math.pi), w)
                                   for w in outs)]
        O[(prev, cur)] = {0: ws[0]} if len(ws) == 1 else (
            {1: ws[0], -1: ws[1]} if len(ws) == 2 else {1: ws[0], 0: ws[1], -1: ws[2]})
CEN = []; RIM = []
for area, cyc in U.gaps(cells):
    if abs(area - 2.9389) < 0.01:
        c = (0, 0, 0, 0)
        for p in cyc: c = U.zadd(c, p)
        if any(len(adj[p]) < 2 for p in cyc): continue
        CEN.append((U.xy(c)[0] / 10, U.xy(c)[1] / 10)); RIM.append(cyc)
ia = min(range(len(CEN)), key=lambda i: abs(math.hypot(*CEN[i]) - 12.0))
ib = min(range(len(CEN)), key=lambda j: abs(math.hypot(CEN[j][0] - CEN[ia][0],
                                                       CEN[j][1] - CEN[ia][1]) - 40))
DA = math.hypot(CEN[ib][0] - CEN[ia][0], CEN[ib][1] - CEN[ia][1])
print("五芒星 %d 個   A=%s  B=%s  実距離 %.4f"
      % (len(CEN), tuple(round(x, 3) for x in CEN[ia]), tuple(round(x, 3) for x in CEN[ib]), DA))
SA = [(RIM[ia][i], RIM[ia][(i + 1) % 10]) for i in range(10)]
SB = set(RIM[ib])

# --- 有向辺の最短距離と、最短経路の全列挙 ---
dist = {e: 0 for e in SA}
q = deque(SA); hit = None
while q:
    e = q.popleft(); d = dist[e]
    if hit is not None and d >= hit: continue
    for x, nxt in O[e].items():
        ne = (e[1], nxt)
        if ne not in dist:
            dist[ne] = d + 1; q.append(ne)
            if nxt in SB and hit is None: hit = d + 1
print("最短の歩数 d = %d" % hit)

ends = [e for e, d in dist.items() if d == hit and e[1] in SB]
pre = defaultdict(list)
for e, d in dist.items():
    for x, nxt in O[e].items():
        ne = (e[1], nxt)
        if dist.get(ne) == d + 1: pre[ne].append(e)
paths = []
def back(e, acc):
    if dist[e] == 0:
        paths.append([e] + acc[::-1] if False else [e] + list(reversed(acc))); return
    for p in pre[e]: back(p, acc + [e])
for e in ends: back(e, [])
paths = [p for p in paths]
print("最短経路 %d 本" % len(paths))

def Sof(path, dig):
    return sum(dig[e[1]] for e in path[1:])
S0 = [Sof(p, DIG) for p in paths]
cnt = Counter(S0); n = len(S0)
w3 = cmath.exp(2j * math.pi / 3)
def amp(Ss): return abs(sum(w3 ** s for s in Ss))
print("桁和 S の分布 %s" % dict(sorted(cnt.items())))
print("本番  |Σω^S| = %.4f   Σ(−1)^S = %d   （打ち消しなしなら %d）"
      % (amp(S0), sum((-1) ** s for s in S0), n))
print("必ず落ちる設定（桁を全部0）  |Σω^S| = %.4f" % amp([0] * n))
print("無位相の目安 √n = %.4f" % math.sqrt(n))
print()

# --- 対照A：桁の付け替え ---
V = list(adj)
labels = [DIG[p] for p in V]
rng = random.Random(20260904)
vals = []
for _ in range(2000):
    rng.shuffle(labels)
    dig = dict(zip(V, labels))
    vals.append(amp([Sof(p, dig) for p in paths]))
vals.sort()
lo = sum(1 for v in vals if v <= amp(S0)) / len(vals)
print("対照A 桁の付け替え（2000回）")
print("  |Σ| 中央 %.4f   5%%点 %.4f   95%%点 %.4f   本番 %.4f は下から %.1f%%"
      % (statistics.median(vals), vals[len(vals) // 20], vals[-len(vals) // 20], amp(S0), 100 * lo))

# --- 対照B：独立引き ---
degc = Counter(len(adj[p]) for p in adj)
pool = [-1] * degc[2] + [0] * degc[3] + [1] * degc[4]
vals2 = []
for _ in range(2000):
    Ss = [sum(rng.choice(pool) for _ in range(hit)) for _ in range(n)]
    vals2.append(amp(Ss))
vals2.sort()
lo2 = sum(1 for v in vals2 if v <= amp(S0)) / len(vals2)
print("対照B 独立引き（2000回）")
print("  |Σ| 中央 %.4f   5%%点 %.4f   95%%点 %.4f   本番 %.4f は下から %.1f%%"
      % (statistics.median(vals2), vals2[len(vals2) // 20], vals2[-len(vals2) // 20], amp(S0), 100 * lo2))
print()

# --- 配置換えグラフは d=39 を作れるか ---
deg = {p: len(adj[p]) for p in adj}
stubs = []
for p, d in deg.items(): stubs += [p] * d
rng.shuffle(stubs)
radj = defaultdict(set)
for i in range(0, len(stubs) - 1, 2):
    u, v = stubs[i], stubs[i + 1]
    if u != v: radj[u].add(v); radj[v].add(u)
src = next(iter(radj))
seen = {(src, x): 0 for x in radj[src]}
q = deque(seen); far = 0
while q:
    e = q.popleft(); d = seen[e]; far = max(far, d)
    for nxt in radj[e[1]]:
        if nxt == e[0]: continue
        ne = (e[1], nxt)
        if ne not in seen: seen[ne] = d + 1; q.append(ne)
print("配置換えグラフ（同じ次数列）の非後退の到達距離 最大 %d 歩（本番は %d 歩）" % (far, hit))
print("→ d=%d の対を作れないので、§6 の指定どおりの比較は構成できない" % hit)
