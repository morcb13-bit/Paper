# CM3 の下ごしらえ：compound.py と同一の手順で adj / O / ST / CEN を作り、pickle に落とす
import json, math, pickle, time
from collections import defaultdict
import b13_chain_units as U

t0 = time.time()
cells = {tuple(int(x) for x in k.split(",")): a
         for k, a in json.load(open("carrier_1245.json"))["cells"].items()}
adj = defaultdict(set)
for q, a in cells.items():
    vs = [U.zadd(q, U.zt(a + 2 * i)) for i in range(5)]
    for i in range(5):
        u, w = vs[i], vs[(i + 1) % 5]
        adj[u].add(w); adj[w].add(u)
XY = {p: U.xy(p) for p in adj}
print("頂点 %d  経過 %.1fs" % (len(adj), time.time() - t0), flush=True)

O = {}
for cur in adj:
    for prev in adj[cur]:
        outs = [w for w in adj[cur] if w != prev]
        if not outs:
            O[(prev, cur)] = {}; continue
        back = math.atan2(XY[prev][1] - XY[cur][1], XY[prev][0] - XY[cur][0])
        ws = [w for _, w in sorted(((math.atan2(XY[w][1] - XY[cur][1],
                                                XY[w][0] - XY[cur][0]) - back) % (2 * math.pi), w)
                                   for w in outs)]
        O[(prev, cur)] = {0: ws[0]} if len(ws) == 1 else (
            {1: ws[0], -1: ws[1]} if len(ws) == 2 else {1: ws[0], 0: ws[1], -1: ws[2]})
print("有向辺 %d  経過 %.1fs" % (len(O), time.time() - t0), flush=True)

ST = []; CEN = []
for area, cyc in U.gaps(cells):
    if abs(area - 2.9389) < 0.01:
        c = (0, 0, 0, 0)
        for p in cyc: c = U.zadd(c, p)
        if any(len(adj[p]) < 2 for p in cyc): continue
        CEN.append(c); ST.append([(p, cyc[(i + 1) % 10]) for i, p in enumerate(cyc)])
print("五芒星 %d 個  経過 %.1fs" % (len(CEN), time.time() - t0), flush=True)

pickle.dump({"O": O, "ST": ST, "CEN": CEN}, open("cm3_base.pkl", "wb"), 4)
print("cm3_base.pkl 書き出し  経過 %.1fs" % (time.time() - t0), flush=True)
