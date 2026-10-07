# 経路の記録：L1 の学習用 960・種 0 の通票ありの模型で、試験画像の数枚について、ロボットごとの経路を書き出す
import sys; sys.argv = [sys.argv[0], '0', '0']
import json, random
import b13_amoeba10_gap as G
A = G.A
fl = A.make_floor(); DS = [A.scent(fl, e) for e in fl['ex']]
Xte, yte = A.make_data(fl, 4000, random.Random(12345))
res = json.load(open('l1_result.json'))
st, sw = next(r['model'] for n, s, r in res if n == 960 and s == 0)
name = lambda p: f"{p}" + ("" if fl['kind'][p] == '五角形' else f"({fl['kind'][p]})")
def trace(img, T=120):
    rob = [[fl['pix'][k], st[k], [f"網膜{fl['pix'][k]}・層{'ABC'[st[k]]}"]] for k in range(len(img)) if img[k]]
    exset = set(fl['ex']); first = None; done = []
    for t in range(T):
        if not rob: break
        occ = set(r[0] for r in rob); claimed = set(); nxt = []
        for pos, lay, log in rob:
            D = DS[lay]; nb = sorted(fl['adj'][pos], key=lambda j: -D[j]); pos2 = None
            for n, c in enumerate(nb[:2]):
                if D[c] <= D[pos]: break
                if c in claimed or (c in occ and c not in exset): continue
                claimed.add(c); pos2 = c; log.append(("待避→" if n else "→") + name(c)); break
            if pos2 is None:
                pos2 = pos
                if sw[pos]: lay = (lay + 1) % 3; log.append(f"待機・乗り換え→層{'ABC'[lay]}")
                else: log.append("待機")
            if pos2 in exset:
                log.append(f"出口{'ABC'[fl['ex'].index(pos2)]}（刻{t+1}）"); done.append(log)
                if first is None: first = (t, fl['ex'].index(pos2))
            else: nxt.append([pos2, lay, log])
        rob = nxt
    return done + [r[2] + ["着かない"] for r in rob]
shown = 0
for i, (img, yy) in enumerate(zip(Xte, yte)):
    if sum(img) < 3: continue
    logs = trace(img)
    print(f"■ 試験 {i}：正解 {['0°','60°','120°'][yy]}（出口{'ABC'[yy]}）")
    for lg in logs: print("   " + " ".join(lg))
    shown += 1
    if shown == 3: break
