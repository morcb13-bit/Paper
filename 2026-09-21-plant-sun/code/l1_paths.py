# 測定（判定なし）：L1 の通票ありの模型（学習用 240 と 960、種 10 ずつ）を固定の試験 4000 枚に当て、経路から数える
import sys; sys.argv = [sys.argv[0], '0', '0']
import json, random
from collections import Counter
import b13_amoeba10_gap as G
A = G.A
fl = A.make_floor(); DS = [A.scent(fl, e) for e in fl['ex']]
Xte, yte = A.make_data(fl, 4000, random.Random(12345))
res = json.load(open('l1_result.json'))
pixset = set(fl['pix'])
def run(img, st, sw, T=120):
    rob = [[fl['pix'][k], st[k]] for k in range(len(img)) if img[k]]
    exset = set(fl['ex']); arr = [0] * 3; first = None; ft = None; waits = swaps = 0
    for t in range(T):
        if not rob: break
        occ = set(r[0] for r in rob); claimed = set(); nxt = []
        for pos, lay in rob:
            D = DS[lay]; nb = sorted(fl['adj'][pos], key=lambda j: -D[j]); pos2 = None
            for c in nb[:2]:
                if D[c] <= D[pos]: break
                if c in claimed or (c in occ and c not in exset): continue
                claimed.add(c); pos2 = c; break
            if pos2 is None:
                pos2 = pos; waits += 1
                if sw[pos]: lay = (lay + 1) % 3; swaps += 1
            if pos2 in exset: arr[fl['ex'].index(pos2)] += 1
            else: nxt.append([pos2, lay])
        if first is None:
            got = {j for j in range(3) if arr[j]}
            if got: first = got.pop() if len(got) == 1 else -1; ft = t + 1
        rob = nxt
    return first, ft, arr, waits, swaps
for n in (240, 960):
    S = Counter(); where = Counter(); ftc = Counter(); conv = Counter()
    for nn, s, r in res:
        if nn != n: continue
        st, sw = r['model']
        S['札A'] += st.count(0); S['札B'] += st.count(1); S['札C'] += st.count(2); S['乗り換え札'] += sum(sw)
        for i, v in enumerate(sw):
            if v: where[('網膜' if i in pixset else '床') + '・' + fl['kind'][i]] += 1
        for img, yy in zip(Xte, yte):
            f, ft, arr, w, sp = run(img, st, sw)
            S['画像'] += 1; S['待機'] += w; S['乗り換え'] += sp; S['ロボット'] += sum(img)
            if f == yy: S['当たり'] += 1; ftc[ft] += 1; S['当たりの刻'] += ft
            else:
                S['外れ'] += 1
                if f == -1: S['同時'] += 1
                elif f is None: S['未着'] += 1
                k = sum(1 for a in arr if a)
                conv[k] += 1
                if f is not None and f >= 0 and arr[f] == sum(arr) and sum(arr) == sum(img): S['外れ・全員が同じ出口'] += 1
                if arr[yy] == 0: S['外れ・正解の出口に誰も着かない'] += 1
    print(f"■ 学習用 {n}（模型 10、試験 4000×10）")
    print(f"  最初の層の札（網膜 40×10）：A {S['札A']}  B {S['札B']}  C {S['札C']}   乗り換えの札 {S['乗り換え札']}")
    print(f"  乗り換えの札の場所：{dict(where.most_common())}")
    print(f"  一枚あたり（×100）：ロボット {100*S['ロボット']//S['画像']}  待機 {100*S['待機']//S['画像']}  乗り換え {100*S['乗り換え']//S['画像']}")
    print(f"  当たり {S['当たり']}：最初に着くまでの刻の合計 {S['当たりの刻']}（一枚あたり×100 {100*S['当たりの刻']//S['当たり']}）  刻の分布 {sorted(ftc.items())[:12]}")
    print(f"  外れ {S['外れ']}：同時着 {S['同時']}  未着 {S['未着']}  全員が同じ出口 {S['外れ・全員が同じ出口']}  正解の出口に誰も着かない {S['外れ・正解の出口に誰も着かない']}  着いた出口の種類数 {sorted(conv.items())}")
