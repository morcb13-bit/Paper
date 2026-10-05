# 測定（判定なし）：通票ありで、待避・待機・乗り換えが何回起きたか。番地ごとにも数える。
# 模型は b13_amoeba10_first.py と同じ種・同じ学習（通票あり）で作り、見せていない試験 400 枚 × 10 回で数える。
import random, sys
from collections import Counter
import b13_amoeba10_first as A

NTR = int(sys.argv[1]) if len(sys.argv) > 1 else 240
fl = A.make_floor(); DS = [A.scent(fl, e) for e in fl['ex']]; N = fl['N']; L = 3
up = [[sum(1 for j in fl['adj'][i] if DS[l][j] > DS[l][i]) for l in range(L)] for i in range(N)]  # 濃い方の隣の数（層ごと）

C = Counter(); cell = {k: Counter() for k in ('step', 'first', 'side', 'wait_no2', 'wait_both', 'swap', 'stuck')}
def run(img, start, swap, T=120):
    rob = [[fl['pix'][k], start[k]] for k in range(len(img)) if img[k]]
    exset = set(fl['ex'])
    for t in range(T):
        if not rob: break
        occ = set(r[0] for r in rob); claimed = set(); nxt = []
        for pos, lay in rob:
            D = DS[lay]; nb = sorted(fl['adj'][pos], key=lambda j: -D[j])
            C['step'] += 1; cell['step'][pos] += 1
            ups = [c for c in nb[:2] if D[c] > D[pos]]
            pos2 = None
            for n, c in enumerate(ups):
                if c in claimed or (c in occ and c not in exset): continue
                claimed.add(c); pos2 = c
                kind = 'first' if n == 0 else 'side'; C[kind] += 1; cell[kind][pos] += 1
                if n == 0 and len(ups) == 2: C['first_with2'] += 1
                break
            if pos2 is None:
                kind = 'wait_no2' if len(ups) == 1 else ('wait_both' if len(ups) == 2 else 'stuck')
                C[kind] += 1; cell[kind][pos] += 1
                if swap[pos]: lay = (lay + 1) % L; C['swap'] += 1; cell['swap'][pos] += 1
                pos2 = pos
            if pos2 in exset: C['arrive'] += 1
            else: nxt.append([pos2, lay])
        rob = nxt
    C['never'] += len(rob)

for s in range(10):
    rng = random.Random(8000 + s)
    Xtr, ytr = A.make_data(fl, NTR, rng); Xte, yte = A.make_data(fl, 400, rng)
    st, sw = A.learn(fl, DS, Xtr, ytr, True, rng)
    C['swaptags'] += sum(sw)
    for img in Xte: C['robots'] += sum(img); run(img, st, sw)

print(f"学習用 {NTR}、試験 4000 枚、ロボット {C['robots']} 体、一体一刻の手 {C['step']}")
for k, name in (('first', '一番濃い隣へ進む'), ('side', '待避（二番目へ）'), ('wait_both', '待機（二本とも塞がる）'),
                ('wait_no2', '待機（濃い隣が一本しかなく塞がる）'), ('stuck', '止まる（濃い隣がない）'), ('swap', 'うち乗り換え')):
    print(f"  {name:22s} {C[k]}")
print(f"  一番濃い隣へ進んだうち、二本目もあった手 {C['first_with2']}")
print(f"  出口に着いた {C['arrive']}、120 刻で着かなかった {C['never']}、乗り換えの札（10 模型の合計） {C['swaptags']}")
# 床の構造：濃い隣の数ごとの番地数（層ごと）
for l in range(L): print(f"  層{l}：濃い隣の数ごとの番地数 {sorted(Counter(min(u[l],2) for u in up).items())}")
for k, name in (('side', '待避'), ('wait_both', '待機（二本とも）'), ('wait_no2', '待機（一本）'), ('swap', '乗り換え')):
    tot = sum(cell[k].values()); top = cell[k].most_common(8)
    print(f"番地ごと {name}：起きた番地 {len(cell[k])}、上位8番地で {sum(v for _, v in top)}/{tot}  "
          + "  ".join(f"{i}(扇{fl['fan'][i]},隣{len(fl['adj'][i])},{'網膜' if i in fl['pix'] else '床'}):{v}" for i, v in top))
