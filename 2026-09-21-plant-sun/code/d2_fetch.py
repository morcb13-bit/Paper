# 検定D2：確かになるまで持ち上げない。荷物に着いたら網膜の位置をずらして同じ絵札を k 回見て、k 回すべて X と読めたときだけ持ち上げる。
# 倉庫・荷物 12 個・注文・試行 4000（脳 10 × 400）は D1 と同じ。脳は L1 の学習用 960・通票ありの模型（学習はやり直さない）。
# 一回の見直し：絵札の中心を x, y それぞれ ±1.0 の範囲でずらして描く（点が 2 未満なら「読めない」＝X ではない）。
# 基準（走らせる前に決めたもの）
#   合格  持ち帰った荷物のうち別の札が 1% 以下（4000 試行で 40 件以下）
#   負の対照  でたらめに読む脳では、別の札の割合が下がらない
#   NG  k = 7 でも別の札が 1% を越える
#   あわせて測る：手ぶらで戻った回数、一試行あたりの歩数と見た回数
import sys; sys.argv = [sys.argv[0], '0', '0']
import random
from collections import Counter
import d1_fetch as D1
A = D1.A; fl = D1.fl; HOME = D1.HOME; NP = D1.NP
r_ret = max(abs(fl['xy'][k][0]) for k in fl['pix'])

def card(cls, rng):
    while True:
        cx, cy = rng.uniform(-r_ret, r_ret), rng.uniform(-r_ret, r_ret)
        if sum(A.draw_bar(fl, cls, cx, cy, 7.0, 1.8)) >= 2: return (cls, cx, cy)

def trial(read, k, rng, grng):
    cells = [c for c in range(fl['N']) if c != HOME]
    spots = rng.sample(cells, NP); cls = [rng.randrange(3) for _ in spots]
    X = rng.randrange(3)
    if X not in cls: cls[rng.randrange(NP)] = X
    cards = [card(c, rng) for c in cls]
    pos = HOME; steps = 0; seen = set(); got = None; looks = 0
    while len(seen) < NP:
        i = min((i for i in range(NP) if i not in seen), key=lambda i: D1.D(spots[i])[pos])
        steps += D1.walk(pos, spots[i]); pos = spots[i]; seen.add(i)
        c, cx, cy = cards[i]; ok = True
        for _ in range(k):
            looks += 1
            img = A.draw_bar(fl, c, cx + grng.uniform(-1, 1), cy + grng.uniform(-1, 1), 7.0, 1.8)
            if sum(img) < 2 or read(img) != X: ok = False; break
        if ok: got = i; break
    steps += D1.walk(pos, HOME)
    return got, (cls[got] == X) if got is not None else None, steps, looks

if __name__ == "__main__":
    for name, mk in (('脳（L1 960）', D1.brain), ('でたらめに読む', D1.randread)):
        for k in (1, 3, 5, 7):
            C = Counter()
            for s in range(10):
                read = mk(s); rng = random.Random(9000 + s); grng = random.Random(7700 + s)
                for _ in range(D1.TRIALS):
                    got, ok, steps, looks = trial(read, k, rng, grng)
                    C['試行'] += 1; C['歩数'] += steps; C['見た回数'] += looks
                    if got is None: C['手ぶら'] += 1
                    else: C['持ち帰り'] += 1; C['別の札'] += int(not ok)
            print(f"{name} k={k}：持ち帰り {C['持ち帰り']}  別の札 {C['別の札']}  手ぶら {C['手ぶら']}  歩数の合計 {C['歩数']}  見た回数の合計 {C['見た回数']}", flush=True)
