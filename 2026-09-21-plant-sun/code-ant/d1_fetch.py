# 検定D1：絵札を読んで、見つけて、持ち帰る（ロボット一体）
# 倉庫＝②の床。家＝中心の五芒星。荷物 12 個を家以外の番地にでたらめに置き、それぞれに絵札（棒の 3 方向、位置はでたらめ）を貼る。
# 注文「札 X の荷物を一つ持ち帰れ」。ロボットは家から出て、まだ見ていない荷物のうちいちばん近いものへ匂い（段数）をたどって行き、
# 着いたら絵札を網膜で見て脳で読む。X と読めたら持ち上げて家の匂いをたどって戻る。違えば次へ。全部見て X が無ければ手ぶらで戻る。
# 脳：L1 の学習用 960・通票ありの模型（種 0〜9）。読み＝脳の中で最初に着いた出口（同時着・未着は「読めない」）。
# 基準（走らせる前に決めたもの）
#   合格  持ち帰った荷物が本当に札 X である割合が、でたらめに読む脳（3 方向を等しく選ぶ）を 5 ポイントを越えて上回る
#   NG    差がない
#   あわせて測る：持ち帰るまでの歩数、見に行った荷物の数、手ぶらで戻った回数
import sys; sys.argv = [sys.argv[0], '0', '0']
import json, random
from collections import Counter
import b13_amoeba10_gap as G
A = G.A
fl = A.make_floor(); DS = [A.scent(fl, e) for e in fl['ex']]
res = json.load(open('l1_result.json'))
HOME = fl['kind'].index('五芒星'); NP = 12; TRIALS = 400

def dist_from(src):
    d = {src: 0}; fr = [src]
    while fr:
        nf = []
        for u in fr:
            for w in fl['adj'][u]:
                if w not in d: d[w] = d[u] + 1; nf.append(w)
        fr = nf
    return d
DIST = [None] * fl['N']
def D(src):
    if DIST[src] is None: DIST[src] = dist_from(src)
    return DIST[src]

def card(cls, rng):
    xs = [fl['xy'][k][0] for k in fl['pix']]; r = max(abs(x) for x in xs)
    while True:
        img = A.draw_bar(fl, cls, rng.uniform(-r, r), rng.uniform(-r, r), 7.0, 1.8)
        if sum(img) >= 2: return img

def walk(pos, target):
    """匂い（段数）を一段ずつ下る。歩数を返す"""
    d = D(target); steps = 0
    while pos != target:
        pos = min(fl['adj'][pos], key=lambda j: d[j]); steps += 1
    return steps

def trial(read, rng):
    cells = [c for c in range(fl['N']) if c != HOME]
    spots = rng.sample(cells, NP); cls = [rng.randrange(3) for _ in spots]
    X = rng.randrange(3)
    if X not in cls: cls[rng.randrange(NP)] = X
    cards = [card(c, rng) for c in cls]
    pos = HOME; steps = 0; seen = set(); got = None
    while len(seen) < NP:
        i = min((i for i in range(NP) if i not in seen), key=lambda i: D(spots[i])[pos])
        steps += walk(pos, spots[i]); pos = spots[i]; seen.add(i)
        if read(cards[i]) == X: got = i; break
    steps += walk(pos, HOME)
    return got, cls[got] == X if got is not None else None, steps, len(seen)

def run(readers):
    out = {}
    for name, mk in readers.items():
        C = Counter()
        for s in range(10):
            read = mk(s); rng = random.Random(9000 + s)
            for _ in range(TRIALS):
                got, ok, steps, nseen = trial(read, rng)
                C['試行'] += 1; C['歩数'] += steps; C['見た荷物'] += nseen
                if got is None: C['手ぶら'] += 1
                else: C['持ち帰り'] += 1; C['正しい'] += int(ok)
        out[name] = C
        print(f"{name}：試行 {C['試行']}  持ち帰り {C['持ち帰り']}  うち本当に札X {C['正しい']}  手ぶら {C['手ぶら']}  "
              f"歩数の合計 {C['歩数']}  見に行った荷物の合計 {C['見た荷物']}", flush=True)
    return out

def brain(s):
    st, sw = next(r['model'] for n, ss, r in res if n == 960 and ss == s)
    def read(img):
        _, f = A.simulate(fl, DS, img, st, sw, True)
        return f if f in (0, 1, 2) else None
    return read
def randread(s):
    r = random.Random(500 + s); return lambda img: r.randrange(3)
if __name__ == "__main__":
    run({'脳（L1 960）': brain, 'でたらめに読む': randread})
