# 検定D4：見つからなければ巣に戻る（ロボット一体）
# 倉庫・巣・荷物 12 個・脳（L1 の 960 枚、種 0〜9）は D1 と同じ。注文「札 X を持ってこい」。注文の 2 割は倉庫に札 X が無い。
# 燃料 E 歩（一回の外出）。巣の匂い＝巣までの段数 h(番地)。
# 脳が決めること（整数の比較だけ）
#   引き返す：次の一歩の先 q で、残りの燃料 − 1 が h(q) を下回るなら、探すのをやめて巣へ戻る（外出が終わる。巣で燃料を満たして次の外出）
#   覚える：見た荷物の番地と読んだ札を記憶に残し、外出をまたいで持つ
#   次の外出：記憶に無い荷物のうち一番近いものへ。全部見て X と読めたものが無ければ、巣で「無い」と返して終わる
# 対照：記憶なし（外出のたびに記憶を消す。外出 30 回で打ち切り）／引き返しなし（燃料が尽きるまで探す）
# 基準（走らせる前に決めたもの）
#   D4-1 燃料が尽きて巣に戻れなかった回数が 4000 注文で 0
#   D4-2 X が無い注文はすべて有限回の外出で「無い」と返して巣に戻る
#   D4-3 記憶ありは記憶なしより、持ち帰るまでの外出の回数が少ない
#   負の対照 引き返しなしでは、戻れない回数が 0 でなくなる
import sys; sys.argv = [sys.argv[0], '0', '0']
import random
from collections import Counter
import d1_fetch as D1
A = D1.A; fl = D1.fl; HOME = D1.HOME; NP = D1.NP
h = D1.D(HOME); E = 2 * max(h.values()) + 4

def order(read, rng, memory=True, turn=True, cap=30):
    cells = [c for c in range(fl['N']) if c != HOME]
    spots = rng.sample(cells, NP); X = rng.randrange(3)
    absent = rng.randrange(5) == 0
    cls = [rng.choice([c for c in range(3) if c != X]) if absent else rng.randrange(3) for _ in spots]
    if not absent and X not in cls: cls[rng.randrange(NP)] = X
    cards = [D1.card(c, rng) for c in cls]
    mem = {}; trips = 0; steps = 0
    while trips < cap:
        trips += 1; fuel = E; pos = HOME
        if not memory: mem = {}
        while True:
            todo = [i for i in range(NP) if i not in mem]
            if not todo: break
            i = min(todo, key=lambda i: D1.D(spots[i])[pos]); d = D1.D(spots[i]); back = False
            while pos != spots[i]:
                q = min(fl['adj'][pos], key=lambda j: d[j])
                if turn and fuel - 1 < h[q]: back = True; break        # 一歩進むと、残りの燃料が巣までの段数を下回るなら引き返す
                if fuel == 0: return dict(r='戻れない', trips=trips, steps=steps)
                pos = q; fuel -= 1; steps += 1
            if back: break
            mem[i] = read(cards[i])
            if mem[i] == X:
                while pos != HOME:
                    if fuel == 0: return dict(r='戻れない', trips=trips, steps=steps)
                    pos = min(fl['adj'][pos], key=lambda j: h[j]); fuel -= 1; steps += 1
                return dict(r='持ち帰り', ok=cls[i] == X, absent=absent, trips=trips, steps=steps)
        while pos != HOME:
            if fuel == 0: return dict(r='戻れない', trips=trips, steps=steps)
            pos = min(fl['adj'][pos], key=lambda j: h[j]); fuel -= 1; steps += 1
        if all(i in mem for i in range(NP)): return dict(r='無い', absent=absent, trips=trips, steps=steps)
    return dict(r='打ち切り', absent=absent, trips=trips, steps=steps)

if __name__ == "__main__":
    print(f"巣から一番遠い番地まで {max(h.values())} 段、燃料 E = {E}")
    for name, kw in (('記憶あり・引き返しあり', {}), ('記憶なし', dict(memory=False)), ('引き返しなし', dict(turn=False))):
        C = Counter()
        for s in range(10):
            read = D1.brain(s); rng = random.Random(9000 + s)
            for _ in range(400):
                o = order(read, rng, **kw); C[o['r']] += 1; C['外出'] += o['trips']; C['歩数'] += o['steps']
                if o['r'] == '持ち帰り':
                    C['持ち帰りの外出'] += o['trips']; C['正しい'] += int(o['ok']); C['無い注文で持ち帰り'] += int(o['absent'])
                if o['r'] == '無い': C['無い（本当に無い）'] += int(o['absent'])
                C['本当に無い注文'] += int(o.get('absent', False))
        print(f"{name}：持ち帰り {C['持ち帰り']}（正しい {C['正しい']}、うち X が無い注文で持ち帰った {C['無い注文で持ち帰り']}、外出の合計 {C['持ち帰りの外出']}）"
              f"  無いと返した {C['無い']}（本当に無い {C['無い（本当に無い）']}）  戻れない {C['戻れない']}  打ち切り {C['打ち切り']}  "
              f"外出の合計 {C['外出']}  歩数の合計 {C['歩数']}", flush=True)
