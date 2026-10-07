# 検定③改b：層1（40 番地）は 3 つの出口で共有、層2（10 番地）だけ出口 A/B/C ごとに写す。ほかは③改と同じ。
# 基準（走らせる前に決めたもの）
#   合格  通票あり − 通票なし ＞ 2 ポイント、かつ点の数だけの読みを上回る
#   負の対照  ラベルを混ぜる → 33% 前後
#   NG  差が 2 ポイント以内、または点の数だけの読み以下
import sys
argv = sys.argv[:]
import b13_amoeba10_warp3 as W3
A = W3.A; L = 3
_warp_floor = W3._warp_floor

def make_floor():
    fl = _warp_floor()
    kinds = fl['kind']
    shared = [v for v in range(fl['N']) if kinds[v] != '層2']            # 層0 と層1
    split = [v for v in range(fl['N']) if kinds[v] == '層2']
    sid = {v: i for i, v in enumerate(shared)}
    S = len(shared); cp = len(split); pos = {v: i for i, v in enumerate(split)}
    m = lambda c, v: sid[v] if v in sid else S + c * cp + pos[v]
    N = S + L * cp
    adj = [set() for _ in range(N)]; own = [None] * N
    for c in range(L):
        for v in range(fl['N']):
            for w in fl['adj'][v]:
                adj[m(c, v)].add(m(c, w))
        for v in split: own[m(c, v)] = c
    twin = {m(c, v): [m(k, v) for k in range(L)] for c in range(L) for v in split}
    xy = [fl['xy'][v] for v in shared] + [fl['xy'][split[i % cp]] for i in range(L * cp)]
    kind = [kinds[v] for v in shared] + ['層2' + 'ABC'[i // cp] for i in range(L * cp)]
    pix = [sid[p] for p in fl['pix']]; ex = [sid[e] for e in fl['ex']]
    fl.update(N=N, adj=[sorted(s) for s in adj], own=own, twin=twin, xy=xy, kind=kind, pix=pix, ex=ex,
              base=sum(1 for k in kind if not k.startswith('層')), cp=cp)
    return fl
A.make_floor = make_floor
W3.make_floor = make_floor

if __name__ == "__main__":
    sys.argv = argv
    src = open(W3.__file__).read().split('if __name__ == "__main__":')[1]
    src = src.replace('argv = sys.argv[:]', '').replace('globals()["STAT"]', 'W3.STAT')
    g = dict(W3.__dict__); g.update(argv=argv, A=A, make_floor=make_floor, W3=W3)
    import textwrap; exec(textwrap.dedent(src), g)
