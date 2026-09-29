"""通票（タブレット）方式のアメーバ
匂いの濃さは max(0, t-d)。登り先は「餌までの距離が一つ少ない隣」で、刻によらず決まる。
単線区間＝登り先が一つしかない番地が続く並び。区間は、登り先が二つ以上ある番地
（分岐＝駅）か自分の餌で終わる。
規則
 - 各番地に札を一枚。持ち主は0か1体。自分の居る番地の札は常に自分が持つ。
 - 駅（または区間の外）から区間に入るときは、区間の全番地の札が空いているときだけ
   まとめて受け取って入る。空いていなければ別の登り先の区間を試し、全部だめなら駅で待つ。
 - 区間の中では、札を持っている先へ一歩ずつ進み、抜けた番地の札を返す。
判定は整数の比較だけ。
"""
import random, json
from amoeba_multi import build_carrier, bfs
from amoeba_multi_v2 import make_map, summarize, KEYS

def run(pos, nbr, K, per, seed, T=400, loop=False):
    rng, blocked, free = make_map(pos, nbr, seed)
    pick = rng.sample(free, K + K * per)
    foods, starts = pick[:K], pick[K:]
    dist = [bfs(nbr, f, blocked) for f in foods]
    upc = {}
    def up(k, v):
        key = (k, v)
        if key not in upc:
            upc[key] = [w for w in nbr[v] if w not in blocked and dist[k].get(w, 10**9) == dist[k][v] - 1]
        return upc[key]
    am = [dict(k=i // per, at=s, steps=0, wsmell=0, wtraffic=0, done=None, last=0,
               d0=dist[i // per][s], ahead=[]) for i, s in enumerate(starts)]
    hold = {a['at']: i for i, a in enumerate(am)}
    for t in range(1, T + 1):
        live = [a for a in am if a['done'] is None]
        if not live or t - max(a['last'] for a in am) > 500: break
        order = list(range(len(am))); rng.shuffle(order)
        for i in order:
            a = am[i]
            if a['done'] is not None: continue
            u, k = a['at'], a['k']
            if t <= dist[k][u]:           # 匂いがまだ届いていない
                a['wsmell'] += 1; continue
            if not a['ahead']:
                ch = up(k, u)[:]; rng.shuffle(ch)
                for c in ch:
                    seq = [c]; x = c
                    while x != foods[k] and len(up(k, x)) == 1:
                        x = up(k, x)[0]; seq.append(x)
                    if all(hold.get(v, i) == i for v in seq):
                        for v in seq: hold[v] = i
                        a['ahead'] = seq; break
                if not a['ahead']:
                    a['wtraffic'] += 1
                    if loop:
                        # 駅の二本目の線：自分の番地を通りたい相手が駅で待っていて、
                        # 相手の方が餌に近い（同じなら番号の小さい方が先）なら、横の空き番地へ一歩よける
                        me = (dist[k][u], i)
                        yield_me = False
                        for j, b in enumerate(am):
                            if j == i or b['done'] is not None or b['ahead']: continue
                            if t <= dist[b['k']][b['at']]: continue
                            if me <= (dist[b['k']][b['at']], j): continue
                            for c in up(b['k'], b['at']):
                                seq = [c]; x = c
                                while x != foods[b['k']] and len(up(b['k'], x)) == 1:
                                    x = up(b['k'], x)[0]; seq.append(x)
                                if u in seq: yield_me = True; break
                            if yield_me: break
                        if yield_me:
                            side = [w for w in nbr[u] if w not in blocked and w not in hold]
                            if side:
                                w = rng.choice(side)
                                del hold[u]; hold[w] = i; a['at'] = w; a['steps'] += 1; a['last'] = t; a.setdefault('log', []).append((t, a['at']))
                    continue
            mv = a['ahead'].pop(0)
            del hold[u]; a['at'] = mv; a['steps'] += 1; a['last'] = t; a.setdefault('log', []).append((t, a['at']))
            if mv == foods[k]:
                a['done'] = t
                for v in [mv] + a['ahead']:
                    if hold.get(v) == i: del hold[v]
                a['ahead'] = []
    for a in am:
        a['end'] = 'own' if a['done'] is not None else ('wrong' if a['at'] in foods else 'none')
    return am

if __name__ == '__main__':
    import amoeba_siding as S
    pos, nbr = build_carrier(7)
    res = {}
    for name, fn in (('待避線', lambda K, s: S.run(pos, nbr, K, 2, s, True)),
                     ('通票', lambda K, s: run(pos, nbr, K, 2, s)),
                     ('通票＋駅の二本目', lambda K, s: run(pos, nbr, K, 2, s, loop=True))):
        for K in (2, 4, 8, 16, 32):
            tot = dict.fromkeys(KEYS, 0)
            for seed in range(20):
                for key, v in summarize(fn(K, seed)).items(): tot[key] += v
            res[f'{name} K={K}'] = tot
            print(f'{name} K={K}', tot, flush=True)
    json.dump(res, open('amoeba_tablet_result.json', 'w'), ensure_ascii=False, indent=1)
