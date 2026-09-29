"""複線の手前で待機（待避線）を足したアメーバ
一本道（単線）は地図ではなく匂いが作る：自分の登り先が一つしかなく、そこに居る相手の
登り先も自分の番地しかない ＝ 向かい合い。このとき片方が横の空き番地（待避線）へ一歩よけ、
相手を通してから匂いに従って戻る。よける側は「自分の餌の匂いが薄い方（遠い方）」、
同じなら番号の大きい方。判定は整数の比較だけ。
"""
import random, json
from amoeba_multi import build_carrier, bfs
from amoeba_multi_v2 import make_map, summarize, KEYS

def run(pos, nbr, K, per, seed, siding=True, T=400):
    rng, blocked, free = make_map(pos, nbr, seed)
    pick = rng.sample(free, K + K * per)
    foods, starts = pick[:K], pick[K:]
    dist = [bfs(nbr, f, blocked) for f in foods]
    conc = lambda k, v, t: (t - dist[k][v]) if t > dist[k][v] else 0
    am = [dict(k=i // per, at=s, steps=0, wsmell=0, wtraffic=0, done=None, last=0,
               d0=dist[i // per][s], side=0) for i, s in enumerate(starts)]
    where = {a['at']: i for i, a in enumerate(am)}
    def ups(i, t):
        a = am[i]; u = a['at']; cu = conc(a['k'], u, t)
        return [v for v in nbr[u] if v not in blocked and conc(a['k'], v, t) > cu]
    for t in range(1, T + 1):
        live = [a for a in am if a['done'] is None]
        if not live or t - max(a['last'] for a in am) > 500: break
        order = list(range(len(am))); rng.shuffle(order)
        for i in order:
            a = am[i]
            if a['done'] is not None: continue
            u, k = a['at'], a['k']
            up_all = ups(i, t)
            up = [v for v in up_all if v not in where]
            mv = None
            if up:
                m = max(conc(k, v, t) for v in up)
                mv = rng.choice([v for v in up if conc(k, v, t) == m])
            elif not up_all:
                a['wsmell'] += 1
            else:
                a['wtraffic'] += 1
                if siding:
                    # 向かい合いか：登り先に居る相手の登り先が自分の番地しかない
                    for v in up_all:
                        j = where.get(v)
                        if j is None: continue
                        if ups(j, t) != [u]: continue
                        b = am[j]
                        mine, theirs = conc(k, u, t), conc(b['k'], v, t)
                        yield_me = mine < theirs or (mine == theirs and i > j)
                        if yield_me:
                            side = [w for w in nbr[u] if w not in blocked and w not in where and w != v]
                            if side:
                                mv = rng.choice(side); a['side'] += 1
                        break
            if mv is None: continue
            del where[u]; a['at'] = mv; a['steps'] += 1; a['last'] = t; a.setdefault('log', []).append((t, a['at']))
            if mv == foods[k]:
                a['done'] = t
            else:
                where[mv] = i
    for a in am:
        a['end'] = 'own' if a['done'] is not None else ('wrong' if a['at'] in foods else 'none')
    return am

if __name__ == '__main__':
    pos, nbr = build_carrier(7)
    res = {}
    for siding in (False, True):
        for K in (2, 4, 8, 16):
            tot = dict.fromkeys(KEYS + ('side',), 0)
            for seed in range(20):
                am = run(pos, nbr, K, 2, seed, siding)
                for key, v in summarize(am).items(): tot[key] += v
                tot['side'] += sum(a['side'] for a in am)
            name = f"{'待避あり' if siding else '待避なし'} K={K}"
            res[name] = tot
            print(name, tot, flush=True)
    json.dump(res, open('amoeba_siding_result.json', 'w'), ensure_ascii=False, indent=1)
