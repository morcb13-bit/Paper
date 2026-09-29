"""混じらずに違う匂いを追うアメーバ（v2）― ペンローズ担体上の試験
v1 からの変更
 - 他の餌の番地は通り抜けられる。着いたと数えるのは自分の餌だけ。
   最後に他の餌の上で止まっていたら「取り違え」。
 - 待ちを二つに分ける：匂いがまだ届かず動けない待ち／隣が塞がっている待ち（渋滞）。
 - 匂いの層の濃さは「届いた回数」＝ max(0, t - d)。v1 の逐次シミュレーションと
   全番地で一致することを確認済み（不一致 0）なので、この式で数える。整数のみ。
"""
import random, json
from amoeba_multi import build_carrier, bfs

def make_map(pos, nbr, seed):
    rng = random.Random(seed)
    n = len(pos)
    inside = [v for v in range(n) if abs(pos[v]) < 0.85]
    inset = set(inside)
    blocked = set(v for v in range(n) if v not in inset)
    blocked |= set(rng.sample(inside, len(inside) * 15 // 100))   # 棚
    free = [v for v in inside if v not in blocked]
    best = {}
    for v in free:
        if v not in best:
            cc = bfs(nbr, v, blocked)
            if len(cc) > len(best): best = cc
    blocked |= set(v for v in inside if v not in best)
    return rng, blocked, list(best)

def run(pos, nbr, K, per, mode, seed, T=400):
    rng, blocked, free = make_map(pos, nbr, seed)
    pick = rng.sample(free, K + K * per)
    foods, starts = pick[:K], pick[K:]
    dist = [bfs(nbr, f, blocked) for f in foods]
    conc = lambda k, v, t: (t - dist[k][v]) if t > dist[k][v] else 0
    am = [dict(k=i // per, at=s, steps=0, wsmell=0, wtraffic=0, done=None,
               d0=dist[i // per][s], path=[s]) for i, s in enumerate(starts)]
    occ = {a['at'] for a in am}
    for t in range(1, T + 1):
        order = list(range(len(am))); rng.shuffle(order)
        for i in order:
            a = am[i]
            if a['done'] is not None: continue
            u, k = a['at'], a['k']
            if mode == 'own':
                val = lambda v: conc(k, v, t)
            elif mode == 'mixed':
                val = lambda v: sum(conc(j, v, t) for j in range(K))
            opn = [v for v in nbr[u] if v not in blocked]
            cand = [v for v in opn if v not in occ]
            if mode == 'blind':
                mv = rng.choice(cand) if cand else None
                if mv is None: a['wtraffic'] += 1
            else:
                cu = val(u)
                up_all = [v for v in opn if val(v) > cu]
                up = [v for v in up_all if v not in occ]
                if up:
                    m = max(val(v) for v in up)
                    mv = rng.choice([v for v in up if val(v) == m])
                else:
                    mv = None
                    if up_all: a['wtraffic'] += 1
                    else: a['wsmell'] += 1
            if mv is None: continue
            occ.discard(u); a['at'] = mv; a['steps'] += 1; a['path'].append(mv); a.setdefault('log', []).append((t, mv))
            if mv == foods[k]:
                a['done'] = t          # 自分の餌に着いたら抜ける
            else:
                occ.add(mv)
    for a in am:
        if a['done'] is not None: a['end'] = 'own'
        elif a['at'] in foods: a['end'] = 'wrong'
        else: a['end'] = 'none'
    return foods, am, blocked

KEYS = ('n', 'own', 'wrong', 'none', 'excess', 'wsmell', 'wtraffic', 'delay')

def summarize(am):
    s = dict.fromkeys(KEYS, 0)
    s['n'] = len(am)
    for a in am:
        s[a['end']] += 1
        if a['end'] == 'own':
            s['excess'] += a['steps'] - a['d0']      # 最短より余分に歩いた歩数
            s['wsmell'] += a['wsmell']
            s['wtraffic'] += a['wtraffic']
            s['delay'] += a['done'] - a['d0']         # 着いた刻 − 最短距離
    return s

if __name__ == '__main__':
    pos, nbr = build_carrier(7)
    res = {}
    print('担体の頂点', len(pos))
    for K in (2, 4, 8):
        for mode in ('own', 'mixed', 'blind'):
            tot = dict.fromkeys(KEYS, 0)
            for seed in range(20):
                _, am, _ = run(pos, nbr, K, 2, mode, seed)
                for key, v in summarize(am).items(): tot[key] += v
            res[f'K{K}_{mode}'] = tot
            print(f'K={K} {mode:5s}', tot, flush=True)
    json.dump(res, open('amoeba_multi_v2_result.json', 'w'), ensure_ascii=False, indent=1)
