"""避難：全員が同じ「出口の匂い」の層を読む。出口（搬送用エレベーター）は複数でもよい。
濃さ＝経過刻−最寄りの出口までの距離。濃い隣へ一歩、同じ番地に二体は入れない。
出口に着いたら運び出される。一つの出口が運び出せるのは一刻に一体まで。
一つの層を全員が下るので向かい合いは起きず、通票は要らない。"""
import cmath, math, random
from collections import deque
from amoeba_multi import build_carrier
from amoeba_multi_v2 import make_map

def msbfs(nbr, srcs, blocked):
    d = {s: 0 for s in srcs}; q = deque(srcs)
    while q:
        u = q.popleft()
        for v in nbr[u]:
            if v not in d and v not in blocked:
                d[v] = d[u] + 1; q.append(v)
    return d

def pick_exits(pos, free, n_rim, center=True):
    fs = set(free); ex = []
    targets = ([0j] if center else []) + [cmath.rect(0.78, math.pi/2 + 2*math.pi*j/5) for j in range(n_rim)]
    for z in targets:
        ex.append(min(fs, key=lambda v: abs(pos[v] - z)))
    return ex

def run(pos, nbr, seed, n_am, exits_rim, center=True, T=3000):
    rng, blocked, free = make_map(pos, nbr, seed)
    exits = pick_exits(pos, free, exits_rim, center)
    rest = [v for v in free if v not in exits]
    starts = rng.sample(rest, n_am)
    d = msbfs(nbr, exits, blocked)
    am = [dict(k=i // 2, at=s, d0=d[s], done=None, log=[], steps=0) for i, s in enumerate(starts)]
    occ = {a['at'] for a in am}; exset = set(exits)
    for t in range(1, T + 1):
        live = [a for a in am if a['done'] is None]
        if not live: break
        order = list(range(len(am))); rng.shuffle(order)
        used = set()                                    # この刻に運び出しを使った出口
        for i in order:
            a = am[i]
            if a['done'] is not None: continue
            u = a['at']
            if t <= d[u]: continue                      # 出口の匂いがまだ届いていない
            up = [v for v in nbr[u] if v in d and d[v] == d[u] - 1 and v not in occ and v not in used]
            if not up: continue
            v = rng.choice(up)
            occ.discard(u); a['at'] = v; a['steps'] += 1; a['log'].append((t, v))
            if v in exset: a['done'] = t; used.add(v)  # 一つの出口は一刻に一体を運び出す
            else: occ.add(v)
    for a in am: a['end'] = 'own' if a['done'] else 'none'
    return exits, starts, am

if __name__ == '__main__':
    pos, nbr = build_carrier(7)
    for n_am, rim, c in ((192, 0, True), (192, 5, True), (512, 0, True), (512, 5, True), (1024, 5, True)):
        res = []
        for seed in range(20):
            ex, st, am = run(pos, nbr, seed, n_am, rim, c)
            out = sum(1 for a in am if a['done'])
            last = max(a['done'] or 0 for a in am)
            lb = max(a['d0'] for a in am)
            res.append((out == n_am, last, lb))
        print(f'{n_am}体 出口{len(ex)}: 全員出た地図 {sum(r[0] for r in res)}/20  全員出るまで 最大{max(r[1] for r in res)} 刻'
              f'（渋滞なしの下限の最大 {max(r[2] for r in res)} 刻の2倍: 匂いが届く分を含む）  平均 {sum(r[1] for r in res)//20}')
