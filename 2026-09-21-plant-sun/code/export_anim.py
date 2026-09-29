import json
from amoeba_multi import build_carrier
from amoeba_multi_v2 import make_map
import amoeba_multi_v2 as V2, amoeba_siding as SD, amoeba_tablet as TB
pos, nbr = build_carrier(7)
inside = [v for v in range(len(pos)) if abs(pos[v]) < 0.85]
idx = {v: i for i, v in enumerate(inside)}
verts = [[round(pos[v].real * 1000), round(pos[v].imag * 1000)] for v in inside]
edges = sorted({tuple(sorted((idx[u], idx[w]))) for u in inside for w in nbr[u] if w in idx})
out = dict(verts=verts, edges=[list(e) for e in edges], maps={}, scenes=[])

def scene(key, title, note, seed, K, fn, T):
    rng, blocked, free = make_map(pos, nbr, seed)
    pick = rng.sample(free, K + K * 2)
    foods, starts = pick[:K], pick[K:]
    res = fn()
    am = res[1] if isinstance(res, tuple) else res
    mk = str(seed)
    if mk not in out['maps']:
        out['maps'][mk] = sorted(idx[v] for v in blocked if v in idx)
    A = []
    for i, a in enumerate(am):
        log = a.get('log', [])
        flat = []
        for t, v in log: flat += [t, idx[v]]
        A.append(dict(k=a['k'], s=idx[starts[i]], m=flat, d=a['done'] or 0, e=a['end']))
    tend = max([a['d'] for a in A] + [max(a['m'][-2::-2][:1] or [0]) for a in A]) + 20
    out['scenes'].append(dict(key=key, title=title, note=note, map=mk, K=K,
                              foods=[idx[f] for f in foods], am=A, T=min(tend, T)))
    own = sum(1 for a in A if a['e'] == 'own')
    print(key, len(A), own, tend)

S = 3
scene('mixed', '匂いを混ぜる', '全部の匂いを足した濃さを読む（負の対照）。多くが匂いの山の間で止まる。',
      S, 8, lambda: V2.run(pos, nbr, 8, 2, 'mixed', S, T=300), 300)
scene('plain', '自分の匂いだけ（待避なし）', '混じらずに追うが、細い所で向かい合うと固まる。',
      S, 8, lambda: SD.run(pos, nbr, 8, 2, S, False, T=700), 700)
scene('siding', '待避線', '向かい合ったら、餌から遠い方が横へ一歩よける。',
      S, 8, lambda: SD.run(pos, nbr, 8, 2, S, True, T=700), 700)
scene('tablet', '通票＋駅の二本目', '区間の札を全部受け取れたときだけ区間に入る。駅では横へよけて通す。',
      S, 8, lambda: TB.run(pos, nbr, 8, 2, S, T=700, loop=True), 700)
scene('limit', '上限：192体', '使える番地約1990に192体。20枚の地図すべてで全員が着く最大の数。',
      0, 96, lambda: TB.run(pos, nbr, 96, 2, 0, T=3000, loop=True), 3000)
scene('over', '上限超え：256体', '20枚中この1枚で24体が固まる。赤い輪が固まった個体。',
      1, 128, lambda: TB.run(pos, nbr, 128, 2, 1, T=3000, loop=True), 3000)
s = json.dumps(out, separators=(',', ':'), ensure_ascii=False)
open('anim_data.json', 'w').write(s); print(len(s))

# ---- 避難の場面 ----
import evac as EV
def evac_scene(key, title, note, seed, n_am, rim):
    exits, starts, am = EV.run(pos, nbr, seed, n_am, rim, True)
    rng, blocked, free = make_map(pos, nbr, seed)
    mk = str(seed)
    if mk not in out['maps']:
        out['maps'][mk] = sorted(idx[v] for v in blocked if v in idx)
    A = []
    for i, a in enumerate(am):
        flat = []
        for t, v in a['log']: flat += [t, idx[v]]
        A.append(dict(k=a['k'], s=idx[starts[i]], m=flat, d=a['done'] or 0, e=a['end']))
    T = max(a['d'] for a in A) + 20
    out['scenes'].append(dict(key=key, title=title, note=note, map=mk, K=len(exits), foods=[],
                              exits=[idx[e] for e in exits], am=A, T=T, evac=1))
    print(key, len(A), sum(1 for a in A if a['e'] == 'own'), T)

evac_scene('ev1', '避難：出口1つ', '全員が「出口の匂い」の層に切り替える。中心の搬送用エレベーター1基が一刻に一体ずつ運び出す。', 0, 192, 0)
evac_scene('ev6', '避難：出口6つ', '中心と外周5か所のエレベーター。同じ192体が約4分の1の刻で出る。', 0, 192, 5)
evac_scene('ev6b', '避難：1024体', '番地の約半分を埋めた混み具合でも、出口6つで全員が出る。', 0, 1024, 5)
s = json.dumps(out, separators=(',', ':'), ensure_ascii=False)
open('anim_data.json', 'w').write(s); print(len(s))
