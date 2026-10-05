# 検定③改：PVP のワープの層（層1・層2）を出口 A/B/C ごとに一組ずつ写す。層0（②の床 281 番地）は全員で共有。
# 匂いの層 l のロボットは、層0と写し l だけを使う（匂い D_l は層0＋写し l の網で回す。ほかの写しには匂いがない）。
# 乗り換え：写し l の番地で待たされたら、同じ位置の写し l+1 の番地へ移って層を替える。層0では層だけ替える（③と同じ）。
# 規則・読み（最初に着いた出口）・種・枚数・判定は③と同じ。
# 基準（走らせる前に決めたもの）
#   合格  通票あり − 通票なし ＞ 2 ポイント、かつ点の数だけの読みを上回る
#   負の対照  ラベルを混ぜる → 33% 前後
#   NG  差が 2 ポイント以内、または点の数だけの読み以下
import sys
argv = sys.argv[:]; sys.argv = [sys.argv[0], '0', '0']
import random
from collections import Counter
import b13_amoeba10_warp as W
A = W.A
_warp_floor = A.make_floor
L = 3

def make_floor():
    fl = _warp_floor()
    base = sum(1 for k in fl['kind'] if not k.startswith('層'))       # 層0 の番地数（281）
    cp = fl['N'] - base                                                # 写し一組の番地数（50）
    m = lambda c, v: v if v < base else base + c * cp + (v - base)
    N = base + L * cp
    adj = [set() for _ in range(N)]; own = [None] * N
    for c in range(L):
        for v in range(fl['N']):
            for w in fl['adj'][v]:
                if v < base and w < base: adj[v].add(w)
                else: adj[m(c, v)].add(m(c, w))
        for v in range(base, fl['N']): own[m(c, v)] = c
    twin = [[m(c, v) for c in range(L)] if v >= base else None for v in range(fl['N'])]
    twin = {m(c, v): twin[v] for c in range(L) for v in range(base, fl['N'])}
    xy = fl['xy'][:base] + [fl['xy'][base + (i % cp)] for i in range(L * cp)]
    kind = fl['kind'][:base] + [fl['kind'][base + (i % cp)] + 'ABC'[i // cp] for i in range(L * cp)]
    fl.update(N=N, adj=[sorted(s) for s in adj], own=own, twin=twin, xy=xy, kind=kind, base=base, cp=cp)
    return fl

def scent_l(fl, l):
    sub = dict(fl); sub['adj'] = [[w for w in a if fl['own'][w] in (None, l)] if fl['own'][i] in (None, l) else []
                                  for i, a in enumerate(fl['adj'])]
    return A.scent(sub, fl['ex'][l])

STAT = None
def simulate(fl, DS, img, start, swap, tokens=True, T=120, waits=None):
    rob = [[fl['pix'][k], start[k]] for k in range(len(img)) if img[k]]
    arrived = [0] * L; exits = fl['ex']; exset = set(exits); first = None
    for t in range(T):
        if not rob: break
        occ = set(r[0] for r in rob) if tokens else set()
        claimed = set(); nxt = []
        for pos, lay in rob:
            D = DS[lay]; g = fl['kind'][pos][:2] if fl['kind'][pos].startswith('層') else '層0'
            if STAT is not None: STAT[g]['手'] += 1
            nb = sorted(fl['adj'][pos], key=lambda j: -D[j]); moved = False
            for c in nb[:2]:
                if D[c] <= D[pos]: break
                if tokens and (c in claimed or (c in occ and c not in exset)): continue
                claimed.add(c); pos2 = c; moved = True; break
            if not moved:
                pos2 = pos
                if STAT is not None: STAT[g]['待機'] += 1
                if waits is not None: waits.add(pos)
                if tokens and swap[pos]:
                    lay = (lay + 1) % L
                    if fl['own'][pos] is not None: pos2 = fl['twin'][pos][lay]
            if pos2 in exset: arrived[exits.index(pos2)] += 1
            else: nxt.append([pos2, lay])
        if first is None:
            got = {j for j in range(L) if arrived[j] > 0}
            if got: first = got.pop() if len(got) == 1 else -1
        rob = nxt
    return arrived, first

A.make_floor = make_floor; A.simulate = simulate

if __name__ == "__main__":
    seeds = int(argv[1]); NTR = int(argv[2]); NTE = 400
    fl = A.make_floor(); DS = [scent_l(fl, l) for l in range(L)]
    print(f"床：番地 {fl['N']}（層0 {fl['base']}、写し {fl['cp']}×3）", flush=True)
    # 一体だけのときの歩数
    for l in range(L):
        st = []
        for p in fl['pix']:
            pos = p; n = 0
            while pos != fl['ex'][l] and n < 200:
                nb = max(fl['adj'][pos], key=lambda j: DS[l][j])
                if DS[l][nb] <= DS[l][pos]: break
                pos = nb; n += 1
            st.append(n if pos == fl['ex'][l] else None)
        print(f"  一体だけ 出口{'ABC'[l]}：着かない {st.count(None)}  歩数 {min(x for x in st if x is not None)}〜{max(x for x in st if x is not None)}")
    tot = Counter(); STATS = {t: {g: Counter() for g in ('層0', '層1', '層2')} for t in (True, False)}
    for s in range(seeds):
        rng = random.Random(8000 + s)
        Xtr, ytr = A.make_data(fl, NTR, rng); Xte, yte = A.make_data(fl, NTE, rng)
        res = []
        for tok in (True, False):
            st, sw = A.learn(fl, DS, Xtr, ytr, tok, rng)
            globals()["STAT"] = STATS[tok]; res.append(A.score(fl, DS, Xte, yte, st, sw, tok)[0]); globals()["STAT"] = None
        ysh = list(ytr); rng.shuffle(ysh)
        stc, swc = A.learn(fl, DS, Xtr, ysh, True, rng)
        ctl = A.score(fl, DS, Xte, yte, stc, swc, True)[0]
        nv = A.naive(Xtr, ytr, Xte, yte)
        tot['tok'] += res[0]; tot['free'] += res[1]; tot['nv'] += nv; tot['ctl'] += ctl
        print(f"  {s}：通票あり {res[0]}  通票なし {res[1]}  点の数 {nv}  対照 {ctl}", flush=True)
    print(f"学習用 {NTR}（試験 {NTE*seeds}）：通票あり {tot['tok']}  通票なし {tot['free']}  点の数 {tot['nv']}  対照 {tot['ctl']}")
    for tok in (True, False):
        print(f"通票{'あり' if tok else 'なし'}：" + "  ".join(f"{g} 手 {c['手']} 待機 {c['待機']}（{100*c['待機']//max(1,c['手'])}%）" for g, c in STATS[tok].items()))
