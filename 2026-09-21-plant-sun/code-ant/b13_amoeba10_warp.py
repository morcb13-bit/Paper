# 検定③：PVP のワープ。②の床（層0）に、同じ図形を中心の五芒星のまわりで φ^(3k) 倍した写し（層1・層2）を重ねる。
# 担体は写しに必要なだけ使う（全 7200 枚から取る）。層の中の繋がりは層0と同じ規則（辺の距離 φ、扇をまたぐ φ²）。
# ワープ：層 k の五角形 q と、それを含む層 k+1 の五角形（中心までのノルムが最小。Z[φ] の整数比較）を繋ぐ。
# 規則・読み（最初に着いた出口）・種・枚数・判定は②と同じ。
# 基準（走らせる前に決めたもの）
#   合格  通票あり − 通票なし ＞ 2 ポイント、かつ点の数だけの読みを上回る
#   負の対照  ラベルを混ぜる → 33% 前後
#   NG  差が 2 ポイント以内、または点の数だけの読み以下
import sys, math, pickle, random
from collections import Counter
sys.path.insert(0, '/home/claude/Paper/2026-09-21-plant-sun/code')
import b13_chain_units as U
import b13_amoeba10_gap as G
A = G.A
cells, z0 = pickle.load(open('/home/claude/icolearn/cells10.pkl', 'rb'))
cx, cy = U.xy(z0)
LAYERS = int(sys.argv[3]) if len(sys.argv) > 3 else 2

def zpow_phi(z, k):
    for _ in range(k): z = U.zmul(z, U.PHI)
    return z
def nlt(a, b): return U.phi_lt(U.norm2(a), U.norm2(b))   # |a| < |b|（整数で）

def layer_cells(k, rho):
    """層 k：担体の五角形 q を z0 + φ^(3k)(q − z0) に写す。床の半径 rho を覆う分だけ"""
    s = 1.618034 ** (3 * k)
    qs = [q for q in cells if math.hypot(U.xy(q)[0] - cx, U.xy(q)[1] - cy) * s < rho + 2.5 * s]
    pos = {q: U.zadd(z0, zpow_phi(U.zsub(q, z0), 3 * k)) for q in qs}
    return qs, pos

_base = A.make_floor
def make_floor():
    fl = _base()
    adj = [set(a) for a in fl['adj']]; xy = list(fl['xy']); kind = list(fl['kind']); fan = list(fl['fan'])
    # 層0の五角形の厳密な位置
    key = lambda x, y: (round(x, 2), round(y, 2))
    loc = {key(*fl['xy'][n]): n for n in range(len(fl['xy']))}
    prev = {}
    for q in cells:
        k0 = key(U.xy(q)[0] - cx, U.xy(q)[1] - cy)
        if k0 in loc and kind[loc[k0]] == '五角形': prev[loc[k0]] = q          # 番地 → 厳密な位置
    PHI = 1.618034
    for k in range(1, LAYERS + 1):
        qs, pos = layer_cells(k, 20.0)
        base = len(adj); idx = {q: base + i for i, q in enumerate(qs)}
        for q in qs:
            adj.append(set()); x, y = U.xy(pos[q]); xy.append((x - cx, y - cy)); kind.append(f'層{k}'); fan.append([])
        cw = {}
        for a in qs:
            for b in qs:
                if b <= a: continue
                dd = math.hypot(U.xy(a)[0] - U.xy(b)[0], U.xy(a)[1] - U.xy(b)[1])
                if abs(dd - PHI) < 0.01 or abs(dd - PHI * PHI) < 0.01:
                    adj[idx[a]].add(idx[b]); adj[idx[b]].add(idx[a])
        # ワープ：下の層の各五角形 → それを含む（中心が最も近い）この層の五角形
        for n, p in prev.items():
            best = None
            for q in qs:
                d = U.zsub(p, pos[q])
                if best is None or nlt(d, best[1]): best = (q, d)
            adj[n].add(idx[best[0]]); adj[idx[best[0]]].add(n)
        prev = {idx[q]: pos[q] for q in qs}
    fl.update(N=len(adj), adj=[sorted(s) for s in adj], xy=xy, kind=kind, fan=fan)
    return fl
A.make_floor = make_floor

if __name__ == "__main__":
    seeds = int(sys.argv[1]); NTR = int(sys.argv[2]); NTE = 400
    fl = A.make_floor(); DS = [A.scent(fl, e) for e in fl['ex']]
    used = Counter(fl['kind'][i] for i in range(fl['N']) if fl['adj'][i])
    print(f"床：番地 {fl['N']}（{dict(used)}）、隣の数の最大 {max(len(a) for a in fl['adj'])}", flush=True)
    tot = Counter()
    for s in range(seeds):
        rng = random.Random(8000 + s)
        Xtr, ytr = A.make_data(fl, NTR, rng); Xte, yte = A.make_data(fl, NTE, rng)
        res = []
        for tok in (True, False):
            st, sw = A.learn(fl, DS, Xtr, ytr, tok, rng)
            res.append(A.score(fl, DS, Xte, yte, st, sw, tok)[0])
        ysh = list(ytr); rng.shuffle(ysh)
        stc, swc = A.learn(fl, DS, Xtr, ysh, True, rng)
        ctl = A.score(fl, DS, Xte, yte, stc, swc, True)[0]
        nv = A.naive(Xtr, ytr, Xte, yte)
        tot['tok'] += res[0]; tot['free'] += res[1]; tot['nv'] += nv; tot['ctl'] += ctl
        print(f"  {s}：通票あり {res[0]}  通票なし {res[1]}  点の数 {nv}  対照 {ctl}", flush=True)
    print(f"学習用 {NTR}（試験 {NTE*seeds}）：通票あり {tot['tok']}  通票なし {tot['free']}  点の数 {tot['nv']}  対照 {tot['ctl']}")
