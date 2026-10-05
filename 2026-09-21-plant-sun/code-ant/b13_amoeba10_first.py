# b13_amoeba10.py — 扇10枚の担体（gen10.py の回転対称図形）の上で、譲り合うアメーバに棒の向き 3 通りを読ませる
#
# 床：geo10.json の五角形のうち、中心から半径 20 以内（主の扇 5 枚・ずらしの扇 5 枚が全部入る）。
#   繋がり：辺の距離 φ で接する五角形どうし ＋ 扇をまたぐ φ² の繋がり（共有する扇のない組。CS1/AD3 と同じ）。
# 網膜：中心に近い五角形。棒（向き 0°・60°・120°）の上にある五角形が点く。
# 出口：床の縁で、向き 90°・210°・330° にいちばん近い五角形を一つずつ（答え 3 つ＝出口 3 つ・匂いの層 3 枚）。
# 規則は b13_amoeba_ai.py と同じ（通票・駅の二本目・待機所・乗り換え）。乗り換えは層を一つ先へ回す（0→1→2→0）。
# 学ぶもの：網膜の番地ごとの最初の層（0/1/2）、床の番地ごとの乗り換え（0/1）。
# 読み：最初に着いた出口（刻で数える。同じ刻の中に順番は付けない）。同じ刻に別の出口へ着いたら外れ、どこにも着かなければ外れ。
#
# 画像を作るところ（棒と五角形の中心の比較）だけ浮動小数。題材のデータで、装置の外。
#
# 基準（走らせる前に決めたもの）
#   合格  通票ありが、通票なし（重なってすり抜ける＝多数決）と点の数だけの読みの両方を上回り、33% を大きく越える
#   負の対照  ラベルを混ぜて学ばせ、本当のラベルで試験 → 33% 前後
#   NG  通票ありと通票なしの差が 2 ポイント以内、または点の数だけの読み以下
import json, math, random, sys, collections
from collections import Counter

GEO = "/home/claude/Paper/2026-09-21-plant-sun/code/geo10.json"
PHI = (1 + 5 ** 0.5) / 2

def make_floor(rho=20.0, ret=7.5):
    d = json.load(open(GEO)); R = d['R']; P = d['P']; cx, cy = d['center']
    cw = collections.defaultdict(set)
    for v in R:
        for p in v[4]: cw[p].add(v[5])
    S = [i for i in range(len(P)) if math.hypot(P[i][0] - cx, P[i][1] - cy) < rho]
    loc = {p: n for n, p in enumerate(S)}
    adj = [set() for _ in S]
    for a in S:
        for b in S:
            if b <= a: continue
            dd = math.hypot(P[a][0] - P[b][0], P[a][1] - P[b][1])
            if abs(dd - PHI) < 0.01 or (abs(dd - PHI * PHI) < 0.01 and not (cw[a] & cw[b])):
                adj[loc[a]].add(loc[b]); adj[loc[b]].add(loc[a])
    adj = [sorted(s) for s in adj]
    xy = [(P[p][0] - cx, P[p][1] - cy) for p in S]
    pix = [n for n in range(len(S)) if math.hypot(*xy[n]) < ret]
    rim = [n for n in range(len(S)) if math.hypot(*xy[n]) > rho - 3]
    ex = []
    for ang in (90, 210, 330):
        ux, uy = math.cos(math.radians(ang)), math.sin(math.radians(ang))
        ex.append(max(rim, key=lambda n: (xy[n][0] * ux + xy[n][1] * uy) - 0.01 * math.hypot(*xy[n])))
    fan = [sorted(cw[p]) for p in S]
    return dict(N=len(S), adj=adj, pix=pix, ex=ex, xy=xy, fan=fan)

def scent(fl, e):
    N = fl['N']; D = [0] * N; prev = D
    K = max(len(a) for a in fl['adj']) + 1
    T = 2 * N
    # 段数（最短）を先に数え、着くまで回す（第19章の式そのもの）
    dist = {e: 0}; fr = [e]
    while fr:
        nf = []
        for u in fr:
            for w in fl['adj'][u]:
                if w not in dist: dist[w] = dist[u] + 1; nf.append(w)
        fr = nf
    T = max(dist.values()) + 2
    for t in range(T + 1):
        nd = [sum(D[j] for j in fl['adj'][i]) for i in range(N)]
        nd[e] += K ** t
        prev, D = D, nd
    return [a + b for a, b in zip(D, prev)]

def simulate(fl, DS, img, start, swap, tokens=True, T=120, waits=None):
    L = len(DS)
    rob = [[fl['pix'][k], start[k]] for k in range(len(img)) if img[k]]
    arrived = [0] * L; exits = fl['ex']; exset = set(exits); first = None
    for t in range(T):
        if not rob: break
        occ = set(r[0] for r in rob) if tokens else set()
        claimed = set(); nxt = []
        for pos, lay in rob:
            D = DS[lay]
            nb = sorted(fl['adj'][pos], key=lambda j: -D[j])
            moved = False
            for c in nb[:2]:
                if D[c] <= D[pos]: break
                if tokens and (c in claimed or (c in occ and c not in exset)): continue
                claimed.add(c); pos2 = c; moved = True; break
            if not moved:
                pos2 = pos
                if waits is not None: waits.add(pos)
                if tokens and swap[pos]: lay = (lay + 1) % L
            if pos2 in exset:
                arrived[exits.index(pos2)] += 1
            else:
                nxt.append([pos2, lay])
        if first is None:
            got = {j for j in range(L) if arrived[j] > 0}
            if got: first = got.pop() if len(got) == 1 else -1   # 同じ刻に別の出口へ着いたら -1（NG）
        rob = nxt
    return arrived, first

def score(fl, DS, X, y, start, swap, tokens, waits=None):
    hit = 0; marg = 0
    for img, yy in zip(X, y):
        a, f = simulate(fl, DS, img, start, swap, tokens, waits=waits)
        other = max(a[j] for j in range(len(a)) if j != yy)
        hit += int(f == yy); marg += a[yy] - other     # 当たり＝最初に着いた出口（同時・未到着は外れ）。同点の比べ方は前のまま
    return hit, marg

def learn(fl, DS, X, y, tokens, rng, passes=8):
    P = len(fl['pix']); N = fl['N']; L = len(DS)
    start = [rng.randrange(L) for _ in range(P)]; swap = [0] * N
    W = set(); cur = score(fl, DS, X, y, start, swap, tokens, W)
    for _ in range(passes):
        moved = False
        # 待ちの起きない番地の乗り換え札は結果を変えないので、待ちが起きた番地だけ試す
        moves = [("s", k, v) for k in range(P) for v in range(L)] + ([("w", i, 0) for i in sorted(W)] if tokens else [])
        rng.shuffle(moves)
        for kind, k, v in moves:
            s2, w2 = start, swap
            if kind == "s":
                if start[k] == v: continue
                s2 = start.copy(); s2[k] = v
            else:
                w2 = swap.copy(); w2[k] ^= 1
            W2 = set(); sc = score(fl, DS, X, y, s2, w2, tokens, W2)
            if sc > cur: start, swap, cur, moved = s2, w2, sc, True; W |= W2
        if not moved: break
    return start, swap

ANG = (0, 60, 120)
def draw_bar(fl, cls, cx, cy, length, width):
    c, s = math.cos(math.radians(ANG[cls])), math.sin(math.radians(ANG[cls]))
    out = []
    for n in fl['pix']:
        x, y = fl['xy'][n]; dx, dy = x - cx, y - cy
        u = dx * c + dy * s; v = -dx * s + dy * c
        out.append(1 if abs(u) <= length / 2 and abs(v) <= width / 2 else 0)
    return out

def make_data(fl, n, rng, length=7.0, width=1.8, lo=2):
    r = max(math.hypot(*fl['xy'][k]) for k in fl['pix'])
    X, y = [], []
    while len(X) < n:
        cls = rng.randrange(3)
        cx, cy = rng.uniform(-r, r), rng.uniform(-r, r)
        img = draw_bar(fl, cls, cx, cy, length, width)
        if sum(img) < lo: continue
        X.append(img); y.append(cls)
    return X, y

def naive(Xtr, ytr, Xte, yte):
    tab = collections.defaultdict(Counter)
    for img, yy in zip(Xtr, ytr): tab[sum(img)][yy] += 1
    glob = Counter(ytr).most_common(1)[0][0]
    return sum(int((tab[sum(img)].most_common(1)[0][0] if tab[sum(img)] else glob) == yy) for img, yy in zip(Xte, yte))

if __name__ == "__main__":
    seeds = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    NTR = int(sys.argv[2]) if len(sys.argv) > 2 else 120
    NTE = 400
    fl = make_floor()
    DS = [scent(fl, e) for e in fl['ex']]
    print(f"床：五角形 {fl['N']}、網膜 {len(fl['pix'])}、出口 {fl['ex']}（扇 {[fl['fan'][e] for e in fl['ex']]}）", flush=True)
    tot = Counter()
    for s in range(seeds):
        rng = random.Random(8000 + s)
        Xtr, ytr = make_data(fl, NTR, rng); Xte, yte = make_data(fl, NTE, rng)
        res = []
        for tok in (True, False):
            st, sw = learn(fl, DS, Xtr, ytr, tok, rng)
            res.append((score(fl, DS, Xtr, ytr, st, sw, tok)[0], score(fl, DS, Xte, yte, st, sw, tok)[0], sum(sw)))
        ysh = list(ytr); rng.shuffle(ysh)
        stc, swc = learn(fl, DS, Xtr, ysh, True, rng)
        ctl = score(fl, DS, Xte, yte, stc, swc, True)[0]
        nv = naive(Xtr, ytr, Xte, yte)
        lit = round(sum(map(sum, Xte)) / NTE, 1)
        tot['tok'] += res[0][1]; tot['free'] += res[1][1]; tot['nv'] += nv; tot['ctl'] += ctl
        print(f"  {s}：通票あり 学習 {res[0][0]}/{NTR} 試験 {res[0][1]}/{NTE}（乗り換え {res[0][2]}） ｜ 通票なし 学習 {res[1][0]}/{NTR} 試験 {res[1][1]}/{NTE} ｜ 点の数 {nv} ｜ 対照 {ctl} ｜ 点の平均 {lit}", flush=True)
    n = NTE * seeds
    print(f"合計（試験 {n}）：通票あり {tot['tok']}  通票なし {tot['free']}  点の数 {tot['nv']}  対照 {tot['ctl']}")
