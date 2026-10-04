# b13_ricci.py — ペンタゴン・リッチフローAIプロセッサ：伸び縮みするセル
#
# セル＝五角形12枚を捻りながら押し出す立体。捻り θ は 36° 刻みの番号 j（20 で一周＝表裏つき）。
#   j が偶数（θ が 72° の倍数）：正12面体に閉じる → 信号を受けない・出さない
#   j が奇数（θ が 36° ずれ）　：切頂20面体（フラーレン）に開く → 届いた指しを扇に着地させ、捻り j を足して隣へ渡す
# 各セルは周波数 ω を持ち、一刻ごとに j ← j + ω（ω = 1 なら一刻ごとに開閉する＝一回転で単振動5回の離散版）。
# 信号は、送り手がその刻に開いていて、受け手が次の刻に開いているときだけ届く。
# 学習の一手：セルの位相 p を ±1、表裏を入れ替える（+10）、周波数 ω を ±1。外れが減るときだけ採る。
#
# 基準（走らせる前に決めたもの）
#   R1  見せていない入力の当たりが 6 割以上（前回までの網は半分前後）
#   対照 札を混ぜると半分前後
#   R2  学んだあと、入力 i を裏返すと出力が変わる例があるか（効いている入力）を数える
#       予言：「同符号」は入力 1・2 だけが効き、「一つ目が正」は入力 1 だけが効く
#   NG  R1 で半分前後
import itertools, random, sys
import numpy as np
import b13_fannet as F
import b13_spinnet as S

def run(net, edges, p, w, X, T):
    N = net['N']; Ex = X.shape[0]
    src = np.array([a for a, b in edges], dtype=np.int64)
    dst = np.array([b for a, b in edges], dtype=np.int64)
    inp = np.array(net['inputs']); out = net['out']
    xin = np.where(X > 0, 0, np.where(X < 0, 10, -1))
    last = -np.ones(Ex, dtype=np.int64)
    Pm = -np.ones((Ex, N), dtype=np.int64)              # 各セルがこの刻に送り出す指し
    for t in range(T):
        j = (p + w * t) % 20                             # 捻り
        opn = (j % 2 == 1)
        # 入力のセル：開いていれば入力の指しに捻りを足して送る
        jin = j[inp]
        Pm[:, inp] = np.where((xin >= 0) & opn[inp][None, :], (xin + jin[None, :]) % 20, -1)
        # 届く：送り手が今開いていて送り出し、受け手が次の刻に開いている
        jn = (p + w * (t + 1)) % 20; opn_next = (jn % 2 == 1)
        C = np.zeros((Ex, N, 10), dtype=np.int64); Z = np.zeros((Ex, N), dtype=np.int64)
        ok = (Pm[:, src] >= 0) & opn_next[dst][None, :]
        e_idx, k_idx = np.nonzero(ok)
        pj = Pm[e_idx, src[k_idx]]
        np.add.at(C, (e_idx, dst[k_idx], pj % 10), 1)
        np.add.at(Z, (e_idx, dst[k_idx]), pj // 10)
        K = F.land(C)
        J = np.where(K >= 0, K + 10 * (Z % 2), -1)
        J = np.where(opn_next[None, :], J, -1)
        last = np.where(opn_next[out] & (J[:, out] >= 0), J[:, out], np.where(opn_next[out], -1, last))
        Pm = np.where(J >= 0, (J + jn[None, :]) % 20, -1)
    return last

CLOSE = {"on": False}

def tuned_phase(net, rng):
    """出力からの段数の偶奇で位相の偶奇をそろえる（隣どうしが交互に開く＝全部が共鳴して通る状態）"""
    adj = {}
    for a, b in net['nb']: adj.setdefault(a, []).append(b); adj.setdefault(b, []).append(a)
    d = {net['out']: 0}; fr = [net['out']]
    while fr:
        nf = []
        for u in fr:
            for v in adj.get(u, []):
                if v not in d: d[v] = d[u] + 1; nf.append(v)
        fr = nf
    return np.array([2 * rng.randrange(10) + (d.get(v, 0) % 2) for v in range(net['N'])], dtype=np.int64)

def learn(net, edges, X, y, T, rng, passes=12):
    N = net['N']
    p = tuned_phase(net, rng)
    w = np.ones(N, dtype=np.int64)
    cur = S.score(run(net, edges, p, w, X, T), y)
    for _ in range(passes):
        moved = False
        kinds = (("p", 1), ("p", -1), ("p", 10), ("w", 1), ("w", -1)) + ((("close", 0),) if CLOSE["on"] else ())
        moves = [(v, k, d) for v in range(N) for k, d in kinds]
        rng.shuffle(moves)
        for v, k, d in moves:
            p2, w2 = p, w
            if k == "p": p2 = p.copy(); p2[v] = (p2[v] + d) % 20
            elif k == "w": w2 = w.copy(); w2[v] = (w2[v] + d) % 20
            else:   # 閉じたまま止める（正12面体で振動をやめる）
                if w[v] == 0 and p[v] % 2 == 0: continue
                p2 = p.copy(); w2 = w.copy(); p2[v] -= p2[v] % 2; w2[v] = 0
            s2 = S.score(run(net, edges, p2, w2, X, T), y)
            if s2 > cur or (k == "close" and s2 == cur):
                p, w, moved = p2, w2, (moved or s2 > cur); cur = s2
        if not moved: break
    return p, w

def effective(net, edges, p, w, X, T):
    base = run(net, edges, p, w, X, T)
    eff = []
    for i in range(X.shape[1]):
        X2 = X.copy(); X2[:, i] = -X2[:, i]
        eff.append(int((run(net, edges, p, w, X2, T) != base).sum()))
    return eff

if __name__ == "__main__":
    net = F.make_floor(); m = len(net['inputs'])
    ALL0 = np.array([s for s in itertools.product((-1, 0, 1), repeat=m) if any(s)])
    NETS = {"A 一方向": (net['fwd'], net['L'] + 1), "C 隣＋φ²": (F.directed(net['nb'] + net['phi2']), net['L'] + 5)}
    seeds = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    NTR = int(sys.argv[2]) if len(sys.argv) > 2 else 100
    only = sys.argv[3] if len(sys.argv) > 3 else None
    for tname, (f, sel) in S.TASKS.items():
        if only and tname != only: continue
        A = sel(ALL0); y_all = f(A)
        print(f"題「{tname}」 入力 {len(A)} 通り（はい {int(y_all.sum())}）、学習用 {NTR}、{seeds} 回の合計", flush=True)
        for nname, (edges, T) in NETS.items():
            trh = teh = ctl = tot = 0; effs = []; ws = []
            for s in range(seeds):
                rng = random.Random(4000 + s); idx = list(range(len(A))); rng.shuffle(idx)
                tr, te = np.array(idx[:NTR]), np.array(idx[NTR:])
                p, w = learn(net, edges, A[tr], y_all[tr], T, rng)
                trh += S.score(run(net, edges, p, w, A[tr], T), y_all[tr])
                teh += S.score(run(net, edges, p, w, A[te], T), y_all[te])
                effs.append(effective(net, edges, p, w, A, T)); ws.append(int((w != 1).sum()))
                perm = list(range(len(A))); rng.shuffle(perm); ysh = y_all[perm]
                pc, wc = learn(net, edges, A[tr], ysh[tr], T, rng)
                ctl += S.score(run(net, edges, pc, wc, A[te], T), ysh[te]); tot += len(te)
            print(f"  {nname}：学習用 {trh}/{NTR*seeds}  試験 {teh}/{tot}  ｜ 対照 {ctl}/{tot}  ｜ "
                  f"入力ごとに裏返すと出力が変わった例の数 {effs}  ｜ 周波数を変えたセル {ws}", flush=True)
