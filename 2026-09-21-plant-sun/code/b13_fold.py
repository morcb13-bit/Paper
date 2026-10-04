# b13_fold.py — 折りたたまれた次元：繋がりの強さを五芒星の内側の桁で持たせる
#
# 辺 e の強さ：w_e = Σ_{k<K} d_k φ^(−3k)、桁 d_k ∈ {−2,…,+2}（平衡5進）。一段内側へ潜るごとに φ³ 分の 1。
# 整数で持つため全体を φ^(3(K−1)) 倍して W_e = Σ d_k φ^(3(K−1−k)) ∈ Z[φ] とする（比べるだけなので倍率は効かない）。
# 届いた指しは「1 個」ではなく W_e 個ぶん、その扇に足される。着地は Z[φ] の内積の比較だけ（浮動小数なし）。
# 表裏は、強さ 0 でない届き方のうち裏の数の偶奇。
# 学習：粗い桁から。段 k を開くと、その段の桁を ±1 動かす手が加わる（位相・周波数・閉じる手はそのまま）。
#
# 基準（走らせる前に決めたもの）
#   F1  段を増やす（K = 0 → 1 → 2 → 3）ほど、見せていない棒の当たりが上がること。K=0 は今の網（約 69%）
#   F2  学習用の当たりも上がること（器が広がったか）
#   対照 札を混ぜると半分前後
#   NG  K を増やしても 69% から動かない
# 当たり数だけで手を選ぶと、桁を一つ動かしても当たり数が変わらず、桁は一つも動かなかった（K=1・K=3 とも）。
# そこで当たり数が同じときは、出力での着地の余白で比べる。
import random, sys
import numpy as np
import b13_fannet as F, b13_spinnet as S, b13_ricci as R
import b13_retina_ai as RA

CA = F.CA; CB = F.CB
MARGIN = {"on": True}   # 当たり数が同じなら、出力の着地の余白（一番と二番の差）が増える手を採る

def zmul(a1, b1, a2, b2):  # (a1 + b1 φ)(a2 + b2 φ)
    return a1 * a2 + b1 * b2, a1 * b2 + a2 * b1 + b1 * b2

MCA = np.stack([np.roll(CA, k) for k in range(10)], 1)   # [j, k] = 2cos(36(j−k)) の整数部
MCB = np.stack([np.roll(CB, k) for k in range(10)], 1)   # φ の係数
def land_phi(Ca, Cb):
    """扇ごとの強さ（Z[φ]）→ 着地した扇（並んだら −1）。内積は整数の行列の掛け算（足し算の繰り返し）"""
    A = Ca @ MCA + Cb @ MCB
    B = Ca @ MCB + Cb @ (MCA + MCB)
    best = np.zeros(Ca.shape[:-1], dtype=np.int64)
    ba, bb = A[..., 0].copy(), B[..., 0].copy()
    sa = np.full_like(ba, -10**15); sb = np.zeros_like(bb)          # 二番目（始めはとても小さい）
    for k in range(1, 10):
        s1 = F.vsign(A[..., k] - ba, B[..., k] - bb)
        up = s1 > 0
        s2 = F.vsign(A[..., k] - sa, B[..., k] - sb)
        up2 = (~up) & (s2 > 0)
        sa = np.where(up, ba, np.where(up2, A[..., k], sa)); sb = np.where(up, bb, np.where(up2, B[..., k], sb))
        best = np.where(up, k, best)
        ba = np.where(up, A[..., k], ba); bb = np.where(up, B[..., k], bb)
    ma, mb = ba - sa, bb - sb                                       # 余白＝一番 − 二番（Z[φ]、0 なら並び）
    tie = F.vsign(ma, mb) <= 0
    LAST_MARGIN["a"], LAST_MARGIN["b"] = np.where(tie, 0, ma), np.where(tie, 0, mb)
    return np.where(tie, -1, best)

LAST_MARGIN = {"a": None, "b": None}

PHI3 = [(1, 0), (1, 2), (5, 8), (21, 34)]   # φ^0, φ^3 = 1+2φ, φ^6 = 5+8φ, φ^9 = 21+34φ

def weights(D, K):
    """D: (辺, 3) の桁 → (Wa, Wb)"""
    Wa = np.zeros(D.shape[0], dtype=np.int64); Wb = np.zeros(D.shape[0], dtype=np.int64)
    for k in range(K):
        a, b = PHI3[K - 1 - k]
        Wa += D[:, k] * a; Wb += D[:, k] * b
    return Wa, Wb

def run(net, edges, p, w, Wa, Wb, X, T):
    N = net['N']; Ex = X.shape[0]
    src = np.array([a for a, b in edges], dtype=np.int64)
    dst = np.array([b for a, b in edges], dtype=np.int64)
    nz = (Wa != 0) | (Wb != 0)
    inp = np.array(net['inputs']); out = net['out']
    xin = np.where(X > 0, 0, np.where(X < 0, 10, -1))
    last = -np.ones(Ex, dtype=np.int64)
    lma = np.zeros(Ex, dtype=np.int64); lmb = np.zeros(Ex, dtype=np.int64)
    Pm = -np.ones((Ex, N), dtype=np.int64)
    for t in range(T):
        j = (p + w * t) % 20; opn = (j % 2 == 1)
        Pm[:, inp] = np.where((xin >= 0) & opn[inp][None, :], (xin + j[inp][None, :]) % 20, -1)
        jn = (p + w * (t + 1)) % 20; opn_next = (jn % 2 == 1)
        Ca = np.zeros((Ex, N, 10), dtype=np.int64); Cb = np.zeros((Ex, N, 10), dtype=np.int64)
        Z = np.zeros((Ex, N), dtype=np.int64)
        ok = (Pm[:, src] >= 0) & opn_next[dst][None, :] & nz[None, :]
        e_idx, k_idx = np.nonzero(ok)
        pj = Pm[e_idx, src[k_idx]]
        np.add.at(Ca, (e_idx, dst[k_idx], pj % 10), Wa[k_idx])
        np.add.at(Cb, (e_idx, dst[k_idx], pj % 10), Wb[k_idx])
        np.add.at(Z, (e_idx, dst[k_idx]), pj // 10)
        K = land_phi(Ca, Cb)
        mA, mB = LAST_MARGIN["a"][:, out], LAST_MARGIN["b"][:, out]
        J = np.where(K >= 0, K + 10 * (Z % 2), -1)
        J = np.where(opn_next[None, :], J, -1)
        upd = opn_next[out]
        last = np.where(upd & (J[:, out] >= 0), J[:, out], np.where(upd, -1, last))
        if upd: lma, lmb = np.where(J[:, out] >= 0, mA, 0), np.where(J[:, out] >= 0, mB, 0)
        Pm = np.where(J >= 0, (J + jn[None, :]) % 20, -1)
    OUT_MARGIN["a"], OUT_MARGIN["b"] = lma, lmb
    return last

OUT_MARGIN = {"a": None, "b": None}
def mscore(out, y):
    """(当たり数, 余白の和)：当たった例の余白は足し、外れた例の余白は引く（Z[φ]）"""
    pred = ((out >= 0) & (out < 10)).astype(np.int64)
    ok = pred == y
    sg = np.where(ok, 1, -1)
    return int(ok.sum()), (int((sg * OUT_MARGIN["a"]).sum()), int((sg * OUT_MARGIN["b"]).sum()))

def learn(net, edges, X, y, T, rng, Kmax, passes=6):
    N = net['N']; E = len(edges)
    p = R.tuned_phase(net, rng); w = np.ones(N, dtype=np.int64)
    D = np.zeros((E, 3), dtype=np.int64); D[:, 0] = 1           # 外側の桁 1＝今の網
    K = max(Kmax, 1)
    def ev(p, w, D):
        Wa, Wb = weights(D, K)
        return mscore(run(net, edges, p, w, Wa, Wb, X, T), y) if MARGIN["on"] else (S.score(run(net, edges, p, w, Wa, Wb, X, T), y), (0, 0))
    cur = ev(p, w, D); hist = []
    for stage in range(0, Kmax + 1 if Kmax > 0 else 1):
        # stage 0：位相・周波数・閉じる手だけ／stage k≥1：桁 k−1 を動かす手を加える
        for _ in range(passes):
            moved = False
            moves = [("n", v, kd) for v in range(N) for kd in (("p", 1), ("p", -1), ("p", 10), ("w", 1), ("w", -1), ("close", 0))]
            if stage >= 1:
                moves += [("e", e, d) for e in range(E) for d in (1, -1)]
            rng.shuffle(moves)
            for kind, a, d in moves:
                p2, w2, D2 = p, w, D
                if kind == "n":
                    k, dd = d
                    if k == "p": p2 = p.copy(); p2[a] = (p2[a] + dd) % 20
                    elif k == "w": w2 = w.copy(); w2[a] = (w2[a] + dd) % 20
                    else:
                        if w[a] == 0 and p[a] % 2 == 0: continue
                        p2 = p.copy(); w2 = w.copy(); p2[a] -= p2[a] % 2; w2[a] = 0
                else:
                    nd = D[a, stage - 1] + d
                    if nd < -2 or nd > 2: continue
                    D2 = D.copy(); D2[a, stage - 1] = nd
                s2 = ev(p2, w2, D2)
                if F.better(s2, cur) or (kind == "n" and d[0] == "close" and s2 == cur):
                    p, w, D, moved = p2, w2, D2, (moved or F.better(s2, cur)); cur = s2
            if not moved: break
        hist.append(cur[0])
    return p, w, D, K, hist

if __name__ == "__main__":
    Kmax = int(sys.argv[1]); seed = int(sys.argv[2]); NTR = int(sys.argv[3]) if len(sys.argv) > 3 else 120
    NTE = 400
    net, cells, verts, d = RA.make_retina_net()
    T = net['L'] + 3; E = net['fwd']
    rng = random.Random(6000 + seed)
    Xtr, ytr = RA.make_data(cells, NTR, rng); Xte, yte = RA.make_data(cells, NTE, rng)
    p, w, D, K, hist = learn(net, E, Xtr, ytr, T, rng, Kmax)
    Wa, Wb = weights(D, K)
    a_tr = S.score(run(net, E, p, w, Wa, Wb, Xtr, T), ytr)
    a_te = S.score(run(net, E, p, w, Wa, Wb, Xte, T), yte)
    ysh = ytr.copy(); rng.shuffle(ysh)
    pc, wc, Dc, Kc, _ = learn(net, E, Xtr, ysh, T, rng, Kmax)
    Wac, Wbc = weights(Dc, Kc)
    yte_sh = np.array([rng.randrange(2) for _ in yte])
    ctl = S.score(run(net, E, pc, wc, Wac, Wbc, Xte, T), yte_sh)
    used = [int((D[:, k] != (1 if k == 0 else 0)).sum()) for k in range(3)]
    print(f"K={Kmax} 試行{seed}：学習用 {a_tr}/{NTR}（段ごと {hist}）  試験 {a_te}/{NTE}  ｜ 対照 {ctl}/{NTE}  ｜ 動いた桁の数（段0,1,2） {used}", flush=True)
