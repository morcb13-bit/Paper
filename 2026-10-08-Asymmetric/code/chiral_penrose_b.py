# 追加：P1を正しい鏡像（＋をe ↔ −をconj(e)、符号反転と鏡映の組）で取り直す。
#       P3の最終の片寄り（同側の割合）を整数の段で数え、担体でない隣接（輪・方眼）と比べる。
#  事前登録
#  検定P1' ＋をeに足した結果と、−をconj(e)に足した結果で (N+,N-) が入れ替わる → 全番地でOKなら鏡像再現
#  検定P5 同側の割合の段（20*同側 >= q*占有 を満たす最大の q）の分布を、
#         担体420・輪420・方眼20x21 で比べる。種と足す番地の扱いは揃える。
#      担体だけが高い段に寄る → 担体の隣接が混ざりの一部を供給している
#      輪や方眼と変わらない   → 担体固有の効果は無い
import io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    exec(open("chiral_penrose.py").read().split('deg = {}')[0])

def level(h):
    R, L, D = h[-1]; occ = R + L; same = R if D >= 0 else L
    q = 0
    while q < 20 and 20 * same >= (q + 1) * occ: q += 1
    return q

def sweep(nbs, mirror, seeds_pairs, label):
    global nb, N
    nb = nbs; N = len(nbs)
    def sd(e=None, s=1):
        x = [0] * N
        for i in seeds_pairs: x[i] = 1; x[mirror[i]] = -1
        if e is not None: x[e] = s
        return x
    D0, o0, h0 = run(sd())
    fr = [i for i in range(N) if sd()[i] == 0]
    dist = {}; mir = 0
    for e in fr:
        _, _, hp = run(sd(e, 1)); _, _, hm = run(sd(mirror[e], -1))
        if hp[-1][0] == hm[-1][1] and hp[-1][1] == hm[-1][0]: mir += 1
        q = level(hp); dist[q] = dist.get(q, 0) + 1
    print(f"{label}: 対称初期 D=0 {'OK' if all(t[2]==0 for t in h0) else 'NG'}  鏡像 {mir}/{len(fr)}  "
          f"段の分布(20*同側>=q*占有 の最大q) {dict(sorted(dist.items()))}")

pen_nb = nb
sweep(pen_nb, conj, pairs, "担体420")
# 輪420：鏡映 i -> 419-i、種は同じ添字
M = 420
ring_nb = [[(i - 1) % M, (i + 1) % M] for i in range(M)]
sweep(ring_nb, [M - 1 - i for i in range(M)], pairs, "輪420  ")
# 方眼 20x21：鏡映は左右反転
W, H = 20, 21
g_nb = []
for i in range(W * H):
    r, c = divmod(i, W); l = []
    if c > 0: l.append(i - 1)
    if c < W - 1: l.append(i + 1)
    if r > 0: l.append(i - W)
    if r < H - 1: l.append(i + W)
    g_nb.append(l)
gm = [(i // W) * W + (W - 1 - i % W) for i in range(W * H)]
gp = [p for p in pairs if gm[p] != p]
sweep(g_nb, gm, gp, "方眼420")
