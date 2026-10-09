# 対象：60環担体の五角形420枚を番地とし、辺を共有する五角形（中心間の二乗長 = N_CELL）を隣とする。
#       その上で、＋−の偏り1個が並べ替えなしに全体を片側で満たすか。
#
#  事前登録（整数と比較だけで判定）
#  検定P0 対称初期（＋をc、−を鏡映conj(c)に置く対を3組）で D=0 が最後まで保たれる
#      OK なら：規則と担体が鏡映について偏りを持たない（実装確認）
#      NG なら：担体か規則に向きの偏りがある。以後の結果は読まない
#  検定P1 鏡像：＋を1個足した場合と−を同じ番地に1個足した場合で、最終の (N+,N-) が入れ替わる
#      OK なら：結果は型の名前に依らない
#  検定P2 負の対照（写さず、組だけ除く）で F1 が NG
#      OK(=NGが出る) なら：増幅は写す段が担っている
#  検定P3 本題：写す＋組を除く、並べ替えなし。1個を足す番地を空いた全番地で掃く
#      F2強（20*同側 >= 19*占有、D!=0、占有 >= 初期占有）が通った番地の数を数える
#      全番地で通る  → 担体の隣接だけで「混ざり」が供給される
#      全番地で落ちる → 担体の隣接は混ざりを供給しない（輪と同じ）
#      一部だけ通る  → 足す番地（着地）が結果を決める
#  検定P4 正の対照：P3 に型を見ない並べ替え（番地をk飛び）を加えると通る
#      NG なら：P3 の NG は規則か担体の大きさの問題で、混ざりの問題ではない
import json, sys, io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    import b13_two_tilings as T

d = json.load(open("rings_integer.json"))
cells = [tuple(int(t) for t in k.split(",")) for k in d["cells"].keys()]
cells.sort()
idx = {c: i for i, c in enumerate(cells)}
N = len(cells)

nb = [[] for _ in range(N)]
for i in range(N):
    for j in range(i + 1, N):
        if T.norm2(T.zsub(cells[i], cells[j])) == T.N_CELL:
            nb[i].append(j); nb[j].append(i)

# 鏡映（共役）で担体が閉じているか
conj = [idx.get(tuple(T.zconj(c))) for c in cells]
closed = all(c is not None for c in conj)

def step(x, copy=True, pair=True, k=0):
    y = list(x)
    if copy:
        for i in range(N):
            if x[i] == 0:
                s = 0
                for j in nb[i]: s += x[j]
                if s > 0: y[i] = 1
                elif s < 0: y[i] = -1
    if pair:
        z = list(y)
        for i in range(N):
            if y[i] != 0:
                for j in nb[i]:
                    if y[j] == -y[i]: z[i] = 0; break
        y = z
    if k:
        w = [0] * N; j = 0
        for i in range(N):
            w[i] = y[j]; j += k
            while j >= N: j += -N
        y = w
    return y

def cnt(x):
    p = sum(1 for v in x if v > 0); m = sum(1 for v in x if v < 0)
    return p, m

def run(x, steps=400, **kw):
    p, m = cnt(x); D0, occ0 = p + -m, p + m
    h = []
    for _ in range(steps):
        x = step(x, **kw); p, m = cnt(x); h.append((p, m, p + -m))
    return D0, occ0, h

def judge(D0, occ0, h):
    a0 = D0 if D0 >= 0 else -D0
    Ds = [t[2] for t in h]; f1 = False
    if a0 > 0:
        for t in range(len(Ds)):
            if all((d if d >= 0 else -d) >= a0 + a0 for d in Ds[t:]): f1 = True; break
    R, L, D = h[-1]; occ = R + L; same = R if D >= 0 else L
    f2 = occ > 0 and D != 0 and 20 * same >= 19 * occ and occ >= occ0
    return f1, f2

# 種：鏡映で写り合う対を3組（自分自身に写る番地は避ける）
cand = [i for i in range(N) if conj[i] is not None and conj[i] != i]
pairs = []; used = set()
for i in cand[::37]:
    if i in used or conj[i] in used: continue
    pairs.append(i); used |= {i, conj[i]}
    if len(pairs) == 3: break

def seed(extra=None, sign=1):
    x = [0] * N
    for i in pairs: x[i] = 1; x[conj[i]] = -1
    if extra is not None: x[extra] = sign
    return x

deg = {}
for i in range(N): deg[len(nb[i])] = deg.get(len(nb[i]), 0) + 1
print(f"番地 {N}  隣の数の分布 {dict(sorted(deg.items()))}  鏡映で閉じる {closed}")
print(f"種の対 {[(i, conj[i]) for i in pairs]}")

D0, o0, h = run(seed())
print(f"P0 対称初期: 最終 {h[-1]}  D=0 が全歩で保たれる {'OK' if all(t[2]==0 for t in h) else 'NG'}")

free = [i for i in range(N) if seed()[i] == 0]
ok3 = mir = f1n = 0; ok2 = 0; ok4 = 0; finals = {}
for e in free:
    Dp, op, hp = run(seed(e, 1))
    Dm, om, hm = run(seed(e, -1))
    f1, f2 = judge(Dp, op, hp)
    if f2: ok3 += 1
    if f1: f1n += 1
    if hp[-1][0] == hm[-1][1] and hp[-1][1] == hm[-1][0]: mir += 1
    finals[hp[-1]] = finals.get(hp[-1], 0) + 1
    Dn, on, hn = run(seed(e, 1), copy=False)
    if judge(Dn, on, hn)[0]: ok2 += 1
nf = len(free)
print(f"P1 鏡像の入れ替わり: {mir}/{nf}")
print(f"P2 負の対照（写さない）で F1 が通った番地: {ok2}/{nf}  （0 なら OK）")
print(f"P3 並べ替えなし: F1 通過 {f1n}/{nf}  F2強 通過 {ok3}/{nf}")
top = sorted(finals.items(), key=lambda t: -t[1])[:6]
print(f"   最終 (N+,N-,D) の多いもの {top}")
for k in (7, 11, 13, 17):
    if any(N % q == 0 and k % q == 0 for q in (2, 3, 5, 7)): continue
    ok4 = 0
    for e in free:
        Dp, op, hp = run(seed(e, 1), k=k)
        if judge(Dp, op, hp)[1]: ok4 += 1
    print(f"P4 並べ替え k={k}: F2強 通過 {ok4}/{nf}")
