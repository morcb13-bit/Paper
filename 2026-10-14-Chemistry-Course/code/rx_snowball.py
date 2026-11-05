# 余った電子が雪だるま式に増える ── 形の勘定（整数だけ）
# rx_wedge.py：O2 を引き離す楔は電子のギャスケット 3 個（O1 の席 A・O2 の席 A・O1 の席 B）。
# 各 O の手の席は A と B（電子 1 個ずつ、計 4 個）。楔が入った席は電子 2 個（楔＋手）になる。
# H が 4 個、席 A・A・B・B に入ると、席にいた電子を押し出す：2・2・2・1 で 7 個。
# 1 組（2H2 + O2）で 楔 3 個を使い、電子 7 個が出る。出た電子は次の O2 の楔になる。
def seats():
    hand = {('O1','A'):1, ('O2','A'):1, ('O1','B'):1, ('O2','B'):1}
    wedge = [('O1','A'), ('O2','A'), ('O1','B')]
    for s in wedge: hand[s] += 1
    H = [('O1','A'), ('O2','A'), ('O1','B'), ('O2','B')]
    return len(wedge), [hand[s] for s in H]
def gens(need, emit, start, G=10):
    pool, out = start, []
    for _ in range(G):
        n, r = 0, pool
        while r >= need: r -= need; n += 1         # 割り算を使わず、引けるだけ引く
        out.append(n)
        if n == 0: break
        pool = r
        for _ in range(n): pool += emit
    return out
o = []; P = lambda *a: o.append(' '.join(map(str, a)))
need, pushed = seats(); emit = sum(pushed)
P('楔（要る電子）', need, '  H が押し出す電子', pushed, '計', emit)
P('S1 出る電子が要る電子より多い', emit, '>', need, 'OK' if emit > need else 'NG')
g = gens(need, emit, need); P('S2 火花 3 から、世代ごとの O2 の数', g, 'OK' if all(g[i] < g[i+1] for i in range(len(g)-1)) else 'NG')
c = gens(need, need, need); P('負の対照 H が楔だけを押し出す（出る 3）:', c[:10], '増えない' if max(c) == 1 else '対照として失格')
open('rx_snowball_out.txt', 'w').write('\n'.join(o) + '\n'); print('\n'.join(o))
