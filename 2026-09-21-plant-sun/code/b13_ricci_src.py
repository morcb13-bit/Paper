# b13_ricci_src.py — 伸び縮みするセルの網に「源」（いつも +1 を出す入力セル）を一つ足す
# 見えたこと（ric_probe.py）：入力1・2だけが効く網でも、関係ない入力が 0（黙り）だと信号が出力まで届かず、
#   関係ない入力が 4 つとも鳴っているときだけ「同符号」が正しく出ていた。関係ない入力が、運び手の役をしていた。
# そこで、題と無関係にいつも鳴る源を一つ置く（PVP の源 B と同じ役）。それ以外は b13_ricci.py と同じ。
# 学習の一手に「閉じたまま止める」（ω=0・位相を偶数に）を足す。外れが同じなら止める方を採る。
# 基準：R1 試験 6 割以上／対照 半分前後／R2 効く入力が 1・2 だけ（同符号）、1 だけ（一つ目が正）
#       R3 全問正解（試験 224/224）が出る試行があるか
import itertools, random, sys
import numpy as np
import b13_fannet as F, b13_spinnet as S, b13_ricci as R
R.CLOSE['on'] = True

if __name__ == "__main__":
    net = F.make_floor(m=7)
    ALL0 = np.array([s for s in itertools.product((-1, 0, 1), repeat=6) if any(s)])
    NETS = {"A 一方向": (net['fwd'], net['L'] + 1), "C 隣＋φ²": (F.directed(net['nb'] + net['phi2']), net['L'] + 5)}
    seeds = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    NTR = int(sys.argv[2]) if len(sys.argv) > 2 else 100
    task = sys.argv[3]; only_net = sys.argv[4] if len(sys.argv) > 4 else None
    f, sel = S.TASKS[task]; A = sel(ALL0); y = f(A)
    A7 = np.hstack([A, np.ones((len(A), 1), dtype=A.dtype)])      # 7 番目＝源
    print(f"題「{task}」 源つき、入力 {len(A)} 通り、学習用 {NTR}、{seeds} 回", flush=True)
    for nname, (E, T) in NETS.items():
        if only_net and not nname.startswith(only_net): continue
        for s in range(seeds):
            rng = random.Random(5000 + s); idx = list(range(len(A))); rng.shuffle(idx)
            tr, te = np.array(idx[:NTR]), np.array(idx[NTR:])
            p, w = R.learn(net, E, A7[tr], y[tr], T, rng)
            a_tr = S.score(R.run(net, E, p, w, A7[tr], T), y[tr]); a_te = S.score(R.run(net, E, p, w, A7[te], T), y[te])
            eff = R.effective(net, E, p, w, A7, T)[:6]
            # 関係ない入力を全部 0 にしたときと、そのままとで出力が変わる例の数（運び手になっていないか）
            Az = A7.copy(); Az[:, 2:6] = 0
            dep = int((R.run(net, E, p, w, Az, T) != R.run(net, E, p, w, A7, T)).sum())
            perm = list(range(len(A))); rng.shuffle(perm); ysh = y[perm]
            pc, wc = R.learn(net, E, A7[tr], ysh[tr], T, rng)
            ctl = S.score(R.run(net, E, pc, wc, A7[te], T), ysh[te])
            print(f"  {nname} {s}：学習用 {a_tr}/{NTR}  試験 {a_te}/{len(te)}  ｜ 対照 {ctl}/{len(te)}  ｜ "
                  f"効く入力 {eff}  関係ない入力を0にすると変わる {dep}/{len(A)}", flush=True)
