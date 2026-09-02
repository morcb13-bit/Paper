#  検定WK2  1a/1b を「折れ角」で選ぶ
#
#      WK1（角の数で右左を選ぶ）は 20本とも 2〜6歩で4またに当たって止まった。
#      次数は 2:260 / 3:620 / 4:200 で、4または18.5%。例外ではない。
#
#      そこで選び方を、角の数ではなく**折れ角**にする。確定事項から一つだけ借りる：
#          奇数の鎖は連続接続だけ・折れ角は全箇所108°（V字）。V字＝1周期。
#      これは b13-verify §鎖の平衡3進 にある確定で、新しいつまみではない。
#      スキルには「108°で曲がり続ける経路だけを選ぶ規則で覆えるか」が未測定として載っている。
#
#      1a = 108° を右へ折る／1b = 108° を左へ折る。歩くのは円環中心のグラフ。
#
#      測る量：閉じるか・歩数・通った接続の種類（連続／2つ飛ばし）の列
#      OK なら：二本が止まらずに進み、閉じるか規則的に分かれる
#      NG なら：108° が選べない場所で止まる（＝面には持ち越せない）
#      必ず落ちる設定：右左を入れ替えると二本の役が入れ替わるだけであること
#      負の対照：108° ではなく 144° で同じことをすると別の結果になること

import json, math
from collections import Counter, defaultdict
import b13_chain_units as U

d = json.load(open("rings_integer.json"))
R = [tuple(r) for r in d["rings"]]
XY = {r: U.xy(r) for r in R}

N_CONT = U.norm2(U.CONT[0]); N_SKIP = U.norm2(U.SKIP[0])
adj = defaultdict(dict)
for i, a in enumerate(R):
    for b in R[i + 1:]:
        n = U.norm2(U.zsub(a, b))
        if n == N_CONT:   adj[a][b] = "連"; adj[b][a] = "連"
        elif n == N_SKIP: adj[a][b] = "飛"; adj[b][a] = "飛"
print("円環 %d 個／次数 %s" % (len(R), dict(sorted(Counter(len(v) for v in adj.values()).items()))))
print("接続 連 %d 本・飛 %d 本\n"
      % (sum(1 for a in adj for b in adj[a] if adj[a][b] == "連") // 2,
         sum(1 for a in adj for b in adj[a] if adj[a][b] == "飛") // 2))

def turn(prev, cur, nxt):
    """prev→cur→nxt の折れ角（度）と向き（+が反時計）"""
    ax, ay = XY[cur][0] - XY[prev][0], XY[cur][1] - XY[prev][1]
    bx, by = XY[nxt][0] - XY[cur][0], XY[nxt][1] - XY[cur][1]
    cr = ax * by - ay * bx; dt = ax * bx + ay * by
    a = math.degrees(math.atan2(cr, dt))
    return 180.0 - abs(a), (1 if cr > 0 else -1)

def walk(prev, cur, side, want=108.0, limit=200):
    path = [prev, cur]; kinds = [adj[prev][cur]]
    for _ in range(limit):
        cand = []
        for w in adj[cur]:
            if w == prev: continue
            t, s = turn(prev, cur, w)
            if abs(t - want) < 1e-6 and s == side:
                cand.append(w)
        if not cand:
            return path, kinds, "選べない"
        nxt = cand[0]
        kinds.append(adj[cur][nxt])
        prev, cur = cur, nxt
        path.append(cur)
        if len(path) > 3 and cur == path[1] and prev == path[0]:
            return path, kinds, "閉じた"
    return path, kinds, "打ち切り"

def report(want):
    print("■ 折れ角 %.0f° で選ぶ" % want)
    out = []
    starts = sorted(R, key=lambda r: math.hypot(*XY[r]))[:5]   # 内側の5環から出す
    for a in starts:
        for b in adj[a]:
            for side, nm in ((+1, "1a"), (-1, "1b")):
                p, k, e = walk(a, b, side, want)
                out.append((nm, len(p), "".join(k), e))
    ends = Counter(e for _, _, _, e in out)
    print("  出発 %d 通り／終わり方 %s" % (len(out), dict(ends)))
    lens = [n for _, n, _, e in out]
    print("  歩数 最小%d 最大%d 中央%d" % (min(lens), max(lens), sorted(lens)[len(lens) // 2]))
    for nm, n, k, e in out[:6]:
        print("    %s 歩数%3d  %s  %s" % (nm, n, k[:28], e))
    same = sum(1 for i in range(0, len(out), 2) if out[i][2] == out[i + 1][2])
    print("  1a と 1b が同じ列になった出発 %d / %d\n" % (same, len(out) // 2))

report(108.0)
report(144.0)      # 負の対照
