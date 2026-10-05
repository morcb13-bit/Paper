# 検定D5：荷物に貼られた札（QR の代わり）の場所を探し当てて読む。学習した脳は使わない。判定は番地の白黒の一致と数え上げだけ。
# 荷物の面：扇10枚の担体のうち、中心から半径 40 以内（五芒星 36 個）。どの番地も白か黒。汚れとして 3 割の番地をでたらめに黒くする。
# 札：面の五芒星の一つ g に貼る（向きは 72° 刻みでどれでも）。番地は g からの相対位置（Z[ζ] の整数）で決まる。
#   目印    距離 φ の 5 枚 ＝ すべて黒
#   静かな帯 距離 φ² の 5 枚 ＝ すべて白
#   符号    距離 φ³ の 5 枚 ＝ 白黒の巡回語（5 ビット。72° 回しても裏返しても同じ番号になる巡回の型 8 通り）
#   写し    距離 5.626 の 10 枚のうち片方の 5 枚 ＝ 符号と同じ並び（読みの検査）
# 探し方：視線（網膜）を中心の五芒星から近い順に五芒星へ跳ばし、そこに目印・静かな帯・写しの一致がそろったら符号を読む。
#         面の五芒星を全部見ても見つからなければ「見つからない」として巣へ戻る（迷ったら戻る）。
# 基準（走らせる前に決めたもの）
#   合格  4000 回のうち、違う番号を読んだ回数 0
#   あわせて測る：見つけるまでの視線の跳躍の回数、見つからなかった回数
#   負の対照  静かな帯と写しの検査を外す（目印だけで読む）と、違う番号を読む回数が 0 でなくなる
#   NG  違う番号を 1 回でも読む
import sys, pickle, math, random
from collections import Counter
sys.path.insert(0, '/home/claude/Paper/2026-09-21-plant-sun/code')
import b13_chain_units as U
cells, z0 = pickle.load(open('/home/claude/icolearn/cells10.pkl', 'rb'))
S = set(cells); cx, cy = U.xy(z0)
rad = lambda z: math.hypot(U.xy(z)[0] - cx, U.xy(z)[1] - cy)
stars = {}
for q in S:
    for m in range(10):
        g = U.zsub(q, U.zmul(U.PHI, U.zt(m)))
        if g in stars or g in S: continue
        for k in (0, 1):
            if all(U.zadd(g, U.zmul(U.PHI, U.zt(k + 2 * i))) in S for i in range(5)): stars[g] = k
FACE = [q for q in S if rad(q) < 40]
G = sorted([g for g in stars if rad(g) < 40 - 7], key=lambda g: (round(rad(g), 6), U.xy(g)))
k0 = stars[z0]
V = {'目印': (1, 0, 1, 0), '帯': (-2, 1, -1, 2), '符号': (-2, -1, -2, 0), '写しA': (-4, 1, -2, 4)}
def ring(g, name, turn):
    """g のまわりの name の 5 枚。turn は 72° 刻みの貼る向き（0〜4）"""
    k = stars[g] - k0; v = V[name]
    return [U.zadd(g, U.zmul(v, U.zt(k + 2 * ((i + turn) % 5)))) for i in range(5)]
def necklace(bits):
    rots = [tuple(bits[i:] + bits[:i]) for i in range(5)]
    rots += [tuple(reversed(r)) for r in rots]
    return min(rots)
CODES = sorted({necklace([(n >> i) & 1 for i in range(5)]) for n in range(32)})

def make_face(rng, noise=3):
    black = {q for q in FACE if rng.randrange(10) < noise}
    g = rng.choice(G); turn = rng.randrange(5); code = list(rng.choice(CODES))
    s = rng.randrange(5); code = code[s:] + code[:s]                     # 巡回語をどこから書くかもでたらめ
    for q in ring(g, '目印', turn): black.add(q)
    for q in ring(g, '帯', turn): black.discard(q)
    for q, b in zip(ring(g, '符号', turn), code): (black.add if b else black.discard)(q)
    for q, b in zip(ring(g, '写しA', turn), code): (black.add if b else black.discard)(q)
    return black, g, necklace(code)

def read_at(black, g, strict=True):
    if not all(q in black for q in ring(g, '目印', 0)): return None
    if strict and any(q in black for q in ring(g, '帯', 0)): return None
    code = [int(q in black) for q in ring(g, '符号', 0)]
    if strict:
        copy = [int(q in black) for q in ring(g, '写しA', 0)]
        if copy != code: return None
    return necklace(code)

def search(black, strict=True):
    for jumps, g in enumerate(G, 1):
        r = read_at(black, g, strict)
        if r is not None: return r, g, jumps
    return None, None, len(G)

if __name__ == "__main__":
    print(f"面の五角形 {len(FACE)}、札を貼れる五芒星 {len(G)}、番号（巡回の型）{len(CODES)} 通り")
    for name, strict in (('目印＋静かな帯＋写し', True), ('目印だけ（対照）', False)):
        rng = random.Random(4242); C = Counter()
        for _ in range(4000):
            black, g, code = make_face(rng)
            r, gg, jumps = search(black, strict)
            C['跳躍'] += jumps
            if r is None: C['見つからない'] += 1
            elif r == code and gg == g: C['正しい'] += 1
            else: C['違う番号'] += int(r != code); C['違う場所'] += int(gg != g)
        print(f"{name}：正しい {C['正しい']}  違う番号 {C['違う番号']}  違う場所で読んだ {C['違う場所']}  見つからない {C['見つからない']}  視線の跳躍の合計 {C['跳躍']}")
