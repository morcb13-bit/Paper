# 気体の中を転がる電子 ── 初期値だけ乱数、後はオートマトン（整数と加算だけ）
#
# 格子：一辺 6 の立方体（番地 i,j,k）。転がる電子（正四面体）は一歩で面の隣の立方体へ移る。
#       向きは番地の偶奇で決まる（rx_roll_parity.py：377/377）。
# 分子：O2 を 2×2×2 の立方体の塊に置く（塊の角の番地は偶数）。
#       塊の中で奇の立方体 4 つは O の正四面体（A の向き）で、電子と同じ向きで重なるので入れない。
#       偶の立方体 4 つは席（B の向き）で、電子がぴったり重なる。席は通った経路の分を溜め、電子 1 個分に届くと楔が入って埋まる（以後は入れない）。
# 反応：楔 3 個がそろった O2 は反応し、楔の入った 3 席から 2 個ずつ（楔＋手）、残る 1 席から 1 個（手）、
#       計 7 個を、それぞれの席から転がり出す（rx_snowball.py の勘定）。反応した塊は H が席を埋めるので、
#       8 つの立方体すべてに入れなくなる。
# 初期値：塊を置くかどうかだけを乱数（足し算だけの遅れフィボナッチ）で決める。
# 経路：転がりは幾何が許す倒れ方をすべてたどる。一つの経路の重みは 6^{-歩}（割り算はせず、全体を 6^t 倍で持つ）。
#       入れない立方体へは倒れられず、倒れかけて元の立方体に戻る（経路は残る）。枠の外へ出た経路は失われる。
#
# 事前に決めた基準：
#  R1 火花 1 組で、置いた O2 のすべてが反応する濃さがある
#  R2 反応数が山まで増える（加速）
#  R3 余る電子の行き先を数える（席に入った・枠の外へ出た）。入った分 ＋ 外へ出た分 ＋ 転がり中 ＝ 放った分（勘定が閉じる）
#  負の対照：放つ電子を楔の 3 個だけ（手を押し出さない）にすると、どの濃さでも全体には広がらない
class Fib:
    def __init__(self, seed):
        a, b = 1, 2 + seed; self.s = []
        while b >= 6: b -= 6
        for _ in range(55):
            a, b = b, a + b                      # フィボナッチを 6 を引いて 0..5 に保つ
            if b >= 6: b -= 6
            self.s.append(b)
        self.i = 0
    def next(self):
        s, i = self.s, self.i
        v = s[(i + 31) % 55] + s[i]
        if v >= 6: v -= 6
        s[i] = v; self.i = i + 1 if i < 54 else 0
        return v
D6 = [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]
def x6(x): return x+x+x+x+x+x
def build(B, dens, seed):
    """B：一辺に並ぶ塊の数（立方体は 2B 個）。dens：6 のうちいくつで塊を置くか"""
    rng = Fib(seed); mol = {}; c = B // 2
    for bx in range(B):
        for by in range(B):
            for bz in range(B):
                if (bx, by, bz) == (c, c, c) or rng.next() < dens:
                    mol[(bx, by, bz)] = len(mol)
    return mol
def run(B, dens, seed, wedge_out=2, hand_out=1, T=150):
    N = 2 * B; mol = build(B, dens, seed); c = B // 2
    owner = {}; seat = {}; solid = set()
    for (bx, by, bz), mid in mol.items():
        for a in (0, 1):
            for b in (0, 1):
                for d in (0, 1):
                    q = (2*bx+a, 2*by+b, 2*bz+d); owner[q] = mid
                    if (a + b + d) % 2: solid.add(q)
                    else: seat[q] = 0                      # 席の楔（6^t 倍）
    seats_of = {}
    for q, mid in owner.items():
        if q in seat: seats_of.setdefault(mid, []).append(q)
    state = {mid: 0 for mid in mol.values()}              # 0 待つ, 1 反応した
    full = set()                                          # 楔が 1 個に届いた席
    scale = 1; U = {}; hist = []; emitted = 0; captured = 0; lost = 0
    def ignite(mid):
        nonlocal emitted
        state[mid] = 1
        ss = sorted(seats_of[mid])
        wq = [q for q in ss if q in full]; hq = [q for q in ss if q not in full]
        if len(wq) < 3: wq, hq = ss[:3], ss[3:]           # 火花の組：楔 3 個が入ったとみなす
        for q in ss: blocked.add(q)
        for q, k in [(q, wedge_out) for q in wq] + [(q, hand_out) for q in hq[:1]]:
            for _ in range(k): U[q] = U.get(q, 0) + scale; emitted += 1
    blocked = set(solid)
    ignite(mol[(c, c, c)]); hist.append(1)
    for t in range(T):
        U2 = {}
        for q, u in U.items():
            for d in D6:
                r = (q[0]+d[0], q[1]+d[1], q[2]+d[2])
                if not all(0 <= x < N for x in r): lost += 1 if False else 0; U2[('out',)] = U2.get(('out',), 0) + u; continue
                if r in blocked: U2[q] = U2.get(q, 0) + u; continue   # 倒れられない：倒れかけて元に戻る
                U2[r] = U2.get(r, 0) + u
        scale = x6(scale)
        out_w = U2.pop(('out',), 0)
        for q in list(seat): seat[q] = x6(seat[q])
        # 席に入る（楔）：経路はそこで止まる
        for q, u in U2.items():                 # 席を通った経路の分を溜める（経路は止まらず進む：すべての分子に届く）
            if q in seat and q not in blocked: seat[q] += u
        U = U2
        for q, v in seat.items():
            if q not in full and v >= scale: full.add(q); blocked.add(q)
        fire = [mid for mid in mol.values() if state[mid] == 0 and sum(1 for q in seats_of[mid] if q in full) >= 3]
        for mid in fire: ignite(mid)
        hist.append(len(fire))
        if not U: break
    done = sum(state.values())
    # 勘定（電子の個数、6^t で割らずに比べる）
    cap = 0
    for q, v in seat.items(): cap += v
    moving = sum(U.values())
    return done, len(mol), hist, scale, cap, moving, emitted
out = []; P = lambda *a: out.append(' '.join(map(str, a)))
B = 7
for dens in (6, 5, 4, 3):
    res = []
    for seed in range(2):
        d, tot, h, sc, cap, mv, em = run(B, dens, seed)
        while h and h[-1] == 0: h.pop()
        res.append((d, tot, max(h), len(h)))
    P(f'濃さ {dens}/6：', '  '.join(f'種{s}: {d}/{tot} 山{m} 歩{n}' for s, (d, tot, m, n) in enumerate(res)))
open('rx_roll_gas_out.txt', 'w').write('\n'.join(out) + '\n'); print('\n'.join(out))
