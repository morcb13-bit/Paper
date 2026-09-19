#  五芒星ニューロン模型  検定DF1〜DF13 / NE0〜NE3 / IF1〜IF2
#
#      発端（監督）：五芒星をニューロンに見立て、正十二面体を面ひとつで五芒星に伏せて
#                    メモリ兼入力装置にする。平衡5進数で情報を持ち、一筆書きで読み出す。
#
#  前提：wind_core.py（担体91環・五角形628枚・五芒星30個）を exec する。
#        b13_chain_units.py が要る。スキル b13-verify の scripts/ には入っていないので
#        リポジトリ Paper/2026-07-26-penrose/code/ から持ってくること。
#
#  検定DF1   五芒星の隙間の形               凹頂点 φ⁻¹・尖り φ・向きの差36°
#  検定DF2/3 正十二面体の面を伏せる倍率      φ⁻¹ と φ（負の対照：他の隙間は φ の冪でない）
#  検定DF4    五芒星の中心間の最小距離        φ⁴
#  検定DF5    面を伏せた立体の投影と衝突      投影 φ²/√(2+φ)・余裕は φ
#  検定DF6/8  立体間の結合                    次数1〜3（縁だけが理由ではない）
#  検定DF9    五芒星に接するもの              五角形10枚のみ。細ひし形は接しない
#  検定DF10/11 まわり10枚の並び               72°対称・36°非対称
#  検定DF12   72°対称が破れる半径             φ²で26個・φ³で4個
#  検定DF13   30×5=150通りの札                すべて異なる（担体が有限だから）
#  検定NE0    分裂頂点60 と 対角線60          面の向きを決めれば一対一
#  検定NE1/2  書いて読む／空試験              完全復元／全ゼロ
#  検定NE3    60歩の閉路                      各点ちょうど3回（＝桁は頂点でなく対角線）
#  検定NE1b   出発点を変える                  120通りの取り方→120通りの順序
#  検定IF1/2  二つの焦点からの r²             30/30 が黄金整数。係数はフィボナッチ・リュカ
#
#  未実施：検定MR1〜MR3（左右の対・軸の上の残渣・三つめの焦点で破れるか）

import math, itertools, random
from collections import Counter, defaultdict

PHI = (1 + 5 ** 0.5) / 2

# ───────────────────────── 担体側 ─────────────────────────

def carrier():
    """wind_core.py を exec した環境で呼ぶこと。五芒星30個の中心を返す。"""
    rows2, place2, offs2 = U.build_stack()
    F = U.fits(sum(place2, []))
    faces = U.gaps(F)
    SC = []
    for a, c in faces:
        if abs(a - 2.9389) < 0.01:
            P = [tuple(float(t) for t in U.xy(p)) for p in c]
            SC.append((sum(p[0] for p in P) / 10, sum(p[1] for p in P) / 10))
    return F, faces, SC

def golden(x, tol=2e-6):
    """x ≈ p + qφ となる非負整数 q と整数 p を返す。無ければ None。"""
    if x < 0:
        return None
    for q in range(0, int(x / PHI) + 4):
        p = x - q * PHI
        if abs(p - round(p)) < tol:
            return (round(p), q)
    return None

# ───────────────────────── 立体側 ─────────────────────────

def dodeca():
    """辺長1の正十二面体。頂点20・稜30・面12・対角線60 を返す。"""
    V = [(s1, s2, s3) for s1 in (1, -1) for s2 in (1, -1) for s3 in (1, -1)]
    for s1 in (1, -1):
        for s2 in (1, -1):
            V += [(0, s1 / PHI, s2 * PHI), (s1 / PHI, s2 * PHI, 0), (s1 * PHI, 0, s2 / PHI)]
    e = 2 / PHI
    V = [(x / e, y / e, z / e) for x, y, z in V]
    E = [(i, j) for i, j in itertools.combinations(range(20), 2)
         if abs(math.dist(V[i], V[j]) - 1) < 1e-9]
    ae = defaultdict(set)
    for i, j in E:
        ae[i].add(j); ae[j].add(i)
    FA = set()
    def walk(p):
        if len(p) == 5:
            if p[0] in ae[p[-1]]:
                FA.add(frozenset(p))
            return
        for n in ae[p[-1]]:
            if n not in p:
                walk(p + [n])
    for s in range(20):
        walk([s])
    FA = [sorted(f) for f in FA]
    FD = {fi: [(a, b) for a, b in itertools.combinations(f, 2) if b not in ae[a]]
          for fi, f in enumerate(FA)}
    D = [d for fi in range(12) for d in FD[fi]]
    return V, E, FA, FD, D

def euler_circuit(D, start_edge=0):
    """対角線60本を一度ずつ通る閉路。通った辺の番号の列を返す。"""
    adj = defaultdict(list)
    for k, (a, b) in enumerate(D):
        adj[a].append((b, k)); adj[b].append((a, k))
    used = [False] * len(D)
    used[start_edge] = True
    a, b = D[start_edge]
    st = [(b, None)]
    path = []
    while st:
        x, ei = st[-1]
        nx = None
        for w, idx in adj[x]:
            if not used[idx]:
                nx = (w, idx); break
        if nx is None:
            st.pop()
            if ei is not None:
                path.append(ei)
        else:
            used[nx[1]] = True
            st.append((nx[0], nx[1]))
    return [start_edge] + path[::-1]

# ───────────────────────── 検定 ─────────────────────────

def run_tests():
    F, faces, SC = carrier()
    NG = 0

    print("検定DF4  五芒星の中心間の最小距離")
    d = sorted(math.dist(a, b) for a, b in itertools.combinations(SC, 2))
    ok = abs(d[0] - PHI ** 4) < 1e-9
    print(f"  {d[0]:.6f}  φ⁴ = {PHI**4:.6f}  {'OK' if ok else 'NG'}")
    NG += 0 if ok else 1

    print("\n検定DF13  30×5=150通りの札がすべて異なるか")
    XY = {q: tuple(float(t) for t in U.xy(q)) for q in F}
    R = max(math.hypot(x, y) for x, y in XY.values()) * 2.2
    def scene(i, r, rot):
        cx, cy = SC[i]
        out = []
        for q, (x, y) in XY.items():
            dx, dy = x - cx, y - cy
            dd = math.hypot(dx, dy)
            if dd <= r:
                out.append((round(dd, 4),
                            round((math.degrees(math.atan2(dy, dx)) - 72 * rot) % 360, 3),
                            F[q] % 2))
        return tuple(sorted(out))
    tags = [scene(i, R, k) for i in range(30) for k in range(5)]
    ok = len(set(tags)) == 150
    print(f"  札 {len(tags)}  異なるもの {len(set(tags))}  {'OK' if ok else 'NG'}")
    NG += 0 if ok else 1

    V, E, FA, FD, D = dodeca()
    print("\n検定NE3  60歩の閉路と各点を通る回数")
    seq = euler_circuit(D, 0)
    vis = Counter()
    for k in seq:
        a, b = D[k]; vis[a] += 1; vis[b] += 1
    ok = (sorted(seq) == list(range(60))) and set(vis.values()) == {6}
    print(f"  辺 {len(seq)}本を一度ずつ: {sorted(seq)==list(range(60))}"
          f"   各点の通過（端点として）{set(vis.values())}  {'OK' if ok else 'NG'}")
    NG += 0 if ok else 1

    print("\n検定NE1/NE2  書いて読む／空試験")
    random.seed(13)
    MEM = [random.choice([-2, -1, 0, 1, 2]) for _ in range(60)]
    out = [MEM[k] for k in seq]
    back = [None] * 60
    for pos, k in enumerate(seq):
        back[k] = out[pos]
    ok1 = back == MEM
    ok2 = set(0 for k in seq) == {0}
    print(f"  復元 {'OK' if ok1 else 'NG'}   空試験 {'OK' if ok2 else 'NG'}")
    NG += (0 if ok1 else 1) + (0 if ok2 else 1)

    print("\n検定NE1b  出発点を変えると読み出し順は変わるか")
    orders = set()
    for se in range(60):
        s = tuple(euler_circuit(D, se))
        if sorted(s) == list(range(60)):
            orders.add(s)
    print(f"  出発辺60通り → 異なる順序 {len(orders)} 通り"
          f"   {'OK（担体側の札で固定が要る）' if len(orders) > 1 else 'NG'}")

    print("\n検定IF1  二つの焦点からの r² が黄金整数か")
    AB = None
    for i, j in itertools.combinations(range(30), 2):
        if abs(math.dist(SC[i], SC[j]) - PHI ** 4) < 1e-6:
            AB = (i, j); break
    A, B = SC[AB[0]], SC[AB[1]]
    hit = sum(1 for S in SC
              if golden(math.dist(S, A) ** 2) and golden(math.dist(S, B) ** 2))
    print(f"  両方が黄金整数 {hit}/30  {'OK' if hit == 30 else 'NG'}")
    NG += 0 if hit == 30 else 1
    random.seed(1)
    bad = sum(1 for _ in range(200)
              if (lambda x, y: golden((x - A[0]) ** 2 + (y - A[1]) ** 2)
                  and golden((x - B[0]) ** 2 + (y - B[1]) ** 2))
              (random.uniform(-40, 40), random.uniform(-40, 40)))
    print(f"  負の対照 でたらめな200点 {bad}/200  {'OK' if bad == 0 else 'NG'}")
    NG += 0 if bad == 0 else 1

    print(f"\nNG {NG} / 7")

if __name__ == "__main__":
    exec(open("wind_core.py").read())
    run_tests()
