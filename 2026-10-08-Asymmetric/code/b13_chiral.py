"""第20章　右と左の一個差：不斉の符号と、偏りが広がる条件。整数と足し算と比較だけ。"""
from itertools import permutations, product

Vec = tuple[int, int, int]


# ── 右か左か：四面体の向きの符号 ──────────────────────────
def sub(u: Vec, v: Vec) -> Vec:
    return (u[0] + -v[0], u[1] + -v[1], u[2] + -v[2])


def det3(u: Vec, v: Vec, w: Vec) -> int:
    return (u[0] * (v[1] * w[2] + -(v[2] * w[1]))
            + -(u[1] * (v[0] * w[2] + -(v[2] * w[0])))
            + u[2] * (v[0] * w[1] + -(v[1] * w[0])))


# 硤合反応の生成物の中心の炭素から見た4本の腕の先（0.01Å を 1 とする整数）
ARMS = {"H": (33, -43, -95), "OH": (27, -99, 100), "iPr": (83, 127, 30), "環": (-148, 23, -4)}


def hand(m: dict) -> int:
    """四面体の向き det[OH−H, iPr−H, 環−H]。符号が右か左かを決める。"""
    return det3(sub(m["OH"], m["H"]), sub(m["iPr"], m["H"]), sub(m["環"], m["H"]))


def mirror(m: dict) -> dict:
    return {k: (-v[0], v[1], v[2]) for k, v in m.items()}


def turns() -> tuple[list, list]:
    """軸の入れ替えと向きの反転からなる整数の行列48通りを、回転24と鏡映を含む24に分ける。"""
    rot, ref = [], []
    for pm in permutations(range(3)):
        for sg in product((1, -1), repeat=3):
            M = [[0, 0, 0] for _ in range(3)]
            for i in range(3):
                M[i][pm[i]] = sg[i]
            d = det3(*[(M[0][j], M[1][j], M[2][j]) for j in range(3)])
            (rot if d > 0 else ref).append(M)
    return rot, ref


def apply(M, v: Vec) -> Vec:
    return tuple(M[i][0] * v[0] + M[i][1] * v[1] + M[i][2] * v[2] for i in range(3))


def side(m: dict) -> int:
    """受け口 OH・iPr・環 の面に対して、H がどちらの側にあるか（符号）。"""
    return det3(sub(m["iPr"], m["OH"]), sub(m["環"], m["OH"]), sub(m["H"], m["OH"]))


def triangle(m: dict) -> tuple:
    """受け口の三角形の辺の二乗。"""
    sq = lambda u: u[0] * u[0] + u[1] * u[1] + u[2] * u[2]
    return (sq(sub(m["OH"], m["iPr"])), sq(sub(m["iPr"], m["環"])), sq(sub(m["環"], m["OH"])))


# ── 偏りが広がる条件：器 ─────────────────────────────────
def vessel(R: int, L: int, food: int, grow: bool, pair: bool, steps: int = 60) -> tuple:
    """どの分子も互いに出会える器。grow：自分と同じ型を写す。pair：右と左が組になって抜ける。"""
    for _ in range(steps):
        if grow:
            gr = gl = 0
            first_R = R >= L                       # 多い側から1個ずつ交互に配る（割り算を使わない）
            while food > 0 and (gr < R or gl < L):
                for side_R in ((True, False) if first_R else (False, True)):
                    if food > 0 and side_R and gr < R:
                        gr += 1; food += -1
                    if food > 0 and not side_R and gl < L:
                        gl += 1; food += -1
            R, L = R + gr, L + gl
        if pair:
            m = R if R < L else L
            R, L = R + -m, L + -m
    return R, L


# ── 偏りが広がる条件：輪 ─────────────────────────────────
N = 60


def ring_seed(extra: int) -> list:
    x = [0] * N
    for k in range(6):
        x[10 * k] = 1 if k % 2 == 0 else -1
    if extra:
        x[5] = extra
    return x


def ring_step(x: list, copy: bool, pair: bool, k: int) -> list:
    y = list(x)
    if copy:                                       # 写す：空いた席を両隣の和の符号で埋める
        for i in range(N):
            if x[i] == 0:
                s = x[i - 1] + x[(i + 1) % N]
                if s > 0: y[i] = 1
                elif s < 0: y[i] = -1
    if pair:                                       # 組で消える：隣り合った右と左を空席に戻す
        z = list(y)
        for i in range(N):
            if y[i] != 0 and (y[i - 1] == -y[i] or y[(i + 1) % N] == -y[i]):
                z[i] = 0
        y = z
    if k:                                          # 混ぜる：席を k 飛びに読み直す（型は見ない）
        w, j = [0] * N, 0
        for i in range(N):
            w[i] = y[j]; j += k
            if j >= N: j += -N
        y = w
    return y


def ring(extra: int, copy: bool, pair: bool, k: int, steps: int = 200) -> tuple:
    x = ring_seed(extra)
    for _ in range(steps):
        x = ring_step(x, copy, pair, k)
    return sum(1 for v in x if v > 0), sum(1 for v in x if v < 0)


if __name__ == "__main__":
    print("== 右か左か：四面体の向きの符号")
    R, S = ARMS, mirror(ARMS)
    print("腕の長さの二乗", {k: v[0] * v[0] + v[1] * v[1] + v[2] * v[2] for k, v in R.items()})
    print("R体", hand(R), "  S体（鏡像）", hand(S))
    rot, ref = turns()
    keep = sum(1 for M in rot if (hand({k: apply(M, v) for k, v in R.items()}) > 0) == (hand(R) > 0))
    flip = sum(1 for M in ref if (hand({k: apply(M, v) for k, v in R.items()}) > 0) != (hand(R) > 0))
    print(f"回転 {len(rot)} 通りで符号が保たれる: {keep}   鏡映を含む {len(ref)} 通りで符号が反転: {flip}")

    print("\n== 型紙：受け口3つと板の表裏")
    print("受け口の三角形（辺の二乗）R体", triangle(R), " S体", triangle(S))
    print("H の側  R体", side(R), "  S体", side(S))

    print("\n== 器（初め 右51・左50、材料 100000）")
    for name, g, p in (("写さない", False, False), ("写すだけ", True, False),
                       ("組で抜けるだけ", False, True), ("写す＋組で抜ける", True, True)):
        print(f"{name:10s} 右51・左50 → {vessel(51, 50, 100000, g, p)}   右50・左51 → {vessel(50, 51, 100000, g, p)}"
              f"   右50・左50 → {vessel(50, 50, 100000, g, p)}")

    print("\n== 輪（60席、右3・左3 に右を1個足す）")
    for name, c, p, k in (("混ぜるだけ", False, False, 7), ("写すだけ", True, False, 0),
                          ("写す＋組で消える", True, True, 0), ("写す＋混ぜる", True, False, 7),
                          ("写す＋組で消える＋混ぜる", True, True, 7)):
        print(f"{name:14s} 右を足す → {ring(1, c, p, k)}   左を足す → {ring(-1, c, p, k)}   足さない → {ring(0, c, p, k)}")
    print("混ぜる飛び幅を変える（写す＋組で消える＋混ぜる）:",
          [(k, ring(1, True, True, k)) for k in (7, 11, 13, 17)])
