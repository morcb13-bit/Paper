"""
cone_fft_complete.py
cone-FFT v28 接続完成版

構成:
  1. BASE=3120 の整数 cos 射影
  2. T/U 偶奇分割
  3. Parseval 本体: 40 再帰
  4. cone 局所補正: C_CONE = 2/phi
  5. 固定小数点整数で合成出力

注意:
- Parseval 本体は厳密整数
- cone 補正 2/phi は無理数なので、最後に BASE 固定小数点で適用
- したがって「コアは整数」「補正は BASE スケール整数近似」という構成
"""

from __future__ import annotations
import math
from dataclasses import dataclass


# ─────────────────────────────────────────
# 定数
# ─────────────────────────────────────────
PHI = (1 + 5**0.5) / 2
BASE = 3120

CENTER = (7, 7)
ROTATIONS = [(3, 4), (4, -3), (-3, -4), (-4, 3)]
ORBIT = [(CENTER[0] + dx, CENTER[1] + dy) for dx, dy in ROTATIONS]

RHO_A = 5**0.5
R = (1 - RHO_A) / (1 + RHO_A)      # -phi^-2
C_CONE = 2 + 2 * R                 # 2/phi

C_CONE_Q = round(BASE * C_CONE)

COS_TABLE = [round(BASE * math.cos(2 * math.pi * k / BASE)) for k in range(BASE)]


# ─────────────────────────────────────────
# 基本演算
# ─────────────────────────────────────────
def dot(a: list[int], b: list[int]) -> int:
    return sum(x * y for x, y in zip(a, b))

def norm2(a: list[int]) -> int:
    return dot(a, a)

def split_TU(f: list[int]) -> tuple[list[int], list[int]]:
    return f[0::2], f[1::2]

def merge_TU(T: list[int], U: list[int]) -> list[int]:
    n = len(T) + len(U)
    out = [0] * n
    out[0::2] = T
    out[1::2] = U
    return out

def r_global(T: list[int], U: list[int]) -> float:
    den = norm2(T) * norm2(U)
    return dot(T, U) / math.sqrt(den) if den else 0.0


# ─────────────────────────────────────────
# 1. BASE=3120 整数 cos 射影
# ─────────────────────────────────────────
def int_cos_project_3120(f: list[int]) -> list[int]:
    """
    3120周期 cos テーブルによる整数射影。
    出力は BASE スケール付き整数。
    """
    N = len(f)
    return [
        sum(f[j] * COS_TABLE[(j * k) % BASE] for j in range(N))
        for k in range(N)
    ]


# ─────────────────────────────────────────
# 2. Parseval 本体: 40 再帰
# ─────────────────────────────────────────
def step_rule_40(T: list[int], U: list[int]) -> tuple[list[int], list[int]]:
    T_next = [6 * t + 2 * u for t, u in zip(T, U)]
    U_next = [2 * t - 6 * u for t, u in zip(T, U)]
    return T_next, U_next

def build_fractal_integer_signal(B: int, T0: list[int], U0: list[int]) -> tuple[list[int], list[int]]:
    T, U = T0[:], U0[:]
    for _ in range(B):
        T, U = step_rule_40(T, U)
    return T, U


# ─────────────────────────────────────────
# 3. cone 局所固有量
# ─────────────────────────────────────────
def sub2(a: tuple[int, int], b: tuple[int, int]) -> tuple[int, int]:
    return (a[0] - b[0], a[1] - b[1])

def dot2(a: tuple[int, int], b: tuple[int, int]) -> int:
    return a[0] * b[0] + a[1] * b[1]

def gram2(a: tuple[int, int], b: tuple[int, int]) -> tuple[tuple[int, int], tuple[int, int]]:
    aa = dot2(a, a)
    ab = dot2(a, b)
    bb = dot2(b, b)
    return ((aa, ab), (ab, bb))

def eigvals_sym2(G: tuple[tuple[int, int], tuple[int, int]]) -> tuple[float, float]:
    a, b = G[0]
    _, d = G[1]
    tr = a + d
    det = a * d - b * b
    disc = tr * tr - 4 * det
    s = math.sqrt(disc)
    return ((tr + s) / 2, (tr - s) / 2)

def rho_from_eigs(lam_max: float, lam_min: float) -> float:
    if lam_min == 0:
        return float("inf")
    return math.sqrt(lam_max / lam_min)

def rho_a_from_eigs(lam_max: float, lam_min: float) -> float:
    if lam_min == 0:
        return float("inf")
    rho_eig = math.sqrt(lam_max / lam_min)
    return rho_eig - 1.0 / rho_eig

def rho_to_R(rho: float) -> float:
    return (1 - rho) / (1 + rho)

def cone_chain() -> dict:
    V = ROTATIONS
    W = [sub2(V[(i + 1) % 4], V[i]) for i in range(4)]
    G = gram2(V[0], W[0])
    lam1, lam2 = eigvals_sym2(G)
    lam_max = max(lam1, lam2)
    lam_min = min(lam1, lam2)
    rho_eig  = rho_from_eigs(lam_max, lam_min)
    rho_a    = rho_a_from_eigs(lam_max, lam_min)
    R_local  = rho_to_R(rho_a)
    C_local  = 2 + 2 * R_local
    a, b = G[0]; _, d = G[1]
    tr = a + d; det = a * d - b * b
    disc = tr * tr - 4 * det
    return {
        "V": V, "W": W, "G": G,
        "trace": tr, "det": det, "disc": disc,
        "lam": (lam_max, lam_min),
        "rho_eig": rho_eig,
        "rho_a":   rho_a,
        "R":       R_local,
        "C_CONE":  C_local,
        "C_CONE_Q": round(BASE * C_local),
    }


# ─────────────────────────────────────────
# 4. 接続設計
# ─────────────────────────────────────────
@dataclass
class ConeFFTResult:
    projected:        list[int]
    T0:               list[int]
    U0:               list[int]
    T:                list[int]
    U:                list[int]
    merged_parseval:  list[int]
    merged_cone:      list[int]
    parseval_norm:    int
    global_corr:      float
    B:                int
    cone_info:        dict

def apply_fixedpoint_scale(x: list[int], scale_q: int, qbase: int = BASE) -> list[int]:
    out = []
    half = qbase // 2
    for v in x:
        num = v * scale_q
        if num >= 0:
            out.append((num + half) // qbase)
        else:
            out.append(-((-num + half) // qbase))
    return out

def cone_fft_transform(signal: list[int], B: int) -> ConeFFTResult:
    """
    pipeline:
      signal
        -> int_cos_project_3120
        -> split_TU
        -> 40再帰 B段
        -> merge
        -> cone補正 (fixed-point, C_CONE_Q / BASE)
    """
    projected       = int_cos_project_3120(signal)
    T0, U0          = split_TU(projected)
    T, U            = build_fractal_integer_signal(B, T0, U0)
    merged_parseval = merge_TU(T, U)
    info            = cone_chain()
    merged_cone     = apply_fixedpoint_scale(merged_parseval, info["C_CONE_Q"], BASE)
    return ConeFFTResult(
        projected=projected, T0=T0, U0=U0, T=T, U=U,
        merged_parseval=merged_parseval,
        merged_cone=merged_cone,
        parseval_norm=norm2(T) + norm2(U),
        global_corr=r_global(T, U),
        B=B,
        cone_info=info,
    )


# ─────────────────────────────────────────
# 5. デモ
# ─────────────────────────────────────────
def demo() -> None:
    print("─" * 60)
    print("cone-FFT v28 接続完成デモ")
    print("─" * 60)

    seed   = [x + y for x, y in ORBIT]   # [21, 15, 7, 13]
    signal = [seed[i % 4] for i in range(16)]

    print("\n[入力信号]")
    print(signal)

    res = cone_fft_transform(signal, B=2)

    print("\n[cone局所鎖]")
    print("G        =", res.cone_info["G"])
    print("disc     =", res.cone_info["disc"])
    print("rho_eig  =", f'{res.cone_info["rho_eig"]:.8f}')
    print("rho_a    =", f'{res.cone_info["rho_a"]:.8f}')
    print("R        =", f'{res.cone_info["R"]:.8f}')
    print("C_CONE   =", f'{res.cone_info["C_CONE"]:.8f}')
    print("C_CONE_Q =", res.cone_info["C_CONE_Q"], f"(BASE={BASE})")

    print("\n[射影 → Parseval本体]")
    print("len(projected)   =", len(res.projected))
    print("len(T0), len(U0) =", len(res.T0), len(res.U0))
    print("B                =", res.B)
    print("parseval_norm    =", res.parseval_norm)
    print("r_global         =", f"{res.global_corr:.8f}")

    print("\n[出力サンプル]")
    print("merged_parseval[:12] =", res.merged_parseval[:12])
    print("merged_cone[:12]     =", res.merged_cone[:12])

    print("\n[意味]")
    print("  merged_parseval : 40^B を担う整数本体")
    print("  merged_cone     : その本体に 2/φ を固定小数点整数で掛けた cone 補正版")

    # cone補正の比率確認
    ratio = res.cone_info["C_CONE_Q"] / BASE
    print(f"\n[cone補正係数]")
    print(f"  C_CONE_Q/BASE = {ratio:.8f}  (2/φ = {2/PHI:.8f})")

    print("\n" + "─" * 60)
    print("完了")
    print("─" * 60)


if __name__ == "__main__":
    demo()
  
