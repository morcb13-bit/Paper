"""
Local cone invariant chain for cone-FFT v27.

From the local V-W Gram matrix, derive:
    disc -> sqrt(5) -> rho_a -> R = -phi^-2 -> C_CONE = 2/phi.

Proposition B:
    V = (3,4) rotations, W = successive differences.
    Local Gram G = ((25,-25),(-25,50)).
    disc = 3125 = 625×5,
    rho_eig = phi²,
    rho_a   = rho_eig - 1/rho_eig = sqrt(5),
    R       = (1 - rho_a)/(1 + rho_a) = -phi^-2,
    C_CONE  = 2 + 2R = 2/phi.
"""
import math

PHI = (1 + 5**0.5) / 2


# ── 2次元整数ベクトル演算 ──────────────────────

def dot2(a: tuple[int, int], b: tuple[int, int]) -> int:
    return a[0]*b[0] + a[1]*b[1]

def norm2_2(a: tuple[int, int]) -> int:
    return dot2(a, a)

def sub2(a: tuple[int, int], b: tuple[int, int]) -> tuple[int, int]:
    return (a[0]-b[0], a[1]-b[1])

def corr2(a: tuple[int, int], b: tuple[int, int]) -> float:
    den = math.sqrt(norm2_2(a) * norm2_2(b))
    if den == 0:
        return 0.0
    return dot2(a, b) / den


# ── 局所 Gram ──────────────────────────────────

def gram2(
    a: tuple[int, int],
    b: tuple[int, int],
) -> tuple[tuple[int, int], tuple[int, int]]:
    """2ベクトルの 2×2 Gram 行列"""
    return ((dot2(a,a), dot2(a,b)),
            (dot2(a,b), dot2(b,b)))

def eigvals_sym2(
    G: tuple[tuple[int, int], tuple[int, int]],
) -> tuple[float, float]:
    """対称 2×2 行列の固有値 (大, 小)"""
    a, b = G[0]
    _, d = G[1]
    tr   = a + d
    det  = a*d - b*b
    disc = tr*tr - 4*det
    s    = math.sqrt(disc)
    return ((tr+s)/2, (tr-s)/2)


# ── ρ_a → R の変換鎖 ──────────────────────────

def rho_from_eigs(lam_max: float, lam_min: float) -> float:
    """ρ_eig = √(λ+/λ-)"""
    if lam_min == 0:
        return float("inf")
    return math.sqrt(lam_max / lam_min)

def rho_a_from_eigs(lam_max: float, lam_min: float) -> float:
    """ρ_a = ρ_eig - 1/ρ_eig  (反対称化)"""
    if lam_min == 0:
        return float("inf")
    rho_eig = rho_from_eigs(lam_max, lam_min)
    return rho_eig - 1.0/rho_eig

def rho_to_R(rho: float) -> float:
    """R = (1-ρ)/(1+ρ)"""
    return (1 - rho) / (1 + rho)


# ── 完全鎖: V,W → C_CONE ──────────────────────

def cone_chain(
    V: list[tuple[int, int]],
    W: list[tuple[int, int]],
) -> dict:
    """
    4点閉軌道の偏差ベクトル V と差分 W から
    C_CONE = 2/φ までの変換鎖を返す。
    """
    G       = gram2(V[0], W[0])
    lam1, lam2 = eigvals_sym2(G)
    lam_max = max(lam1, lam2)
    lam_min = min(lam1, lam2)

    a, b    = G[0]; _, d = G[1]
    tr      = a + d
    det     = a*d - b*b
    disc    = tr*tr - 4*det

    rho_eig = rho_from_eigs(lam_max, lam_min)
    rho_a   = rho_a_from_eigs(lam_max, lam_min)
    R       = rho_to_R(rho_a)
    C_CONE  = 2 + 2*R

    return {
        "G":       G,
        "disc":    disc,
        "lam":     (lam_max, lam_min),
        "rho_eig": rho_eig,
        "rho_a":   rho_a,
        "R":       R,
        "C_CONE":  C_CONE,
    }
  
