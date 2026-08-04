"""
pentagon_processor_base.py
==========================
ペンタゴンプロセッサ シミュレーター — 基本図形（第二層・円環20個）

すべての座標・向き・接続を **整数のみ** で保持する。
浮動小数点は SVG 描画（to_xy）の1箇所でしか使わない。

座標系
------
Z[ζ]、ζ = exp(2πi/10)。円分多項式 Φ₁₀(x) = x⁴ - x³ + x² - x + 1 より

    ζ⁴ = -1 + ζ - ζ² + ζ³

基底 {1, ζ, ζ², ζ³} の整数4成分 (c0, c1, c2, c3) で1点を表す。
黄金比は整数で書ける：

    φ = 1 + ζ² - ζ³
    φⁿ = (F(n+1), 0, F(n), -F(n))        F はフィボナッチ数

    φ¹=(1,0,1,-1)  φ²=(2,0,1,-1)  φ³=(3,0,2,-2)
    φ⁴=(5,0,3,-3)  φ⁵=(8,0,5,-5)  φ⁶=(13,0,8,-8)

角度（b13phase 連携）
--------------------
b13phase の BASE = 3120 を 360° とする位相インデックスに一致する。

    36° = BASE/10  = 312 = 6 × STEP_MIN(52)
    72° = BASE/5   = 624 = STEP_Z5

五角形の「番地」a は ζ^a の指数であり、位相インデックスは a × 312。
番地は mod 10 で意味を持ち、偶奇だけが向きを決める（a と a+2 は同一の五角形）。

構成則
------
    環半径          R  = φ²
    2つ飛ばし接続   φ²·ζ^a·(ζ² - 1)     |·| = φ²·2sin36°
    連続接続        φ²·ζ^a·(ζ⁴ - 1)     |·| = φ³·2sin36°

    第一層 L1  (5個)  : φ³·ζ^(2m)
    軌道   ORB (10個) : L1_m + φ²·ζ^(2m+5)(ζ⁴-1) , L1_m + φ²·ζ^(2m+6)(ζ⁴-1)
    充填   NEW (5個)  : φ⁵·ζ^(2m+1)
"""

from __future__ import annotations
from typing import List, Tuple, Dict

Zeta = Tuple[int, int, int, int]

# ── b13phase 位相グリッド ────────────────────────────────
try:                                     # ライブラリがあれば実値を使う
    from b13phase import BASE, STEP_MIN, STEP_Z5
except Exception:                        # なければ同値を内蔵
    BASE, STEP_MIN, STEP_Z5 = 3120, 52, 624

STEP_Z10 = BASE // 10                    # 312 → 36°
assert STEP_Z10 == 6 * STEP_MIN and STEP_Z5 == 2 * STEP_Z10


# ── Z[ζ] 整数演算 ───────────────────────────────────────
ONE:  Zeta = (1, 0, 0, 0)
ZERO: Zeta = (0, 0, 0, 0)
PHI:  Zeta = (1, 0, 1, -1)


def zmulz(a: Zeta, k: int = 1) -> Zeta:
    """a に ζ^k を掛ける（ζ⁴ を基底に畳む）。"""
    c0, c1, c2, c3 = a
    for _ in range(k % 10):
        c0, c1, c2, c3 = -c3, c0 + c3, c1 - c3, c2 + c3
    return (c0, c1, c2, c3)


def zadd(a: Zeta, b: Zeta) -> Zeta:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3])


def zsub(a: Zeta, b: Zeta) -> Zeta:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2], a[3] - b[3])


def zmul(a: Zeta, b: Zeta) -> Zeta:
    r = ZERO
    for i, c in enumerate(a):
        if c:
            t = zmulz(b, i)
            r = zadd(r, (c * t[0], c * t[1], c * t[2], c * t[3]))
    return r


def zpow(a: Zeta, n: int) -> Zeta:
    r = ONE
    for _ in range(n):
        r = zmul(r, a)
    return r


def zeta_pow(k: int) -> Zeta:
    """ζ^k を整数4成分で返す。"""
    return zmulz(ONE, k)


def phase_index(address: int) -> int:
    """番地 → b13phase 位相インデックス（BASE=3120、36°=312）。"""
    return (address % 10) * STEP_Z10


# ── 構成定数（すべて整数） ───────────────────────────────
R_RING   = zpow(PHI, 2)                                  # (2, 0, 1, -1)
V_SKIP2  = zmul(R_RING, zsub(zeta_pow(2), ONE))          # 2つ飛ばし
V_CONT   = zmul(R_RING, zsub(zeta_pow(4), ONE))          # 連続


# ── 基本図形の生成 ──────────────────────────────────────
def ring_centers() -> List[Tuple[str, int, Zeta]]:
    """円環20個の中心。(層タグ, 5回対称の枝番, Z[ζ]座標)"""
    out: List[Tuple[str, int, Zeta]] = []
    for m in range(5):
        out.append(("L1", m, zmul(zpow(PHI, 3), zeta_pow(2 * m))))
    for m in range(5):
        base = zmul(zpow(PHI, 3), zeta_pow(2 * m))
        for a in (2 * m + 5, 2 * m + 6):
            out.append(("ORB", m, zadd(base, zmul(V_CONT, zeta_pow(a)))))
    for m in range(5):
        out.append(("NEW", m, zmul(zpow(PHI, 5), zeta_pow(2 * m + 1))))
    return out


def pentagons() -> Dict[Zeta, int]:
    """五角形145枚。 中心 Z[ζ] → 番地（0..9、偶奇のみが向きを決める）"""
    table: Dict[Zeta, int] = {}
    for _, _, c in ring_centers():
        for j in range(10):
            pc = zadd(c, zmul(R_RING, zeta_pow(j)))
            prev = table.get(pc)
            if prev is None:
                table[pc] = j
            elif (prev - j) % 2:
                raise AssertionError(f"番地の偶奇が食い違う: {pc} {prev} vs {j}")
    return table


def pentagon_vertices(center: Zeta, address: int) -> List[Zeta]:
    """五角形の頂点5点。すべて Z[ζ] の整数。"""
    return [zadd(center, zeta_pow(address + 2 * i)) for i in range(5)]


def classify(c1: Zeta, c2: Zeta):
    """2環の接続を整数一致で判定。('2つ飛ばし'|'連続'|None, 向きa)"""
    d = zsub(c2, c1)
    for vec, name in ((V_SKIP2, "2つ飛ばし"), (V_CONT, "連続")):
        for a in range(10):
            if d == zmul(vec, zeta_pow(a)):
                return name, a
    return None, None


def connections() -> List[Tuple[int, int, str, int]]:
    """全接続。(環i, 環j, 種別, 向きa)"""
    rc = ring_centers()
    out = []
    for i in range(len(rc)):
        for j in range(i + 1, len(rc)):
            n, a = classify(rc[i][2], rc[j][2])
            if n:
                out.append((i, j, n, a))
    return out


# ── 検算（整数のみ） ────────────────────────────────────
def verify() -> Dict[str, object]:
    rc = ring_centers()
    pt = pentagons()                      # 偶奇の食い違いはここで例外
    conn = connections()
    kinds = {"2つ飛ばし": 0, "連続": 0}
    for *_, n, _a in conn:
        kinds[n] += 1

    # 五角形の頂点も整数のまま集合として一意か
    vset = set()
    for c, a in pt.items():
        vset.update(pentagon_vertices(c, a))

    # 同一中心に異なる向きが来ていないこと（偶奇一致 ⇒ 頂点集合一致）
    for c, a in pt.items():
        assert set(pentagon_vertices(c, a)) == set(pentagon_vertices(c, a + 2))

    return {
        "円環": len(rc),
        "五角形": len(pt),
        "頂点": len(vset),
        "接続": kinds,
        "番地の偶奇食い違い": 0,
        "浮動小数点の使用": "なし（描画時のみ）",
    }


# ── 描画（ここだけ実数） ────────────────────────────────
def to_xy(a: Zeta) -> Tuple[float, float]:
    import math
    z = [(math.cos(2 * math.pi * k / 10), math.sin(2 * math.pi * k / 10))
         for k in range(4)]
    x = sum(c * z[i][0] for i, c in enumerate(a))
    y = sum(c * z[i][1] for i, c in enumerate(a))
    return x, y


def to_svg(scale: float = 26.0, pad: float = 24.0) -> str:
    pt = pentagons()
    rc = ring_centers()
    col = {"L1": "#4fb3ff", "ORB": "#ffb648", "NEW": "#5be0c0"}

    pts = [to_xy(v) for c, a in pt.items() for v in pentagon_vertices(c, a)]
    xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
    w = (max(xs) - min(xs)) * scale + 2 * pad
    h = (max(ys) - min(ys)) * scale + 2 * pad

    def T(p):
        return (pad + (p[0] - min(xs)) * scale, h - pad - (p[1] - min(ys)) * scale)

    owner = {}
    for tag, _m, c in rc:
        for j in range(10):
            owner.setdefault(zadd(c, zmul(R_RING, zeta_pow(j))), tag)

    body = [f'<rect width="{w:.0f}" height="{h:.0f}" fill="#0e1116"/>']
    for c, a in pt.items():
        k = col[owner[c]]
        d = " ".join(f"{T(to_xy(v))[0]:.2f},{T(to_xy(v))[1]:.2f}"
                     for v in pentagon_vertices(c, a))
        body.append(f'<polygon points="{d}" fill="{k}" fill-opacity="0.16" '
                    f'stroke="{k}" stroke-width="1.1"/>')
    for tag, _m, c in rc:
        x, y = T(to_xy(c))
        body.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="3.2" fill="{col[tag]}"/>')
    return (f'<svg viewBox="0 0 {w:.0f} {h:.0f}" width="{w:.0f}" '
            f'height="{h:.0f}" xmlns="http://www.w3.org/2000/svg">'
            + "".join(body) + "</svg>")


if __name__ == "__main__":
    for k, v in verify().items():
        print(f"{k}: {v}")
    print()
    print("円環20個（Z[ζ] 整数4成分 / b13phase 位相）")
    for i, (tag, m, c) in enumerate(ring_centers()):
        print(f"  {i:2d} {tag:3s} m={m}  {str(c):>18s}")
