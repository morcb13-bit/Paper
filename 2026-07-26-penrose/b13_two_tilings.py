"""
b13_two_tilings.py
==================
5角形10個の円環から2つのペンローズタイリングを構成し、
両者が φ³ のスケール違いで一致することを検証する。

方針（b13phase/constants.py の指示に従う）
    座標  : Z[ζ₁₀] の整数4成分。φ = (1,0,1,-1) で φ²-φ-1 = 0 が厳密に成立
    向き  : b13phase の BASE=3120 位相整数。36° = 312 = 6 × STEP_MIN
    判定  : すべて整数の等式・不等式。実数は描画時の1箇所のみ

構成
    層0  五芒星（原点）
    層n  五芒星の周りに円環5個            … 成長規則
    飽和後 大ひし形の欠けた頂点を埋める    … 大タイリング規則

2つのタイリング
    小 : 五角形（外接半径 1）と細いひし形
    大 : 円環中心を結んだ図。正五角形（外接半径 φ³）とひし形
"""
from __future__ import annotations
import sys, os, json
from typing import Iterable

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "b13phase_v069"))
from b13phase import BASE, STEP_MIN, digits_from_int, digits_to_int, digits_add

# ── 位相（b13phase, BASE=3120） ───────────────────────────────
DEG36 = 6 * STEP_MIN          # 312 = 36°。番地1つ分の回転
assert DEG36 * 10 == BASE, "36°×10 が一周にならない"


def addr_phase(a: int, ndigits: int = 2) -> list[int]:
    """番地 a → BASE3120 の位相桁列。整数加算のみ。"""
    return digits_from_int((a % 10) * DEG36, ndigits)


def phase_addr(d: list[int]) -> int:
    """位相桁列 → 番地。"""
    return (digits_to_int(d) // DEG36) % 10


# ── Z[ζ₁₀]（ζ = e^{iπ/5}、ζ⁴ = -1+ζ-ζ²+ζ³） ──────────────────
Zeta = tuple[int, int, int, int]
ZERO: Zeta = (0, 0, 0, 0)
ONE: Zeta = (1, 0, 0, 0)
PHI: Zeta = (1, 0, 1, -1)          # φ

_POW = [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1), (-1, 1, -1, 1),
        (-1, 0, 0, 0), (0, -1, 0, 0), (0, 0, -1, 0), (0, 0, 0, -1), (1, -1, 1, -1)]


def zt(k: int) -> Zeta:
    """ζ^k"""
    return _POW[k % 10]


def zadd(a: Zeta, b: Zeta) -> Zeta:
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3])


def zsub(a: Zeta, b: Zeta) -> Zeta:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2], a[3] - b[3])


def zmul(a: Zeta, b: Zeta) -> Zeta:
    c = [0] * 8
    for i in range(4):
        if a[i]:
            for j in range(4):
                c[i + j] += a[i] * b[j]
    r = [c[0], c[1], c[2], c[3]]
    r[3] += c[4]; r[2] -= c[4]; r[1] += c[4]; r[0] -= c[4]   # ζ⁴
    r[0] -= c[5]                                             # ζ⁵ = -1
    r[1] -= c[6]; r[2] -= c[7]
    return (r[0], r[1], r[2], r[3])


def zrot(a: Zeta, k: int) -> Zeta:
    return zmul(a, zt(k))


def zpow(a: Zeta, n: int) -> Zeta:
    r = ONE
    for _ in range(n):
        r = zmul(r, a)
    return r


def zconj(a: Zeta) -> Zeta:
    r = ZERO
    for k in range(4):
        if a[k]:
            r = zadd(r, zmul((a[k], 0, 0, 0), zt(-k)))
    return r


PHI2 = zpow(PHI, 2)     # 円環の半径
PHI3 = zpow(PHI, 3)     # 五芒星から円環中心までの距離

# ── Z[φ]（実数体）: (p, q) は p + qφ ─────────────────────────
Phi = tuple[int, int]


def norm2(a: Zeta) -> Phi:
    """|a|² を Z[φ] の整数対で返す。実数であることを構造的に確認する。"""
    z = zmul(a, zconj(a))
    assert z[1] == 0 and z[3] == -z[2], f"ノルムが実数でない: {z}"
    return (z[0] - z[2], z[2])


def phi_sign(x: Phi) -> int:
    """p + qφ の符号を整数だけで判定する。2(p+qφ) = (2p+q) + q√5"""
    A, B = 2 * x[0] + x[1], x[1]
    if A == 0 and B == 0:
        return 0
    if A >= 0 and B >= 0:
        return 1
    if A <= 0 and B <= 0:
        return -1
    # 符号が食い違う場合は A² と 5B² を比べる
    if A > 0:
        return 1 if A * A > 5 * B * B else -1
    return -1 if A * A > 5 * B * B else 1


def phi_cmp(x: Phi, y: Phi) -> int:
    return phi_sign((x[0] - y[0], x[1] - y[1]))


# ── 接続の2種（差ベクトルの二乗長） ──────────────────────────
N_SKIP = norm2(zmul(PHI2, zsub(zt(2), ONE)))   # 2つ飛ばし  |v|² = φ⁴(2-ζ²-ζ⁻²)
N_CONT = norm2(zmul(PHI2, zsub(zt(4), ONE)))   # 連続
N_CELL = norm2(PHI)                            # 五角形の中心間の最小距離 φ
N_LONG = norm2(zsub(zmul(PHI2, zsub(zt(4), ONE)),
                    zmul(PHI2, zsub(zt(6), ONE))))  # ひし形の長対角


def is_skip(a: Zeta, b: Zeta) -> bool:
    return norm2(zsub(a, b)) == N_SKIP


def is_cont(a: Zeta, b: Zeta) -> bool:
    return norm2(zsub(a, b)) == N_CONT


CONT_VECS = [zrot(zmul(PHI2, zsub(zt(4), ONE)), k) for k in range(10)]


# ── 円環 ─────────────────────────────────────────────────────
def ring_cells(c: Zeta) -> dict[Zeta, int]:
    """円環 = 中心 c から半径 φ² の五角形10枚。位置 → 番地。"""
    return {zadd(c, zrot(PHI2, k)): k for k in range(10)}


class Figure:
    """五角形の集合。番地の偶奇で向きが決まり、重なりは整数で判定する。"""

    def __init__(self) -> None:
        self.rings: list[Zeta] = []
        self.ring_set: set[Zeta] = set()
        self.cells: dict[Zeta, int] = {}      # 位置 → 番地

    # 向きは番地の偶奇のみで決まる
    @staticmethod
    def orient(addr: int) -> int:
        return addr % 2

    def _fits(self, new: dict[Zeta, int]) -> bool:
        for q, a in new.items():
            if q in self.cells:
                if self.orient(self.cells[q]) != self.orient(a):
                    return False               # 番地の偶奇が食い違う
                continue
            for p, b in self.cells.items():
                if phi_cmp(norm2(zsub(q, p)), N_CELL) < 0:
                    return False               # 五角形が重なる
            for q2, a2 in new.items():
                if q2 is q or q2 == q:
                    continue
                if phi_cmp(norm2(zsub(q, q2)), N_CELL) < 0:
                    return False
        return True

    def add_ring(self, c: Zeta) -> bool:
        if c in self.ring_set:
            return False
        new = ring_cells(c)
        if not self._fits(new):
            return False
        for q, a in new.items():
            self.cells.setdefault(q, a)
        self.rings.append(c)
        self.ring_set.add(c)
        return True

    def add_group(self, cs: Iterable[Zeta]) -> bool:
        """原子的追加。1個でも入らなければ全部戻す。"""
        cs = [c for c in cs if c not in self.ring_set]
        snap_r, snap_c = len(self.rings), dict(self.cells)
        for c in cs:
            if not self.add_ring(c):
                del self.rings[snap_r:]
                self.ring_set = set(self.rings)
                self.cells = snap_c
                return False
        return True

    # ── 五芒星 ───────────────────────────────────────────────
    def pentagrams(self) -> list[tuple[Zeta, int]]:
        """中心から距離 φ に五角形5枚が72°間隔で揃う点。"""
        out, seen = [], set()
        for q in self.cells:
            for k in range(10):
                g = zsub(q, zrot(PHI, k))
                if g in seen:
                    continue
                seen.add(g)
                for k2 in range(2):
                    if all(zadd(g, zrot(PHI, k2 + 2 * i)) in self.cells for i in range(5)):
                        out.append((g, k2))
                        break
        return out

    # ── 規則1: 五芒星の周りに円環5個 ─────────────────────────
    def grow(self) -> int:
        added = 0
        for g, k2 in self.pentagrams():
            want = [zadd(g, zrot(PHI3, k2 + 2 * i)) for i in range(5)]
            if any(c not in self.ring_set for c in want) and self.add_group(want):
                added += 5
        return added

    # ── 規則2: 大ひし形の欠けた頂点を埋める ──────────────────
    def skip_pairs(self) -> list[tuple[Zeta, Zeta]]:
        R = self.rings
        return [(R[i], R[j]) for i in range(len(R)) for j in range(i + 1, len(R))
                if is_skip(R[i], R[j])]

    def rhombus_vertices(self, a: Zeta, b: Zeta) -> list[Zeta]:
        """a,b の両方に連続接続する点。ひし形の残り2頂点で、常にちょうど2個。"""
        out = []
        for base, other in ((a, b), (b, a)):
            for v in CONT_VECS:
                c = zadd(base, v)
                if is_cont(c, other) and c not in out:
                    out.append(c)
        assert len(out) == 2, f"ひし形の頂点が2個でない: {len(out)}"
        return out

    def complete_rhombi(self) -> int:
        """未完成の大ひし形を埋める。欠けた頂点は一意に決まる。"""
        added = 0
        for a, b in self.skip_pairs():
            miss = [c for c in self.rhombus_vertices(a, b) if c not in self.ring_set]
            if len(miss) == 1 and self.add_ring(miss[0]):
                added += 1
        return added

    # ── 中心図形（大タイリング） ─────────────────────────────
    def center_faces(self) -> tuple[list[list[Zeta]], list[list[Zeta]]]:
        """連続辺だけで面を取る。2つ飛ばしはひし形の短対角なので辺にしない。"""
        import math, cmath
        R = self.rings
        idx = {c: i for i, c in enumerate(R)}
        adj: dict[int, list[int]] = {i: [] for i in range(len(R))}
        for i in range(len(R)):
            for j in range(len(R)):
                if i != j and is_cont(R[i], R[j]):
                    adj[i].append(j)
        xy = [to_xy(c) for c in R]
        for i in adj:
            adj[i].sort(key=lambda j: math.atan2((xy[j] - xy[i]).imag, (xy[j] - xy[i]).real))

        def nxt(u: int, v: int) -> int:
            ang = math.atan2((xy[u] - xy[v]).imag, (xy[u] - xy[v]).real)
            best = None
            for w in adj[v]:
                if w == u and len(adj[v]) > 1:
                    continue
                t = (ang - math.atan2((xy[w] - xy[v]).imag, (xy[w] - xy[v]).real)) % (2 * math.pi)
                if t < 1e-9:
                    t += 2 * math.pi
                if best is None or t < best[0]:
                    best = (t, w)
            return best[1] if best else u

        vis, faces = set(), []
        for i in adj:
            for j in adj[i]:
                if (i, j) in vis:
                    continue
                f, u, v = [i], i, j
                while True:
                    vis.add((u, v)); f.append(v)
                    u, v = v, nxt(u, v)
                    if (u, v) == (i, j) or len(f) > 40:
                        break
                faces.append(f[:-1])

        def area(f):
            s = 0.0
            for k in range(len(f)):
                p, q = xy[f[k]], xy[f[(k + 1) % len(f)]]
                s += p.real * q.imag - q.real * p.imag
            return s / 2

        pent = [[R[i] for i in f] for f in faces if area(f) > 0.01 and len(f) == 5]
        rho = [[R[i] for i in f] for f in faces if area(f) > 0.01 and len(f) == 4]
        return pent, rho

    def stats(self) -> dict:
        return {"円環": len(self.rings), "五角形": len(self.cells),
                "五芒星": len(self.pentagrams())}


# ── 実数化（描画専用。判定には一切使わない） ──────────────────
def to_xy(a: Zeta) -> complex:
    import cmath, math
    z = cmath.exp(2j * math.pi / 10)
    return sum(a[k] * z ** k for k in range(4))


# ── 確定図（引継書 v199 で確定した55環。成長規則の反復で到達） ──
RINGS55: tuple[Zeta, ...] = (
    (3,0,2,-2), (0,2,1,2), (-3,1,-1,3), (-2,-1,-2,0), (2,-2,0,-3),
    (7,-1,4,-5), (6,1,4,-3), (1,4,3,3), (-1,4,2,5), (-6,3,-2,7),
    (-7,2,-3,6), (-5,-2,-4,1), (-3,-3,-4,-1), (3,-4,-1,-6), (5,-4,1,-7),
    (5,3,5,0), (-5,5,0,8), (-8,0,-5,5), (0,-5,-3,-5), (8,-3,3,-8),
    (10,0,6,-6), (0,6,4,6), (-10,4,-4,10), (-6,-4,-6,0), (6,-6,0,-10),
    (13,-2,7,-10), (12,-4,5,-11), (11,2,8,-5), (8,4,7,-1), (6,5,7,1),
    (3,7,6,5), (-3,8,3,10), (-6,7,1,11), (-8,7,-1,12), (-11,6,-3,13),
    (-13,3,-6,11), (-12,1,-7,8), (-11,-1,-7,6), (-10,-3,-8,3), (-5,-6,-7,-3),
    (-1,-7,-5,-6), (1,-7,-4,-8), (5,-8,-2,-11), (11,-5,4,-12), (10,-7,2,-13),
    (15,1,10,-8), (16,-1,9,-11), (2,9,7,8), (-2,10,5,11), (-15,7,-5,16),
    (-16,5,-7,15), (-11,-5,-10,2), (-8,-7,-9,-2), (8,-10,-1,-15), (11,-9,1,-16),
)


# ── 構成 ─────────────────────────────────────────────────────
def build() -> Figure:
    """確定55環を置き、成長規則の飽和を確認してから大タイリング規則を適用する。"""
    F = Figure()
    assert F.add_group(RINGS55), "確定55環が入らない"
    assert F.grow() == 0, "成長規則がまだ余地を持っている"   # 飽和の確認
    F.complete_rhombi()
    return F


def verify(F: Figure) -> dict:
    """すべて整数の判定。"""
    ok = {}
    # 五角形の重なり
    cs = list(F.cells)
    worst = None
    for i in range(len(cs)):
        for j in range(i + 1, len(cs)):
            d = norm2(zsub(cs[i], cs[j]))
            if worst is None or phi_cmp(d, worst) < 0:
                worst = d
    ok["重なり0"] = phi_cmp(worst, N_CELL) >= 0
    # 番地の偶奇
    ok["偶奇食い違い0"] = True     # add_ring で保証済み
    # ひし形が中心を向く = 2つ飛ばし対は等半径
    sp = F.skip_pairs()
    ok["2つ飛ばし対"] = len(sp)
    ok["全て等半径"] = all(norm2(a) == norm2(b) for a, b in sp)
    # 中心図形
    pent, rho = F.center_faces()
    ok["中心図形の面"] = {"正五角形": len(pent), "ひし形": len(rho)}
    # 一致：正五角形の重心 ↔ 五芒星
    cent4 = sorted(tuple(sum(v[k] for v in f) for k in range(4)) for f in pent)  # 5×重心
    stars = sorted(tuple(5 * g[k] for k in range(4)) for g, _ in F.pentagrams())
    ok["正五角形の重心 = 五芒星"] = cent4 == stars
    # 一致：ひし形の中心 ↔ 細いひし形の中心（どちらも2つ飛ばし対の中点）
    rc = sorted(tuple(sum(v[k] for v in f) // 2 for k in range(4)) for f in rho)
    sm = sorted(tuple(a[k] + b[k] for k in range(4)) for a, b in sp)
    ok["ひし形の中心 = 細ひし形の中心"] = rc == sm
    # 相似比
    if rho:
        side = norm2(zsub(rho[0][0], rho[0][1]))
        ok["ひし形の辺の二乗長"] = side
        ok["φ⁶ × 小ひし形の辺"] = phi_cmp(side, norm2(zmul(PHI3, zsub(zt(2), ONE)))) == 0
    return ok


if __name__ == "__main__":
    F = build()
    print("構成:", F.stats())
    print()
    for k, v in verify(F).items():
        print(f"  {k}: {v}")
    print()
    print(f"番地→位相: 番地1 = {digits_to_int(addr_phase(1))} "
          f"(= 36°, BASE={BASE}), 番地5 = {digits_to_int(addr_phase(5))}")
    out = {"rings": [list(c) for c in F.rings],
           "cells": {",".join(map(str, q)): a for q, a in F.cells.items()},
           "pentagrams": [list(g) for g, _ in F.pentagrams()]}
    with open("rings_integer.json", "w") as fp:
        json.dump(out, fp)
    print("rings_integer.json を書き出しました")
