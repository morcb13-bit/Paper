"""B13 第n層の生成。b13_two_tilings.py の Figure を層番号つきに一般化する。
規則は既存のまま。加えたのは (1)層番号 (2)起源の記録 (3)殻番号 (4)近傍索引による高速化。"""
from __future__ import annotations
from collections import deque
from typing import Iterable
from b13_two_tilings import (Zeta, ZERO, ONE, PHI, PHI2, PHI3, zt, zadd, zsub, zmul,
                             zrot, zpow, zconj, norm2, phi_cmp, is_skip, is_cont,
                             CONT_VECS, ring_cells, Figure, to_xy, N_CELL)

class LayeredFigure(Figure):
    def __init__(self) -> None:
        super().__init__()
        self.ring_level: dict[Zeta, int] = {}
        self.ring_origin: dict[Zeta, set[str]] = {}
        self._grid: dict[tuple[int,int], list[Zeta]] = {}   # 近傍索引

    # ── 重なり判定を近傍索引で O(1) 近くに ──
    def _key(self, q: Zeta) -> tuple[int,int]:
        p = to_xy(q); return (int(p.real // 2), int(p.imag // 2))
    def _near(self, q: Zeta):
        kx, ky = self._key(q)
        for dx in (-1,0,1):
            for dy in (-1,0,1):
                yield from self._grid.get((kx+dx, ky+dy), ())
    def _fits(self, new: dict[Zeta, int]) -> bool:
        for q, a in new.items():
            if q in self.cells:
                if self.orient(self.cells[q]) != self.orient(a): return False
                continue
            for p in self._near(q):
                if p != q and phi_cmp(norm2(zsub(q, p)), N_CELL) < 0: return False
            for q2 in new:
                if q2 != q and phi_cmp(norm2(zsub(q, q2)), N_CELL) < 0: return False
        return True
    def add_ring(self, c: Zeta) -> bool:
        if c in self.ring_set: return False
        new = ring_cells(c)
        if not self._fits(new): return False
        for q, a in new.items():
            if q not in self.cells:
                self.cells[q] = a
                self._grid.setdefault(self._key(q), []).append(q)
        self.rings.append(c); self.ring_set.add(c)
        return True

    # ── 層 ──
    def _tag(self, c: Zeta, level: int, origin: str) -> None:
        self.ring_level.setdefault(c, level)
        self.ring_origin.setdefault(c, set()).add(origin)

    def pentagram_candidates(self) -> set[Zeta]:
        out = set()
        for g, k2 in self.pentagrams():
            for i in range(5):
                c = zadd(g, zrot(PHI3, k2 + 2*i))
                if c not in self.ring_set: out.add(c)
        return out

    def rhombus_candidates(self) -> set[Zeta]:
        out = set()
        for a, b in self.skip_pairs():
            for v in CONT_VECS:
                c = zadd(a, v)
                if c not in self.ring_set and is_cont(c, b): out.add(c)
        return out

    def grow_one_layer(self, level: int) -> dict[str, set[Zeta]]:
        P = self.pentagram_candidates()
        R = self.rhombus_candidates()
        added_p, added_r = set(), set()
        for c in sorted(P | R):
            if self.add_ring(c):
                if c in P: self._tag(c, level, "pentagram"); added_p.add(c)
                if c in R: self._tag(c, level, "rhombus");   added_r.add(c)
        return {"pentagram_only": added_p - added_r,
                "rhombus_only":   added_r - added_p,
                "both":           added_p & added_r,
                "all":            added_p | added_r}

    def grow_to_layer(self, n: int) -> list[dict]:
        log = []
        for lv in range(1, n + 1):
            a = self.grow_one_layer(lv)
            log.append({"層": lv, "追加": len(a["all"]),
                        "五芒星のみ": len(a["pentagram_only"]),
                        "ひし形のみ": len(a["rhombus_only"]),
                        "両方": len(a["both"]),
                        "累計円環": len(self.rings), "累計五角形": len(self.cells),
                        "五芒星": len(self.pentagrams())})
            if not a["all"]: break
        return log

    def rings_at_layer(self, n: int) -> set[Zeta]:
        return {c for c, lv in self.ring_level.items() if lv == n}

    # ── 殻番号（生成世代とは別。完成配置だけで一意に決まる） ──
    def shell_levels(self, seeds: Iterable[Zeta]) -> dict[Zeta, int]:
        lev = {c: 0 for c in seeds if c in self.ring_set}
        q = deque(lev)
        while q:
            a = q.popleft()
            for b in self.rings:
                if b in lev: continue
                if is_cont(a, b) or is_skip(a, b):
                    lev[b] = lev[a] + 1; q.append(b)
        return lev

    def cells_at_layer(self, n: int) -> dict[Zeta, int]:
        """第n層で初めて現れた五角形だけ"""
        cur, earlier = {}, set()
        for c in self.rings_at_layer(n): cur.update(ring_cells(c))
        for lv in range(n):
            for c in self.rings_at_layer(lv): earlier.update(ring_cells(c))
        return {q: a for q, a in cur.items() if q not in earlier}


def initial_rings() -> list[Zeta]:
    """原点の五芒星を囲む5円環＝第0層"""
    return [zrot(PHI3, 2*i) for i in range(5)]

def build_layers(n: int) -> tuple[LayeredFigure, list[dict]]:
    F = LayeredFigure()
    seed = initial_rings()
    assert F.add_group(seed), "初期層を配置できない"
    for c in seed: F._tag(c, 0, "seed")
    return F, F.grow_to_layer(n)


if __name__ == "__main__":
    F, log = build_layers(8)
    print(f"{'層':>3} {'追加':>5} {'五芒星のみ':>9} {'ひし形のみ':>9} {'両方':>5} "
          f"{'累計円環':>7} {'累計五角形':>9} {'五芒星':>6}")
    for r in log:
        print(f"{r['層']:>3} {r['追加']:>5} {r['五芒星のみ']:>9} {r['ひし形のみ']:>9} "
              f"{r['両方']:>5} {r['累計円環']:>7} {r['累計五角形']:>9} {r['五芒星']:>6}")
    print("\n── 既知の確定値との照合 ──")
    print(f"  55環385枚 が途中に現れるか: "
          f"{any(r['累計円環']==55 and r['累計五角形']==385 for r in log)}")
    print(f"  60環420枚 が途中に現れるか: "
          f"{any(r['累計円環']==60 and r['累計五角形']==420 for r in log)}")
    print(f"  五芒星21個 が現れるか: {any(r['五芒星']==21 for r in log)}")
    sh = F.shell_levels(initial_rings())
    import collections
    print(f"\n  殻ごとの円環数: {dict(sorted(collections.Counter(sh.values()).items()))}")
    print(f"  生成世代ごと:   {dict(sorted(collections.Counter(F.ring_level.values()).items()))}")
    d = sum(1 for c in F.ring_set if sh.get(c, -1) != F.ring_level.get(c, -2))
    print(f"  世代≠殻 の円環: {d} / {len(F.ring_set)}  ← この差が残渣")


# ── インフレーション：加算が飽和した先を φ³ の相似で継ぐ ──────────────
def inflate(rings, n: int = 1):
    """円環中心を φ^{3n} 倍する。Z[ζ₁₀] の整数演算のみ、誤差なし。"""
    s = zpow(PHI3, n)
    return [zmul(s, c) for c in rings]

def layer_by_inflation(rings, n: int) -> dict:
    """第n周の骨格。積み上げ不要で直接出る。"""
    out = inflate(rings, n)
    r2 = [norm2(c) for c in out if c != ZERO]
    return {"周": n, "円環": len(out),
            "半径²の最大": max(r2, key=lambda x: x[0] + 2*x[1]) if r2 else None,
            "中心": out}
