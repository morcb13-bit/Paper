#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
penrose_retina.py ── 自分の描いた図を、描画コードを使わずに検める

このファイル一つで完結します。必要なのは numpy / Pillow / scipy(cKDTree) だけです。

────────────────────────────────────────────────────────────
Claude（または任意の LLM）に渡す場合、次の順で走らせてください
────────────────────────────────────────────────────────────
  1) まず自己検査：
         python penrose_retina.py
     引数なしだと、設計から円環を描き、その PNG だけを装置に渡し、
     必ず落ちる設定と一緒に数値を出します。乗る率が 1.000 に近く、
     必ず落ちる設定が落ちていれば、環境は正しく動いています。

  2) 手元の図を検める：
         python penrose_retina.py --image mydrawing.png
     matplotlib で描いた図など、枠・格子線・目盛りが入っている場合は
     先に切り出してください（それらも「点いたセル」を食うため）：
         python penrose_retina.py --image plot.png --box 142,129,1336,1267

  3) 中心円の比を弁別する（辺で接する φ² か、頂点で触れる cot18° か）：
         python penrose_retina.py --image mydrawing.png --scan-ratio

  ※ 出力の「必ず落ちる設定」が落ちていない場合、その検査は判定に
     なっていません。乗る率だけを報告せず、必ず両方を出してください。
────────────────────────────────────────────────────────────

装置がすること
  ① 画面を読んで、どのセルが点いているかを決める（読み＝明るさの変わった量）
  ② 型紙の各点を最も近いセルへ着地させ、そのセルが点いているかを返す
  判断はしません。返るのは点ごとの二値だけです。

検図の条件
  型紙は「設計」から作り、図からは中心・尺度・向きの三つだけを取ります。
  型紙と図を同じコードから作ると、規則が間違っていても乗る率は高く出ます。
  照合ではなく自分の写しを自分と比べているだけになるので、出所を分けます。

設計（円環）
  中心から半径 φ² の位置に、正五角形10枚。向きは番地の偶奇（0° と 36° が交互）。
  このとき隣の中心距離は φ ＝ 内接半径の2倍になり、10枚が辺で接して閉じます。
"""
from __future__ import annotations
import argparse, math, sys
import numpy as np

PHI = (1.0 + 5.0 ** 0.5) / 2.0
DEG = math.pi / 180.0
RATIOS = {"phi2": PHI ** 2, "cot18": 1.0 / math.tan(18 * DEG)}


# ══ 読み ════════════════════════════════════════════════════
def blur(a: np.ndarray, sigma: float) -> np.ndarray:
    """分離型のガウスぼかし（端は複製）。"""
    if sigma <= 0:
        return a
    r = max(1, int(3 * sigma))
    k = np.exp(-0.5 * (np.arange(-r, r + 1) / sigma) ** 2)
    k /= k.sum()
    out = a
    for ax in (0, 1):
        p = np.pad(out, [(r, r) if i == ax else (0, 0) for i in range(2)], mode="edge")
        out = np.apply_along_axis(lambda v: np.convolve(v, k, mode="valid"), ax, p)
    return out


def to_gray(path: str) -> np.ndarray:
    from PIL import Image
    A = np.asarray(Image.open(path).convert("RGB")).astype(float) / 255.0
    return 0.299 * A[:, :, 0] + 0.587 * A[:, :, 1] + 0.114 * A[:, :, 2]


# ══ 担体（切断投影・ドブラウンのペンタグリッド）══════════════
def penrose_points(radius_units: float, offsets: np.ndarray) -> np.ndarray:
    z = np.exp(2j * np.pi * np.arange(5) / 5)
    K = int(radius_units) + 2
    ks = np.arange(-K, K + 1, dtype=float)
    out = []
    for r in range(5):
        for s in range(r + 1, 5):
            kr, ksd = (v.ravel() for v in np.meshgrid(ks, ks, indexing="ij"))
            a1, b1 = z[r].real, z[r].imag
            a2, b2 = z[s].real, z[s].imag
            det = a1 * b2 - a2 * b1
            x = ((kr - offsets[r]) * b2 - (ksd - offsets[s]) * b1) / det
            y = (a1 * (ksd - offsets[s]) - a2 * (kr - offsets[r])) / det
            m = x * x + y * y < radius_units ** 2
            x, y = x[m], y[m]
            if x.size == 0:
                continue
            w = x + 1j * y
            Kj = np.empty((5, x.size))
            for j in range(5):
                Kj[j] = np.ceil((w * np.conj(z[j])).real + offsets[j])
            Kj[r], Kj[s] = kr[m], ksd[m]
            v = Kj.T @ z
            out.append(np.stack([v.real, v.imag], 1))
    return np.unique(np.round(np.concatenate(out, 0), 6), axis=0)


class Retina:
    """担体・点灯・着地。判断は持たない。"""

    def __init__(self, L: np.ndarray, cell: float = 4.0, quantile: float = 0.85,
                 seed: int = 13):
        from scipy.spatial import cKDTree
        H, W = L.shape
        rg = np.random.default_rng(seed)
        off = rg.uniform(-0.4, 0.4, 5)
        off -= off.sum() / 5                       # 一般の位置（特異点を避ける）
        P = penrose_points(0.5 * math.hypot(W, H) / cell + 3, off) * cell
        P[:, 0] += W / 2.0
        P[:, 1] += H / 2.0
        m = (P[:, 0] >= 1) & (P[:, 0] < W - 1) & (P[:, 1] >= 1) & (P[:, 1] < H - 1)
        self.P, self.cell, self.shape = P[m], cell, (H, W)

        # 読み：明るさの変わった量。着地の粗さぶんだけ窓を広げる
        E = blur(np.abs(L - blur(L, 3.0)), 2.0)
        xi = self.P[:, 0].astype(int)
        yi = self.P[:, 1].astype(int)
        v = E[yi, xi]
        self.lit = v >= np.quantile(v, quantile)
        self.chance = float(self.lit.mean())       # でたらめに置いたときの乗る率

        # 画素 → 最も近いセルが点いているか（一度だけ作る）
        tree = cKDTree(self.P)
        gy, gx = np.mgrid[0:H, 0:W]
        _, idx = tree.query(np.stack([gx.ravel(), gy.ravel()], 1))
        self.litmap = self.lit[idx].reshape(H, W)

        ys, xs = np.nonzero(self.lit[:, None].ravel()[np.arange(len(self.P))][:, None])[0], None
        self.lit_xy = self.P[self.lit]

    def land(self, pts: np.ndarray) -> np.ndarray:
        """型紙の点 → 乗ったかどうか（画面外は乗らない扱い）。"""
        H, W = self.shape
        x = np.clip(np.round(pts[:, 0]).astype(int), 0, W - 1)
        y = np.clip(np.round(pts[:, 1]).astype(int), 0, H - 1)
        ok = (pts[:, 0] >= 0) & (pts[:, 0] < W) & (pts[:, 1] >= 0) & (pts[:, 1] < H)
        return self.litmap[y, x] & ok


# ══ 型紙（設計から作る。図を描いたコードは見ない）═══════════
def unit_ring(step_units: float, ratio: float = RATIOS["phi2"]):
    """外接半径 1・向き 0° の円環。点列と、その点が何番地のものかを返す。"""
    pts, own = [], []
    for k in range(10):
        c = ratio * np.array([math.cos(36 * k * DEG), math.sin(36 * k * DEG)])
        beta = 36 * k                                  # 偶奇で 0/36 になる
        V = [c + np.array([math.cos((72 * j + beta) * DEG),
                           math.sin((72 * j + beta) * DEG)]) for j in range(5)]
        for j in range(5):
            a, b = V[j], V[(j + 1) % 5]
            n = max(2, int(np.linalg.norm(b - a) / step_units))
            for t in np.linspace(0, 1, n + 1):
                pts.append(a + t * (b - a))
                own.append(k)
    return np.array(pts), np.array(own)


def place(unit: np.ndarray, R: float, theta_deg: float, cx: float, cy: float):
    t = theta_deg * DEG
    c, s = math.cos(t), math.sin(t)
    return R * (unit @ np.array([[c, s], [-s, c]])) + np.array([cx, cy])


# ══ 合わせと検図 ════════════════════════════════════════════
def fit(ret: Retina, ratio: float, step_px: float,
        R_range, theta_range, centers):
    best = None
    cache = {}
    for R in R_range:
        key = round(R, 3)
        if key not in cache:
            cache[key] = unit_ring(step_px / R, ratio)
        u, own = cache[key]
        for th in theta_range:
            for (cx, cy) in centers:
                v = ret.land(place(u, R, th, cx, cy)).mean()
                if best is None or v > best[0]:
                    best = (float(v), float(R), float(th), float(cx), float(cy))
    return best


def inspect(path: str, cell: float, quantile: float, ratio_key: str,
            verbose: bool = True, box=None):
    L = to_gray(path)
    if box:
        x0, y0, x1, y1 = box
        L = L[y0:y1, x0:x1]
    H, W = L.shape
    ret = Retina(L, cell=cell, quantile=quantile)
    step = 1.618 * cell                                 # 五角形の間隔より細かくしない
    ratio = RATIOS[ratio_key]

    # 図から取る三つ：中心・尺度・向き（設計値は使わず、三段で掃く）
    cx0, cy0 = ret.lit_xy[:, 0].mean(), ret.lit_xy[:, 1].mean()
    m = min(W, H)
    g = np.linspace(-0.10 * m, 0.10 * m, 7)
    v, R, th, cx, cy = fit(ret, ratio, step,
                           np.geomspace(m * 0.02, m * 0.20, 24),
                           np.arange(0, 36, 3.0),
                           [(cx0 + dx, cy0 + dy) for dx in g for dy in g])
    d = float(g[1] - g[0])
    v, R, th, cx, cy = fit(ret, ratio, step,
                           np.arange(max(2 * cell, R * 0.85), R * 1.16, cell),
                           np.arange(th - 3, th + 3.01, 1.0),
                           [(cx + dx, cy + dy)
                            for dx in np.arange(-d, d + 1, 2 * cell)
                            for dy in np.arange(-d, d + 1, 2 * cell)])
    v, R, th, cx, cy = fit(ret, ratio, step,
                           np.arange(max(2 * cell, R - 3 * cell), R + 3 * cell, cell / 3),
                           np.arange(th - 1.5, th + 1.51, 0.25),
                           [(cx + dx, cy + dy)
                            for dx in np.arange(-2 * cell, 2.01 * cell, cell / 2)
                            for dy in np.arange(-2 * cell, 2.01 * cell, cell / 2)])

    u, own = unit_ring(step / R, ratio)
    on = ret.land(place(u, R, th, cx, cy))
    per = np.array([on[own == k].mean() for k in range(10)])

    if verbose:
        print(f"画像 {path}  {W}×{H}")
        print(f"  担体 {len(ret.P)} セル（最近接 {cell}px・画素の "
              f"{100*len(ret.P)/(W*H):.1f}%）／点いたセル {int(ret.lit.sum())}")
        print(f"  でたらめに置いたときの乗る率（偶然） {ret.chance:.3f}")
        print(f"\n  [合わせ] 中心 ({cx:.1f},{cy:.1f})／外接半径 {R:.1f}px／向き {th:.2f}°")
        print(f"  [検図]   乗る率 {on.mean():.3f}"
              f"（偶然の {on.mean()/max(ret.chance,1e-9):.2f}倍）")
        print("           番地ごと：" + " ".join(f"{p:.2f}" for p in per))
        print(f"           最も落ちた番地 {int(np.argmin(per))}"
              f"（{per.min():.2f}）／最も乗った番地 {int(np.argmax(per))}（{per.max():.2f}）")

        print("\n  [必ず落ちる設定] ── ここが落ちなければ、この検査は判定になっていない")
        for nm, args in (("向きを 9° ずらす", (R, th + 9, cx, cy)),
                         ("向きを 18° ずらす", (R, th + 18, cx, cy)),
                         ("尺度を 1.05 倍", (R * 1.05, th, cx, cy)),
                         ("尺度を φ 倍", (R * PHI, th, cx, cy)),
                         ("尺度を 1/φ 倍", (R / PHI, th, cx, cy)),
                         ("中心を 3セルずらす", (R, th, cx + 3 * cell, cy)),
                         ("中心を 半セルずらす", (R, th, cx + cell / 2, cy))):
            uu, _ = unit_ring(step / args[0], ratio)
            print(f"    {nm:20s} 乗る率 {ret.land(place(uu, *args)).mean():.3f}")
    return dict(score=float(on.mean()), chance=ret.chance, R=R, theta=th,
                cx=cx, cy=cy, per_address=per.tolist())


# ══ 設計から描く（自己検査用）═══════════════════════════════
def draw_ring(path: str, size: int = 600, R: float = 60.0,
              ratio: float = RATIOS["phi2"], width: int = 2):
    from PIL import Image, ImageDraw
    im = Image.new("RGB", (size, size), (255, 255, 255))
    d = ImageDraw.Draw(im)
    c0 = size / 2.0
    cen = []
    for k in range(10):
        c = np.array([c0, c0]) + ratio * R * np.array(
            [math.cos(36 * k * DEG), math.sin(36 * k * DEG)])
        beta = 36 * k
        V = [c + R * np.array([math.cos((72 * j + beta) * DEG),
                               math.sin((72 * j + beta) * DEG)]) for j in range(5)]
        d.line([tuple(v) for v in V] + [tuple(V[0])], fill=(20, 20, 20), width=width)
        cen.append(c)
    im.save(path)
    cen = np.array(cen)
    dd = [np.hypot(*(cen[(i + 1) % 10] - cen[i])) / R for i in range(10)]
    return dict(gap_ratio_min=min(dd), gap_ratio_max=max(dd),
                inradius2=2 * R * math.cos(36 * DEG), center_ring=ratio)


def selftest(cell: float, quantile: float):
    import tempfile, os
    p = os.path.join(tempfile.gettempdir(), "penrose_ring_selftest.png")
    info = draw_ring(p)
    print("[設計側の数]（画像は読んでいない）")
    print(f"  隣の中心距離 ÷ 外接半径：{info['gap_ratio_min']:.6f} 〜 "
          f"{info['gap_ratio_max']:.6f}（設計 φ = {PHI:.6f}）")
    print(f"  中心が乗る円 ÷ 外接半径：{info['center_ring']:.6f}"
          f"（設計 φ² = {PHI**2:.6f}）")
    print(f"  中心距離 {info['gap_ratio_min']*60:.2f}px ／ 内接半径の2倍 "
          f"{info['inradius2']:.2f}px → "
          f"{'ちょうど接する' if abs(info['gap_ratio_min']*60-info['inradius2'])<1e-6 else 'ずれている'}")
    print(f"\n描いた PNG：{p}")
    print("これ以降、描画コードは使わない。画像だけを装置に渡す。\n")
    r = inspect(p, cell, quantile, "phi2")
    print(f"\n判定：乗る率 {r['score']:.3f}／偶然 {r['chance']:.3f}"
          f"／比 {r['score']/max(r['chance'],1e-9):.2f}倍")
    return r


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--image", help="検める図（PNG/JPG）。省略すると自己検査")
    ap.add_argument("--draw", metavar="PATH", help="設計どおりの円環を描いて保存する")
    ap.add_argument("--ratio", choices=list(RATIOS), default="phi2",
                    help="中心円の比。phi2＝辺で接する（設計）／cot18＝頂点で触れる")
    ap.add_argument("--scan-ratio", action="store_true",
                    help="二つの比を両方当てて、どちらで張られているかを見る")
    ap.add_argument("--box", help="先に切り出す（例 --box 142,129,1336,1267）。"
                    "matplotlib の枠・格子・目盛りは切っておく")
    ap.add_argument("--cell", type=float, default=4.0, help="担体の最近接（画素）")
    ap.add_argument("--quantile", type=float, default=0.85, help="点灯を決める上位")
    a = ap.parse_args()

    if a.draw:
        info = draw_ring(a.draw)
        print(f"{a.draw} に保存した（隣の中心距離 ÷ 外接半径 = "
              f"{info['gap_ratio_min']:.6f}）")
        if not a.image:
            return
    if not a.image:
        selftest(a.cell, a.quantile)
        return
    box = tuple(int(v) for v in a.box.split(",")) if a.box else None
    if a.scan_ratio:
        res = {}
        for key in RATIOS:
            print(f"\n═══ 中心円の比 = {key}（{RATIOS[key]:.6f}）═══")
            res[key] = inspect(a.image, a.cell, a.quantile, key, box=box)
        best = max(res, key=lambda k: res[k]["score"])
        print(f"\n乗る率が高いのは {best}（"
              + "／".join(f"{k} {res[k]['score']:.3f}" for k in res) + "）")
    else:
        inspect(a.image, a.cell, a.quantile, a.ratio, box=box)


if __name__ == "__main__":
    main()
