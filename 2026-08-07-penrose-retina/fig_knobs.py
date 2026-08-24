#  図：二本のつまみ ── 膜の緩和τと視線の減衰γ
#
#  検定G6（事前登録）
#    OK なら：左パネルの2本の折れ線が交差する（γの最適が対象の動きで入れ替わる）。
#             右パネルの版面が重ならない。和文が豆腐にならない。
#    NG なら：交差しない／重なる／文字が出ない → 図を出さない
#
#  数値は 検定FL6（膜つき τ=2・合成場面・背景100・ちらつき8・60枚・一歩0.45・各3種）

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Circle
import numpy as np

plt.rcParams["font.family"] = ["Noto Sans CJK JP"]
plt.rcParams["axes.unicode_minus"] = False

X = [0, 1, 2, 3, 4]
XL = ["向きが\n変わらない\n（電車）", "8枚ごと", "4枚ごと", "2枚ごと",
      "毎枚変わる\n（猫じゃらし）"]
G15 = [60, 60, 43, 43, 26]      # γ=1.5（慣性を残す）
G12 = [16, 16, 29, 38, 49]      # γ=12 （慣性を捨てる）

assert (G15[0] > G12[0]) and (G15[-1] < G12[-1]), "交差しない → 図を出さない"

fig = plt.figure(figsize=(12.6, 5.0), dpi=170)
gs = fig.add_gridspec(1, 2, width_ratios=[1.08, 1.0], wspace=0.24,
                      left=0.065, right=0.985, top=0.86, bottom=0.17)

# ---- 左：交差 ------------------------------------------------------------
ax = fig.add_subplot(gs[0, 0])
ax.plot(X, G15, "o-", lw=2.6, ms=8, color="#c0392b", label="γ 小（慣性を残す）")
ax.plot(X, G12, "s-", lw=2.6, ms=8, color="#1f6fb4", label="γ 大（慣性を捨てる）")
ax.fill_between(X, G15, G12, where=np.array(G15) > np.array(G12),
                color="#c0392b", alpha=0.10, interpolate=True)
ax.fill_between(X, G15, G12, where=np.array(G15) <= np.array(G12),
                color="#1f6fb4", alpha=0.10, interpolate=True)
ax.set_xticks(X); ax.set_xticklabels(XL, fontsize=10.5)
ax.set_ylim(0, 68); ax.set_ylabel("見失うまでの枚数（60枚中）", fontsize=11.5)
ax.set_xlabel("対象の動きの読みにくさ →", fontsize=11.5, labelpad=8)
ax.set_title("減衰 γ の最適は、対象の動きで入れ替わる", fontsize=13, pad=12)
ax.legend(loc="lower left", fontsize=10.5, framealpha=0.95)
ax.grid(axis="y", alpha=0.25)
ax.text(0.25, 63, "流し撮り", fontsize=11, color="#c0392b", weight="bold")
ax.text(3.35, 55, "食らいつく", fontsize=11, color="#1f6fb4", weight="bold")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)

# ---- 右：何を固定し、何を動かすか ---------------------------------------
ax2 = fig.add_subplot(gs[0, 1]); ax2.axis("off")
ax2.set_xlim(0, 10); ax2.set_ylim(0, 10)
ax2.set_title("膜の緩和 τ と 視線の減衰 γ", fontsize=13, pad=12)

# 膜（τ）
xs = np.linspace(0.4, 9.6, 400)
def sheet(cx, d, w=1.1):
    return 7.4 - d * np.exp(-((xs - cx) ** 2) / (2 * w * w))

ax2.plot(xs, sheet(6.6, 1.55), color="#333", lw=2.4)
ax2.plot(xs, np.full_like(xs, 7.4), color="#999", lw=1.4, ls=(0, (5, 4)))
ax2.text(2.0, 7.75, "動かない荷 → τ 枚で均される", fontsize=10.5, color="#666")
ax2.annotate("", xy=(2.4, 7.4), xytext=(2.4, 6.55),
             arrowprops=dict(arrowstyle="-|>", color="#666", lw=1.6))
ax2.text(5.35, 5.35, "動く荷", fontsize=10.5, color="#333")
ax2.add_patch(Circle((6.6, 6.05), 0.30, color="#c0392b", zorder=5))
ax2.text(0.35, 8.9, "τ ── 膜が荷に追いつく枚数", fontsize=12, weight="bold")
ax2.text(0.35, 8.35, "動かない背景は凹みにならない", fontsize=10.5, color="#555")

# 視線（γ）
ax2.text(0.35, 4.35, "γ ── 落ちる球に残す慣性", fontsize=12, weight="bold")
for cx, d, col, lab, vlen in ((2.6, 1.35, "#c0392b", "γ 小：井戸に引かれたまま併走", 1.5),
                              (7.4, 1.35, "#1f6fb4", "γ 大：底に留まる", 0.0)):
    ys = 2.5 - d * np.exp(-((xs - cx) ** 2) / (2 * 0.95 ** 2))
    ax2.plot(xs, ys, color="#333", lw=2.2)
    ax2.add_patch(Circle((cx, 2.5 - d + 0.26), 0.26, color=col, zorder=5))
    if vlen:
        ax2.add_patch(FancyArrowPatch((cx + 0.35, 2.5 - d + 0.26),
                                      (cx + 0.35 + vlen, 2.5 - d + 0.26),
                                      arrowstyle="-|>", mutation_scale=15,
                                      color=col, lw=2.0))
    ax2.text(cx, 0.42, lab, fontsize=10, color=col, ha="center")

# 版面の重なりを見る（検定G6）
fig.canvas.draw()
b1 = ax.get_window_extent(); b2 = ax2.get_window_extent()
assert b1.x1 <= b2.x0 + 1, "版面が重なっている → 図を出さない"

fig.savefig("/mnt/user-data/outputs/fig_knobs.png", facecolor="white")
print("検定G6 OK  交差あり／版面の重なりなし")
print("  γ小 端点 %d → %d ／ γ大 端点 %d → %d" % (G15[0], G15[-1], G12[0], G12[-1]))
