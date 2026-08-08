# 走らせ方

```bash
python3 b13_two_tilings.py      # 60環を作る（rings_integer.json を書く）
python3 expand3.py              # 1245環へ展開（expanded.json を書く）
```

この2本を先に走らせると、担体を読む側（closed_walk / direction / still / pointing）が動きます。
担体を要らない側（wave_fit / approach / track_span / band / zoomstage）は単独で走ります。

| スクリプト | 記事 | 担体 |
|---|---|---|
| `qphi.py` | Q(φ)/Z[ζ₁₀] の厳密演算 | ─（他が import） |
| `b13_chain_units.py` | 910環の積み上げ | ─（similar_probe が import） |
| `b13_two_tilings.py` → `expand3.py` | 担体の生成 | 作る側 |
| `closed_walk.py` | §2 開いた歩きと閉じる歩き | 60環 |
| `wave_fit.py` | §3 波の骨格・24通りの族・比の不変性 | 不要 |
| `approach.py` | §2 寄せの階段（単数） | 不要 |
| `direction.py` | §2 方向の読み | 60環 |
| `track_span.py` / `band.py` | §3 追える周期数・帯 | 不要 |
| `zoomstage.py` | §3 一度に跳べる段 | 不要 |
| `still.py` | §4 静止した輪郭 | 60環 |
| `still_symbols.py` | §4 切替回数の独立な数え直し（図の検算） | 引数で受ける |
| `similar_probe.py` | §4 型紙 | 910環（内部で作る） |
| `pointing.py` | §5 局所の見えの分類 | **1245環**（先に expand3.py） |
| `draw_retina_figs.py` | 図1・図2・図3（両眼視） | 内部で build() |

各スクリプトの冒頭に事前登録（何が出たら何が言えるか）があります。
