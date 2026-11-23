# BepiColombo: 8 Years, 9 Flybys, Into Mercury's Orbit

Animated figures and a small integer-only orbit check for the ESA/JAXA BepiColombo mission, which entered orbit around Mercury on 21 November 2026.

- **Figure 1** – The eight-year path around the Sun, with the nine flybys (Earth ×1, Venus ×2, Mercury ×6).
- **Figure 2** – Each flyby as a hyperbola seen from the planet, and how the approach speed v∞ shrinks over the Mercury flybys.
- **Figure 3** – Two years running alongside Mercury (data), then capture and separation of MPO and Mio (projection from planned values).
- **Orbit check** – 20 days of the spacecraft's orbit propagated with integers and addition only (B13 library), compared against JPL Horizons.

Articles: Zenn (Japanese) – *link* / DEV (English) – *link*

## Files

| File | What it is |
|---|---|
| `bepicolombo_path.html` / `_en.html` | Figure 1 (Japanese / English). Single HTML file, no external libraries; paste into CodePen as is. |
| `bepicolombo_flybys.html` / `_en.html` | Figure 2 |
| `bepicolombo_capture.html` / `_en.html` | Figure 3 |
| `bepicolombo_path_en.gif`, `bepicolombo_flybys_en.gif`, `bepicolombo_capture_en.gif` | GIF versions of the English figures |
| `bepicolombo_cover_en.gif` | Cover image (1000 × 420) |
| `b13_orbit_sim.py` | Integer-only orbit check (see below) |
| `helio.csv` | JPL Horizons samples: Julian date, x, y (AU), vx, vy (AU/day), every 20 days, Sun-centred, ecliptic J2000 |
| `bepi_flybys.py` | Table of orbital elements before/after each flyby, speeds, and 13^k addresses (`table.md` is its output) |
| `build_data.py` | Interpolates `helio.csv` to daily steps and adds planet positions (writes `frames.npz`) |
| `build_figs.py` | Builds the integer data for Figures 1 and 2 (`fig1_data.json`, `fig2_data.json`) |
| `build_fig3.py` | Builds the integer data for Figure 3 (`fig3_data.json`) |
| `build_html.py` | Assembles the Japanese Figures 1 and 2 from `fig1_tpl.html`, `fig2_tpl.html` and `css.txt` |
| `build_fig3_html.py` | Assembles Figure 3 in Japanese and English from `fig3_tpl.html` |

The English HTML files for Figures 1 and 2 are the Japanese versions with the text replaced by English.

## How the figures work

All orbital computation is done beforehand in Python. The results are embedded in each HTML file as integers: 1/10000 AU in Figure 1, 1/100 planet radius in Figure 2, and km in Figure 3. Inside the page, the animation only adds 1 to the frame number and looks up the table. Decimals are used only to place points on the screen. Each page has a speed slider; one HTML file holds one figure and stays under 1 MB.

To rebuild:

```bash
pip install numpy astropy
python3 build_data.py      # frames.npz
python3 build_figs.py      # fig1_data.json, fig2_data.json
python3 build_fig3.py      # fig3_data.json
python3 build_html.py      # Japanese Figures 1 and 2
python3 build_fig3_html.py # Figure 3 (Japanese and English)
```

## The integer-only orbit check

`b13_orbit_sim.py` propagates BepiColombo from 8 to 28 December 2024, between Mercury flybys 5 and 6, when the orbital elements barely change (no engine burns). It uses only the Sun's gravity and takes 320 leapfrog steps of 1/16 day. Positions are integers with 1 AU = 10¹². All arithmetic runs on balanced base-13 digit lists from the B13 library (`b13_balanced.py`, `b13_arith.py`), where multiplication, division and square roots are built from addition.

Pass criteria were fixed before running:

- Pass if the result is within 0.0005 AU of the Horizons position on 28 December.
- Negative control: the same run with gravity set to zero must miss by more than 0.05 AU.
- The integer run must agree with the same steps in floating point to within 1e-9 AU.

```bash
B13LIB=/path/to/b13-library/scripts python3 b13_orbit_sim.py --en
```

```
Position on 12-28, Horizons : [-0.395217, -0.057955, 0.020093]
Position on 12-28, B13 int  : [-0.395223, -0.057954, 0.020050]   run time 4.0 s
Difference (B13 - Horizons) : 0.000043 AU = 6501 km
Negative control, GM = 0    : 0.500059 AU
Integer vs floating-point   : 1.24e-10 AU
Pass (< 0.0005): True / control misses (> 0.05): True
```

Without `--en` the output is in Japanese. If `B13LIB` is not set, the script looks for the library in its own folder.

## Data and limits

- Spacecraft positions and orbital elements come from [JPL Horizons](https://ssd.jpl.nasa.gov/horizons/) (target −121, Sun-centred and planet-centred, ecliptic J2000). Planet positions come from the astropy built-in ephemeris.
- Positions between the 20-day samples are interpolated, so the spacecraft–planet distance in the figures carries an error of about 0.005 AU (a few pixels on screen).
- The second half of Figure 3 (capture, separation, final orbits) is a projection. Separation (9–10 December 2026) and the MPO's final orbit (March 2027) had not happened at the time of writing. Orbit sizes and periods are the published planned values; orbit orientation and the capture path are not drawn.
- The orbit check includes only the Sun's gravity. It is a simple check of the orbit's skeleton, not a replacement for mission navigation.

---

## 日本語

2026 年 11 月 21 日に水星の周回軌道に入った日欧共同の水星探査計画 BepiColombo を、動く図 3 枚と、整数と足し算だけで軌道を進める簡単な照合でたどるリポジトリです。

- **図 1**：太陽から見た 8 年の道筋と 9 回のスイングバイ
- **図 2**：惑星から見た 9 本の双曲線と、水星に近づく速さ v∞ の下がり方
- **図 3**：水星と並んで回った 2 年（データ）と、捕獲・分離（計画値による予想図）
- **軌道の照合**：B13 ライブラリ（平衡 13 進、足し算だけ）で 20 日ぶんの軌道を進め、JPL Horizons と突き合わせる

ファイル名に `_en` が付いているのが英語版、付いていないのが日本語版です。HTML は 1 本に図 1 枚で、外部のライブラリを読み込まないので、そのまま CodePen に貼れます。図の中の計算は、コマの番号に 1 を足して表を引くだけで、画面の座標への換算だけ小数を使っています。

照合の結果は、20 日後の位置が Horizons と 0.000043 AU（約 6,500 km）の差で、事前に決めた合格基準（0.0005 AU 以内）に入りました。重力をゼロにした負の対照は 0.50 AU 外れ、整数版と浮動小数版の差は 1.2×10⁻¹⁰ AU でした。計算は普通のパソコンで約 4 秒です。

図 3 の後半は、分離（2026 年 12 月 9〜10 日）と MPO の最終軌道（2027 年 3 月）がまだ先なので、計画値で描いた予想図です。
