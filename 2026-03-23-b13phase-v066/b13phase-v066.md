## b13phase v0.6.6 — 実行性能改善

## 変更概要

`evaluator_proto.py` と `phase_packed_u64.py` の2ファイルを修正。  
他のモジュールは v0.6.5 から無変更。

---

## 背景

v0.6.5 で「B13がFFTより遅い」という観測があった。  
実装を調べると、遅さの原因は理論側ではなく実行パス側にあった。

`fractal_vertex_recurrence.py` にはすでに次の式が証明済みとして書かれている。

```
az_idx(level n) = (az_0 + n × 312) % 3120
```

これは「az_idx の整数加算1回 + COS/SIN テーブル引き1回」で位相が評価できることを意味する。  
ところが `fractal_cos_sin_proto()` は下位桁の補正ループを回しており、この近道を使っていなかった。  
同様に `add_u64_packed()` も、packed u64 をいったん unpack してからループし、また pack するという迂回をしていた。

なお、b13phase はクライオEM実データの B13-FFT 残渣解析パイプライン（`run_b13_fft_v63_final.py`）に組み込まれて実際に動いている。今回の改善はそのパイプラインの上流にある評価器の速度に直接影響する。

---

## 変更内容

`evaluator_proto.py` は `evaluator.py`（production版）と並存している。  
今回は `evaluator_proto.py` 側の実行パスを、`evaluator.py` が既に採用している設計（az_idx の直接計算 + テーブル引き）と整合するように修正した。

### `evaluator_proto.py`

**変更前** — 下位桁ごとに `BASE^l` の冪乗を計算して補正を適用するループ。

```python
coarse = digits[0]
cx = COS_TABLE[coarse]
cy = SIN_TABLE[coarse]

for l in range(1, len(digits)):
    fine = digits[l]
    scale = BASE ** l
    new_cx = (cx * scale - cy * fine) // scale
    new_cy = (cy * scale + cx * fine) // scale
    cx, cy = new_cx, new_cy
```

**変更後** — MSD を az_idx として直接テーブルを引く。

```python
az = digits[0]
return COS_TABLE[az], SIN_TABLE[az]
```

あわせて `fractal_cos_sin_az(az_idx: int)` を追加。  
漸化式から求めた az_idx を整数のまま渡せる。

```python
# 使用例
az = (az_0 + level * 312) % 3120
cx, cy = fractal_cos_sin_az(az)
```

---

### `phase_packed_u64.py`

**変更前** — `add_u64_packed` / `sub_u64_packed` が unpack → 桁ループ → pack の3段階。

**変更後** — 5桁をビット演算で同時展開し、carry / borrow をその場で伝播。

```python
def add_u64_packed(a: int, b: int) -> tuple[int, int]:
    d0 = ((a >> 48) & MASK) + ((b >> 48) & MASK)
    d1 = ((a >> 36) & MASK) + ((b >> 36) & MASK)
    d2 = ((a >> 24) & MASK) + ((b >> 24) & MASK)
    d3 = ((a >> 12) & MASK) + ((b >> 12) & MASK)
    d4 = ( a        & MASK) + ( b        & MASK)
    c = 0
    if d4 >= BASE: d4 -= BASE; c = 1
    d3 += c; c = 0
    if d3 >= BASE: d3 -= BASE; c = 1
    d2 += c; c = 0
    if d2 >= BASE: d2 -= BASE; c = 1
    d1 += c; c = 0
    if d1 >= BASE: d1 -= BASE; c = 1
    d0 += c; carry = 0
    if d0 >= BASE: d0 -= BASE; carry = 1
    return (d0 << 48) | (d1 << 36) | (d2 << 24) | (d3 << 12) | d4, carry
```

---

## 計測結果

v0.6.5 の現物（`b13phase_v065.zip`）に対して直接パッチして計測。  
Python、`timeit`、N=200,000回。

| 関数 | v0.6.5 | v0.6.6 | 倍率 |
|---|---|---|---|
| `add_u64_packed` | 3.48 µs | 1.35 µs | 2.6x |
| `sub_u64_packed` | 3.51 µs | 1.34 µs | 2.6x |
| `fractal_cos_sin_proto`（3桁） | 1.07 µs | 0.22 µs | 4.9x |
| `fractal_cos_sin_az`（新API） | — | 0.16 µs | — |

正確性確認: add 20,000ケース、sub 20,000ケース、proto 全3,120点でエラーなし。

---

## FFT との速度比較について

B13 は FFT の置き換えを目指していない。  
FFT はすでに社会に実装されており、エンジニアは使い慣れている。

B13 の用途として想定しているのは**残渣解析**——FFT を適用した後に残る構造の解釈——であり、そこで B13 が実用的な速度で動くかどうかが問題だった。

この変更後、60頂点1層カーネルの速度を参考値として `numpy.fft.fft(N=60)` と並べると次のようになった。

| | 時間 |
|---|---|
| B13 60頂点 1層（az 更新 + テーブル引き） | 5.3 µs |
| `numpy.fft.fft` (N=60) | 7.7 µs |

v0.6.5 では同じカーネルが遅い実行パスのせいでこの数倍かかっていた。  
v0.6.6 でその差は解消した。

「B13がFFTより速い」とは言わない。  
ただ、**FFT を使い慣れたエンジニアが残渣処理に B13 を試してみようと思えるかどうか**、その速度域には入った。

---

## 今後の課題

- `ACTIVE_MASK` / `BAND_ID` / `ANTIPODE` を B13理論に沿って充填する
- `connection.py` をビルド専用生成器（`compile_tables.py`）と実行時 ROM（`rom_tables.py`）に分離し、実行時の構造再発見コストをゼロにする
- 上記が揃った段階で、FFT との全体変換コスト比較を再度行う
