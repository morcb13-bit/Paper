# 2026-06-20-BpVM / README

reversal_lane harness（v167）と、岩手県沖 M6.9（2026-06-25 07:30）K=3 bundle 実行結果の保管。
**解釈なし。数字と手順だけ置く。** 読み方の規律は §6 に集めた。

---

## 1. 何が置いてあるか

```
2026-06-20-BpVM/
├── README.md                       ← 本ファイル
├── code/                           ← v167 reversal_lane bundle（11ファイル）
│   ├── engine_reversal_eval.py     ← canonical engine
│   ├── engine_b13.py               ← 実エンジン本体（detect_s 含む）
│   ├── engine_admissibility.py     ← 6条件 admissibility gate
│   ├── bundle_b13.py               ← B-3 束定義 v0（M=3 / unanimous）
│   ├── a_prep.py                   ← 実観測入口 harness（K_EXPECT=3）
│   ├── guard_b13.py                ← B-2 棄権ガード（三値 verdict）
│   ├── log_first_real.py
│   ├── test_reversal_directional.py
│   ├── test_reversal_isolation.py  ← P=13, FRAME=256, HOP=256 固定
│   ├── test_reversal_localization.py
│   └── B13_引継書_v167.md
├── data/20260625073000/
│   ├── knet/csv/MYG0012606250730.csv
│   ├── knet/csv/IWT0032606250730.csv
│   └── kik/csv/IWTH032606250730.csv
├── runner/run_realD.py             ← 本 README §4 の最小実行スクリプト
├── result/20260625073000_K3_real_d.log  ← 生実行ログ
└── docs/
    ├── B13_引継書_v168.md           ← gate 完成までの経緯
    ├── B13_v168_bundle_section.md   ← K=3 bundle 名簿（観測点・成分・距離）
    └── B13_引継書_v169.md           ← 本実行の結果と読み
```

---

## 2. canonical 同定（必ず最初に確認）

engine の identity は**名前ではなく sha256** で。

```sh
$ sha256sum code/engine_reversal_eval.py
5d02397953339e7d2d0cf60fb6d4c99d85b8579e7af62c9c09f611ef3ef9feec  code/engine_reversal_eval.py
```

これと一致しなければ canonical ではない（stale engine）。**結果は無効**として扱う。

他ファイルの参考 sha256（v167 zip 展開直後）:

| file | sha256 prefix |
|---|---|
| engine_reversal_eval.py | `5d0239795333…` ← **canonical** |
| engine_admissibility.py | `3c46c1021b28…` |
| engine_b13.py           | `4bfd4055fbcf…` |
| bundle_b13.py           | `ba6f2001054c…` |
| a_prep.py               | `d40680094a98…` |
| guard_b13.py            | `094163912ee1…` |

---

## 3. データ仕様

NIED K-NET / KiK-net CSV。`#` 始まりは全てヘッダ。100Hz・raw 加速度（gal）。

**列の選択**（director 固定）:

| station | net | 採用列（0-indexed）| 意味 |
|---|---|---|---|
| MYG001 | K-NET   | `col=2` | 地表 N-S |
| IWT003 | K-NET   | `col=2` | 地表 N-S |
| IWTH03 | KiK-net | `col=5` | **地表 ch4**（ch1–3 は地中・除外）|

KiK-net 6成分は ch1–3=地中 / ch4–6=地表（IWTH03 では振幅で同定済）。地表と地中は **絶対に混ぜない**。

**震源（波形ヘッダ・3点一致）**: 2026-06-25 07:30:00, 40.2°N, 142.3°E, 50km, M6.9

---

## 4. 再現手順

前処理は **demean のみ**。filter / taper / window 切り出しは入れない。

```python
# runner/run_realD.py（要点だけ）
import sys; sys.path.insert(0, "code/")
import numpy as np
from a_prep import run_bundle

def load_col(path, col):
    vals = []
    with open(path) as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.rstrip("\n").split(",")
            if len(parts) > col: vals.append(float(parts[col]))
    return np.array(vals, dtype=np.float64)

root = "data/20260625073000"
sigs = [
    load_col(f"{root}/knet/csv/MYG0012606250730.csv", 2),
    load_col(f"{root}/kik/csv/IWTH032606250730.csv", 5),
    load_col(f"{root}/knet/csv/IWT0032606250730.csv", 2),
]
sigs = [x - x.mean() for x in sigs]  # demean
res = run_bundle(sigs, target="20260625073000", n_surr=200)  # seed=11 既定
```

`n_surr=200`, `seed=11` は `a_prep.run_bundle` の既定。**変えない**。

---

## 5. 結果（数字のみ・2026-06-27 実行）

| station    | verdict | s | real_d   | gap99    |
|---|---|---|---|---|
| MYG001_NS  | BLIND   | 4 | `0.0000` | `+0.0000` |
| IWTH03_ch4 | BLIND   | 2 | `0.0000` | `+0.0000` |
| IWT003_NS  | BLIND   | 6 | `0.0000` | `+0.0000` |

**BUNDLE = `BUNDLE_ABSTAIN`**（diag: 全本 BLIND）
reason: `real_d=0.000<0.9 (real-collapse)`（全本）

---

## 6. 読み方の規律（重要）

未来の自分・Chappy・他読者が誤読しないための明示。

1. **`BLIND` / `ABSTAIN` は honest output**。失敗ではない。harness が「無い」と言える状態を成功と呼ぶ。
2. **`s = 4, 2, 6` は支配 bin の readout**。`real_d=0` のため**反転署名ではない**。3点で割れている事実だけ記録に残す。
3. **`gate` と `real_d` を混ぜない**。gate は型・admissibility の検査、real_d は測定。
4. **K=2 で real_d を読まない**（`K_EXPECT=3`）。
5. **engine の identity は sha256 のみ**。ファイル名で判定しない。
6. **`FORWARD` が出ても発見扱いしない**——a_prep `PRIOR_DECLARATION`（不変 decoy → 偽 FORWARD の可能性）を先に疑う。
7. **B13 は instrument**。地震予測はしない。本結果は岩手県沖 M6.9 単一・K=3・特定3点に限られる。

---

## 7. このイベントについて言えないこと

- 「岩手県沖 M6.9 に B13 反転署名が**無い**」とは言えない（K=3 / SNR帯 / 単一イベント / 帯域非分離）。
- 「reversal_lane は地震に**使える / 使えない**」も本結果一発では決まらない。
- 他イベントへの外挿はしない（別 bundle で別途測る）。

---

## 8. 関連保留（別枠・本 bundle 不参加）

- **IWT019**（K-NET, 39.803°N 141.658°E, 震央70km, 地表）は IWTH03 とほぼ同位置。
  将来「同地点 K-NET 地表 vs KiK-net 地表/地中」比較実験用に保存。**本 K=3 には混ぜない**（bundle 独立性が濁る）。

---

## 9. 文書系譜

- v167（code 凍結・sha256 確定）
- v168（K=3 gate 完成・admissibility 全点 PASS・real_d 直前で停止）
- **v169（本実行・全 BLIND honest ABSTAIN）** ← この記録
