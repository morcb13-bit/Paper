# B13 engine 契約 v160（実機 `engine_extract` がハネースに結合するための条件）

凍結された harness（v160）に実機 BpVM/B13 を差し込むための **入口契約**。
これを満たさない engine 出力は harness 入口で REJECT される。捏造はしない——
engine が契約を満たした時点で初めて real_d を測り、#5 baseline に重ねる。

---

## 0. プラグ点（唯一）

```python
def engine_extract(signal: np.ndarray) -> tuple[np.ndarray, int]:
    """raw 時系列波形 → (residue_walk, detected_s)"""
    return residue_walk, detected_s
```

harness 側は `evaluate(signal, engine_extract)` を呼ぶだけ。engine 内部構造は不問。
入口で見るのは **入力契約・出力契約・粒度契約** の三つ。

---

## 1. 入力契約（harness が engine に渡すもの・条件1,2）

| 項目 | 仕様 | 根拠 |
|---|---|---|
| 型 | `np.ndarray`, 1次元, `float` | 条件1（residue 列・整数特徴量は REJECT） |
| 長さ | `>= 4*FRAME = 1024` サンプル | 条件1（短い＝特徴量の疑い） |
| 中身 | raw 時系列波形（位相構造を保持） | 条件2（null は raw への phase-shuffle で作る） |

engine は **時系列を必ず経由すること。** residue 化済みの列を受けて返す実装は、
null（phase-shuffle）が壊す対象（時間領域の位相構造）を消すため接続不可（1-B 定理）。

---

## 2. 出力契約（engine が返すもの・条件3,4）

| 要素 | 仕様 | 違反時 |
|---|---|---|
| `residue_walk` | 整数 `ndarray`, 値域 `[0,13)`, 非空 | 条件3 REJECT |
| `detected_s` | `int`, `1..12`（s=0 は readout 対象外） | 条件4 REJECT |
| `detected_s is None` | **不可**（全 s 探索は別 lane・天井再較正が要る） | 条件4 REJECT |

`detected_s` は **engine が信号から検出した値** であること。harness 側で全 s を
走査すると多重比較で天井が甘くなる。s は engine が一本に決めて持ち込む。

---

## 3. 粒度契約 〔v160 新規・最重要 defeater〕

**現 readout（`directional_localization_score`）は toy の run 粒度で較正されている。**
`run_factor = 2/max(2,n)`（n=安定ラン数）と `best_cov`（隣接 (s,−s) ペア一組のみ credit）は、
**真の反転が ≈2 本の支配ランに畳まれた walk** を前提にする。

### 検証済み（合成 walk・捏造でない）

| walk 長 | jitter | n_runs | best_cov | real_d |
|---|---|---|---|---|
| 24 | 0.00 | 2.0 | 1.000 | **1.000** |
| 24 | 0.02 | 2.2 | 0.917 | 0.856 |
| 24 | 0.05 | 2.6 | 0.798 | 0.661 |
| 1200 | 0.00 | 2.0 | 1.000 | **1.000** |
| 1200 | 0.02 | 23.2 | 0.078 | **0.007** |
| 1200 | 0.05 | 47.5 | 0.040 | **0.002** |

**読み:** clean な反転は長さに依らず real_d=1.0。だが長さ1200で residue が 2% ジッタすると
ラン数が 23 に膨れ（run_factor=2/23）かつ best_cov が単一ペアに希釈され、**反転が本物でも
real_d≈0.007 → BLIND_DROP（偽陰性）**。長さ24では同 jitter で 0.856 と優雅に劣化する。

### 帰結（engine 側が選ぶ・二択）

- **(A) toy 粒度で渡す〔推奨〕：** engine は walk を **run-level に畳んで** 渡す
  （連長圧縮 / フレーム mode 化 / cleanup 後）。真の反転が ≈2 本の支配ランになる粒度。
  → 凍結された #5 baseline がそのまま適用できる。**追加作業なし。**
- **(B) 生 length-1200 を渡す：** その場合 **#5 baseline は無効**。
  readout を再導出（全 s→−s 遷移の集約 coverage ＋ 長さ頑健な run_factor）し、
  #4/#5 を当該粒度で再較正してから結合すること。**baseline 流用は禁止。**

**この契約を破ると、本物の反転が静かに偽陰性で落ちる。** engine 投入前に粒度を確定せよ。

---

## 4. 判定意味論（2×2・engine 投入後に harness が出すもの）

`gap99 = real_d − surr_99`（surr=raw phase-shuffle・99pct 天井）。`EPS_SEP=0.05`, `REAL_OK=0.90`。

| | gap99 > ε（分離あり） | gap99 ≤ ε（分離なし） |
|---|---|---|
| **real_d ≥ 0.90** | `ADVANCE` 前進 | `BLIND_SUSPECT_REJECT` 信号盲疑い |
| **real_d < 0.90** | `SIGNAL_WEAK` 薄い陽性 | `BLIND_DROP` 盲目 |

- **ADVANCE** が出て初めて reversal_lane §4 候補が前進。動作 SNR を #5 baseline に重ねて読む。
- **BLIND_SUSPECT_REJECT**：real_d が高くても null が分離しない＝信号依存性が無い。信用しない。
- engine 投入の合格条件は **動作 SNR で ADVANCE を覆うこと**。それ未満は前進と declare しない。

---

## 5. 結合手順（次セッション / リポジトリ連携）

1. engine 側で粒度契約（§3）を確定——(A) run-level 化 or (B) 再較正コミット。
2. `engine_extract` を §1,§2 準拠で実装。
3. `evaluate(real_signal, engine_extract)` を動作 SNR 数点で回す。
4. 出る `lane_verdict` を読む。ADVANCE なら #5 baseline に重ね、本ファイルと引継書 v160 に畳む。
5. **それまで real_d は作らない。入口待ちが正しい停止位置。**
