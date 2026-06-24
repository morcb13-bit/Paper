# B13 引継書 v159（engine-wiring 入口ゲート）

（v158 から継ぐ。v158 は #4 directional readout を engine の凍結 instrument として locked、
#2 を監査ログへ降格、reversal_lane を「正しい物差しで測った薄い機構候補」とした。
本セッションは **engine 配線の admissibility を director 固定**し、6 条件を executable gate に落とし、
**directional baseline（#5）を比較基準として確定**した。**engine 本体は未投入——捏造せず入口待ち。**
これは v159 の入口ゲート文書。engine コードが入った時点で real_d を測り、本ファイルに結果を畳む。）

---

## 0. このセッションの一行

物差しは固定済み（#4）。残るは engine 側一本：real 抽出が動作 SNR で 1.0 近傍を保つか。
**engine-wiring admissibility を 6 条件で固定**し、非準拠 engine を入口で弾くゲートを実装。
**directional baseline（SNR→gap99）を確定**——99pct 天井で全範囲正、-6.5dB でも +0.107、real-collapse 型。
**engine 本体は未投入。捏造しない。入口待ちが正しい停止位置。**

---

## 1. engine-wiring admissibility 〔director 固定・最重要〕

固定仕様:
  `raw time-series` → `engine_extract(signal) -> (residue_walk, detected_s)`
  → `directional_localization_score(residue_walk, detected_s)`（#4 凍結）
  → `phase_shuffle(raw signal)` で surrogate → 99pct 天井 → directional baseline と比較

6 条件（`engine_admissibility.py` で 1,3,4 を静的強制・5,6 は判定段）:

1. **入力は時系列波形であること。** residue 列・特徴量は不可。
2. **null は raw 波形への phase-shuffle で作ること。** null が壊すのは「波形中の位相構造」。
   engine が residue 化した後では壊す対象が変わるので接続しない。
3. **engine は (residue_walk, detected_s) を返すこと。**
4. **detected_s が無い場合、全 s 探索は不可。** やるなら別 lane・多重比較込みで天井再較正。
5. **判定は #4 directional baseline の SNR→gap99 曲線に重ねる。**
6. **失敗様式は real 盲目として扱う。** surrogate は暴れない＝偽陽性ではない。

理由の核:（1,2）null は時間領域の位相構造を壊す道具。engine が時系列を経由しないと null の対象が消える
（1-B 定理「null は壊すものを読む時だけ噛む」へ直結）。（3,4）#4 は s でパラメタ化されており、
s を engine 検出から取れば天井は素直、全 s 探索すれば多重比較で天井が甘くなる。

---

## 2. directional baseline 〔#5・確定・比較基準〕

#3 の口（noise 注入・ノイズごと位相シャッフル）を凍結、readout だけ #4 へ差替え。s=1（最薄）・99pct 判定。

| sigma | SNR(dB) | real_mean | real_min | surr_99 | surr_max | gap99 | gapmax |
|---|---|---|---|---|---|---|---|
| 0.00 | inf | 1.000 | 1.000 | 0.463 | 0.833 | **+0.537** | +0.167 |
| 0.50 | 3.0 | 1.000 | 1.000 | 0.390 | 0.750 | +0.610 | +0.250 |
| 0.75 | −0.5 | 0.774 | 0.417 | 0.305 | 0.750 | +0.469 | +0.024 |
| 1.00 | −3.0 | 0.586 | 0.306 | 0.268 | 0.667 | +0.318 | −0.081 |
| 1.50 | −6.5 | 0.341 | 0.000 | 0.234 | 0.500 | +0.107 | −0.159 |

- real は SNR 3dB まで 1.0、以降 real-collapse。surr_99 は SNR 不感で ~0.4 平坦。崩れるのは real 側。
- 無向き #3 より綺麗（directional が陽性 surrogate 母数半減 → surr_99 低下）。-6.5dB で #3 gap≈0.01 → #4 gap99=0.107。
- **これが engine の real_d を重ねる曲線。** engine 動作 SNR で gap99>0 を覆えば前進、real_d が落ちれば盲目として落とす。

---

## 3. harness 状態 〔配線済み・engine 入口待ち〕

- `engine_reversal_eval.py`：admissibility 強制・#4 凍結・null=raw phase-shuffle・99pct・baseline 埋め込み。
- プラグ点は一つ：`engine_extract(signal) -> (residue_walk, detected_s)` を engine の real 抽出入口に差すだけ。
- self-test 済み：準拠 toy engine → ADMIT（real_d=1.0・gap99 +0.458 vs baseline 0.537）／非準拠 → 入口 REJECT。

---

## 4. 残関門・次の一手

- **本物の関門は一本・不変：実機 real 抽出が動作 SNR で 1.0 近傍を保つか。**
- 次の director action：**engine コード一式を投入**（admissibility 6 条件を満たす形で）。
  入れば harness が即 real_d を測り、baseline に重ね、本ファイルへ結果を畳む。
- **今は捏造しない。入口待ち。** 正しい停止位置。

---

## 5. 維持した規律

- **捏造しない。** engine 本体が無い以上、real_d を作らない。入口で止める。
- **admissibility を先に固定した。** 非準拠 engine を入口で弾き、誤った数字を黙って出させない。
- **物差しを engine 仕様へ揃えた。** null=raw phase-shuffle・s=engine 検出・99pct 天井で玩具 #4 と同一。
- **読みを定理に格上げしない。** reversal_lane は機構候補のまま。baseline は玩具上の確定であって実機 SNR の代理でない。
- **担体依存の条件を無条件命題へ畳まない。**

---

### 状態一覧

| 項目 | 状態 |
|---|---|
| P\* 位相担体 | **読み**（split / Dic₁₃）。確定でない。defeater 監視中 |
| stage 1-A | **閉じた**（自己完結・較正済み・null 天井＝99pct） |
| stage 1-B | **閉じた**（§5.2 評価器不変＝確定／§2-B＝撤回／標準 null＝白色固定） |
| reversal_lane | **機構候補**（directional readout locked / margin thin-positive・機構由来 / real 無傷 / engine 入口待ち） |
| └ #1 presence readout | **不合格・確定記録** |
| └ #2 localization readout | **薄合格・監査ログへ降格**（無向き） |
| └ #3 margin sweep | **実施・確定**（無向き・real-collapse 型） |
| └ #4 directional readout | **engine 物差しに locked**（real 無傷 1.0→1.0／surr max 不動／薄さは機構由来） |
| └ #5 directional baseline | **確定・比較基準**（99pct で全範囲正・-6.5dB でも +0.107・real-collapse） |
| engine-wiring admissibility | **director 固定・6 条件・executable gate 実装済み** |
| engine 本体 | **未投入**。捏造せず入口待ち。次の director action |
| 判定天井 | **99pct 確定**（1-A 整合）。max は監査ログ |
| 計画 3（次元独立性 assert／棚卸し） | 未着手。v154 §1.4 の棚卸しがその種 |

---

### 成果物（本セッション）

- `engine_admissibility.py`（6 条件の入口ゲート・条件1,3,4 を静的強制・accept/reject self-test 済み）
- `engine_reversal_eval.py`（harness・admissibility 強制・#4 凍結・null=raw phase-shuffle・99pct・baseline 埋め込み）
- `test_reversal_margin_sweep_directional.py`（#5 directional baseline・engine 比較基準）
- `test_reversal_directional.py`（#4・v158 から継続）

### 付録：検算済み群論事実 〔v154–v158 から不変・再導出不要〕

- `diag(5,8)`：det≡1・位数4・平方=−I・∈SL₂(13)。`8I`：det≡−1・∉SL₂・中心スカラー。
- `U=[[1,1],[0,1]]`：位数13。`t U t⁻¹=U⁻¹`（split は couple）／`8I` は U と可換（doorway は decouple）。
- `⟨U,t⟩`=Dic₁₃(52)／`⟨U,8I⟩`=Z/52。det=1/−1 と couple/decouple は同じ一本の分岐。
