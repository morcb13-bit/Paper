# B13 引継書 v158

（v157 から継ぐ。v157 は reversal_lane を **機構候補**とタグ付けし engine 配線を authorize、
本物の関門を「実機 real 抽出が動作 SNR で 1.0 近傍を保つか」一本に絞った。
本セッションは **engine に繋ぐ前にロック定義との整合を一回だけ修復**した——#2 の無向き readout を
**方向忠実 readout（#4）**へ差し替え、実機 fixture で margin が開くかを問うた。
結論：**real は無傷（1.0→1.0）／surr max はほぼ不動／#2 の薄さは機構由来であって readout 由来ではない。**
directional readout を engine の凍結 instrument として locked。#1/#2/#3 は凍結・上書きせず #4 を別ファイルで追加。）

---

## 0. このセッションの一行

ロック定義は「stable run s に *続いて* -s mod 13」＝向きを含む。#2 が読んでいたのは無向き antipodal adjacency。
そこで **#4 directional readout** を切り、対蹠判定だけを `va==s かつ vb==(-s)%13`（s 先・-s 後のみ credit）に差替え。
**real は傷つかず（1.0→1.0・全 s・全 seed）／surr max の drop は 0.00〜0.08（しばしば 0.00）。**
**#2 の薄さは「無向き readout の弱さ」ではなく「機構の薄さ」と確定。** directional を engine 物差しとして locked。

---

## 1. 検証の構造 〔#2 から readout 一本だけ差替え〕

- fixture（`make_signal`・`phase_shuffle`・`extract_residue`・`stable_runs`）は #1/#2 から **import で共有**。
  **差し替えたのは対蹠判定の一行だけ**：`vb==(-va)%13`（無向き）→ `va==s and vb==(-s)%13`（方向忠実）。
- 隣接は stable_runs 上の連続＝#2 と同じ「許容ギャップ（min_run=3 で落ちた短ラン）」を含む。新しい tuning は無し。
- 規律: surrogate（位相シャッフル）・null 天井（99pct）・n_surr（400）は **#2 と同一固定**。動かすのは readout のみ。
  pass させるために surrogate・閾値を弄らない。出た数字をそのまま読む。

---

## 2. 実機 anchor の再現 〔噛み合わせ確認・確定〕

- #1 presence readout：s=1 で surr max=1.25>real（漏れ）、real_walk=`[1]×12+[12]×12`。v157 と一致。
- #2 localization readout：real=1.0、surr max **s1=0.75 / s2=0.33 / s3=0.33**、最悪 surr `[(12,10),(1,8)]`。v157 と一致。
- → 私（calculator）が実機の口に正しく噛んでいることを確認。#4 の数字は lock-grade。

---

## 3. #4 directional readout 〔実施・確定・engine 物差し〕

実機 fixture・#2 と同一アンサンブル・readout 一行だけ差替え。

| s | real_u | real_d | surr_max u→d | surr_99 u→d | margin(99) u→d | margin(max) u→d |
|---|---|---|---|---|---|---|
| 1（最薄） | 1.000 | **1.000** | 0.750→0.708 | 0.583→0.542 | +0.417→+0.458 | +0.250→+0.292 |
| 2 | 1.000 | **1.000** | 0.333→0.333 | 0.0025→0.000 | +0.998→+1.000 | +0.667→+0.667 |
| 3 | 1.000 | **1.000** | 0.333→0.292 | 0.292→0.250 | +0.708→+0.750 | +0.667→+0.708 |

- seed 0–5（s=1）で surr max の drop は **0.00〜0.083（しばしば 0.00）**。99pct の drop は ~0.04〜0.06 で一貫。
- s=1 の新最悪 surr は `[(1,8),(12,9)]`＝**正順・被覆 0.708**。文書化の `[(12,10),(1,8)]`（逆順）は殺せたが、
  位相シャッフルは順序対称なので**正順の高被覆 surrogate が同程度の高さで必ず居る**。

所見:
- **real 無傷（確定）**: real の反転は構成上 `s→-s` の正順。向きを要求しても 1.0 を保つ。方向化は real にコストゼロ。
- **surr max ほぼ不動（確定）**: #2 の薄さ（s=1 max 0.75）は **readout 由来ではなく機構由来**。
  文書化の最悪例が逆順だったのは半ば偶然で、被覆分布に正順の高被覆メンバーが居る。
- **99pct は僅かに良化**: 方向化で陽性 surrogate 母数が半減（frac>0 が 0.080→0.043）するぶん 99pct が ~0.04 下がる。

---

## 4. 判定 〔v157 の二分岐の後者・ただし物差しは格上げ〕

- v157 の分岐「方向を入れて margin が開けば reversal_lane はきれいに／開かなければ機構候補の薄さを正直に受け止め engine へ」。
- **結果は後者。margin は大きく開かない。** だが readout は **向き忠実（定義整合）へ格上げされ、real にコスト無し・99pct 僅か良化**。
- → **directional_localization_score を engine の凍結 instrument として locked。** 無向き #2 は監査ログとして残す。
- **本物の関門は不変・一本**：実機 real 抽出が動作 SNR で 1.0 近傍を保つか。
  保てば margin は正（99pct で約 2×）。落ちれば「偽陽性」でなく「盲目」として落とす。

---

## 5. 維持した規律

- **engine 前にロック定義整合を一回だけ修復した。** 「対蹠が隣接」は reversal ではない。向きを入れた。
- **#1/#2/#3 は凍結・上書きせず、#4 を別ファイルで追加。** 較正手順・surrogate・閾値は #2 と同一。
- **緑への最適化をしない。** 動かしたのは readout の対蹠判定一行だけ。margin が開かなかった事実をそのまま記録。
- **読みを定理に格上げしない。** 機構候補のまま。「薄さは機構由来」も玩具上の確定であって実機 SNR の代理ではない。
- **担体依存の条件を無条件命題へ畳まない。** directional の優位は「real にコスト無し＋99pct 僅か良化＋定義整合」に限定。

---

### 状態一覧

| 項目 | 状態 |
|---|---|
| P\* 位相担体 | **読み**（split / Dic₁₃）。確定でない。defeater 監視中 |
| stage 1-A | **閉じた**（自己完結・較正済み・null 天井＝99pct） |
| stage 1-B | **閉じた**（§5.2 評価器不変＝確定／§2-B＝撤回／標準 null＝白色固定） |
| reversal_lane | **機構候補**（toy confirmed / directional readout locked / margin thin-positive・**機構由来** / real 無傷 / engine SNR defeater pending） |
| └ #1 presence readout | **不合格・確定記録**（`|X|` 由来 presence を漏らす） |
| └ #2 localization readout | **薄合格・確定／監査ログへ降格**（無向き・被覆一特徴・s=1 最薄） |
| └ #3 margin sweep | **実施・確定**（real-collapse 型・99pct で SNR −3dB まで正） |
| └ #4 directional readout | **実施・確定・engine 物差しに locked**（real 無傷 1.0→1.0／surr max 不動／薄さは機構由来） |
| 判定天井 | **99pct 確定**（1-A 整合）。max は監査ログ |
| 2-return（二本 lane 弁別） | reversal_lane を engine 実機化して問う段。関門＝real 抽出 SNR |
| 計画 3（次元独立性 assert／棚卸し） | 未着手。v154 §1.4 の棚卸しがその種 |

---

### 成果物（本セッション）

- `test_reversal_directional.py`（#4 directional readout・#1/#2 から fixture/stable_runs を import し対蹠判定だけ差替え）
- **検証ログ**：実機 anchor 再現（#1 max1.25・#2 max 0.75/0.33/0.33）／#4 結果（real 1.0→1.0・s=1 surr max 0.750→0.708・
  新最悪 surr `[(1,8),(12,9)]` 正順 cov0.708）／seed 0–5 で surr max drop 0.00〜0.08

### 付録：検算済み群論事実 〔v154–v157 から不変・再導出不要〕

- `diag(5,8)`：det≡1・位数4・平方=−I・∈SL₂(13)。`8I`：det≡−1・∉SL₂・中心スカラー。
- `U=[[1,1],[0,1]]`：位数13。`t U t⁻¹=U⁻¹`（split は couple）／`8I` は U と可換（doorway は decouple）。
- `⟨U,t⟩`=Dic₁₃(52)／`⟨U,8I⟩`=Z/52。det=1/−1 と couple/decouple は同じ一本の分岐。
