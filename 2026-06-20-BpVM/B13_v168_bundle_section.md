# B13 引継書 v168 — bundle セクション（reversal_lane real_d 準備）

**target event**: `20260625073000` 岩手県沖 M6.9（速報・深さ50km）
**震源（波形ヘッダ）**: 40.2°N, 142.3°E, 50km
**状態**: real_d **保留中**。K=3 **全点 gate-only 完了・全 admissible**。次は prior 再掲 → real_d（director の手番）。

---

## prior 宣言（測定前に固定・据え置き）

> 地震コーダの典型 SNR では **ABSTAIN / 全BLIND（real-collapse）が既定の予想**。
> FORWARD が出ても、a_prep セルフテスト経路4（不変 decoy → 偽 FORWARD）が示す通り
> **prior 依存の偽通過**をまず疑う。「反転を見つけた」と読みすぎない。
> （a_prep PRIOR_DECLARATION: BUNDLE_FORWARD は真の保証ではなく cross-signal stable dominant s の観測。）

---

## K=3 bundle 名簿

| # | station | net | 成分 | 地表/地中 | fs | 長さ | 震央距離 | 震源距離 | gate条1 | gate条3,4 | 状態 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | MYG001 | K-NET | NS | 地表 | 100Hz | 164s | 157.3km | 165.1km | PASS | PASS (s=4※) | **確定** |
| 2 | IWTH03 | KiK-net | **ch4**=NS | **地表** | 100Hz | 221s | 70.7km | 86.6km | PASS | PASS (s=2※) | **確定** |
| 3 | IWT003 | K-NET | NS | 地表 | 100Hz | 185s | 41km | 64.6km | PASS | PASS (s=6※) | **確定** |

※ detected_s は読み出しのみ・発見ではない。real_d 未測定。
**K=3 全点 admissible 確認済（gate 完成）。次は prior 再掲 → real_d（director の手番）。**

---

## 成分の規律（重要）

- **IWTH03（KiK-net）は6成分**。地表/地中は振幅で同定した：
  - **ch4–6 = 地表**（peak 200–246 gal、表層増幅）→ real_d は **ch4（地表NS）**を使う。
  - ch1–3 = 地中（peak 30–45 gal）→ **除外**。地表と混ぜると解釈が濁る。
  - ※初回ゲートで誤って ch1（地中）を使用 → ch4 に訂正済。
- **MYG001 / IWT003（K-NET）は地表単一**。NS 列をそのまま採用、成分の濁りなし。
- 全点 100Hz・raw 加速度（gal）・demean（DC除去）後に admissibility へ。

---

## 同地点検証候補（別枠保存・bundle には入れない）

- **IWT019**（K-NET, 39.803°N 141.658°E, 震央70km）は **IWTH03 とほぼ同位置**。
- 用途: 将来「同地点の地表 K-NET vs 地表/地中 KiK-net」の比較実験。
- 今回 bundle には入れない（同地点だと bundle の独立性が濁るため。director 判断）。

---

## 進行手順（残）

1. ~~IWT003 NS 生波形を受領~~ **済**（K-NET csv、185s）。
2. ~~3点まとめて gate 取り直し~~ **済**（条1・条3,4 とも）。
3. ~~全点 admissible 確認~~ **済 → K=3 gate 完成**。
4. **prior を再掲**（上記）。← 次セッション冒頭、または本日続行時
5. ここで初めて **real_d 測定**（`run_bundle([MYG001, IWTH03_ch4, IWT003], target="20260625073000")`）。director の手番。
6. 出力解釈: BUNDLE_FORWARD でも prior 偽通過を疑う／ABSTAIN・BLIND は honest output として記録（失敗ではない）。

---

## 不変条件（崩さない）

- gate と real_d を混ぜない。gate は型・admissibility のみ、real_d は測定。
- K=2 での real_d は読まない（K_EXPECT=3）。3点揃ってから。
- engine 同定は名前でなく sha256（canonical `engine_reversal_eval.py` = 5d0239795333）。
- ABSTAIN は honest output。FORWARD を成果として急がない。
