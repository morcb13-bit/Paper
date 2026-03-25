# SPEC_MIN.md — coneFFT 実装仕様（最小版）
**対象**: b13phase v0.7.9 / 2026-03-19

---

## 定数一覧

```python
# === 既存定数（変更なし）===
BASE     = 3120   # = 60×52
STEP_Z13 = 240    # = BASE/13
STEP_MIN = 52     # = BASE/60 = ACTIVE
VERTICES = 60
H_ORDER  = 4
Z5_ORDER = 5
Z13_ORDER = 13
SEQ_PERIOD = 28   # Pisano周期 π(13)
ATTRACTOR_PERIOD  = 139
ATTRACTOR_LSD_OFFSET  = 738  # = 14×52+10
ATTRACTOR_MSD_HIGH_N  = 14
ATTRACTOR_MSD_LOW_N   = 46
ATTRACTOR_LSD_MIN     = 10

# === v0.7.9 新定数 ===
TRANSIENT_STD          = 573   # 標準初期条件でのtransient長
TRANSIENT_SEQ_CYCLES   = 20    # = H_ORDER × Z5_ORDER
TRANSIENT_STEP_PER_SEQ = 48    # = (BASE-STEP_Z13) / VERTICES
TRANSIENT_INIT_SHIFT   = 2897  # t=Z13_ORDERでのmin_LSD_shift
CARRY_EXTRA            = 52    # = STEP_MIN = total_macro_sum mod VERTICES
EARLY_CONVERGENCE_T    = 431   # = 3×ATTRACTOR_PERIOD + ATTRACTOR_MSD_HIGH_N
```

---

## 各定数の最短定義式

| 定数 | 値 | 定義式 |
|------|----|--------|
| `TRANSIENT_STD` | 573 | `Z13_ORDER + H_ORDER×Z5_ORDER×SEQ_PERIOD` = 13+20×28 |
| `TRANSIENT_SEQ_CYCLES` | 20 | `H_ORDER × Z5_ORDER` = 4×5 |
| `TRANSIENT_STEP_PER_SEQ` | 48 | `(BASE-STEP_Z13) / VERTICES` = 2880/60 |
| `TRANSIENT_INIT_SHIFT` | 2897 | `(BASE-STEP_Z13) + 17` = 2880+17 |
| `CARRY_EXTRA` | 52 | `STEP_MIN` = `VERTICES-Z5_ORDER-fanout` = 60-5-3 |
| `EARLY_CONVERGENCE_T` | 431 | `3×ATTRACTOR_PERIOD + ATTRACTOR_MSD_HIGH_N` = 3×139+14 |

---

## 依存関係

```
BASE, VERTICES, STEP_Z13, H_ORDER, Z5_ORDER, Z13_ORDER, SEQ_PERIOD
  └─→ TRANSIENT_STEP_PER_SEQ = (BASE-STEP_Z13)/VERTICES
  └─→ TRANSIENT_SEQ_CYCLES   = H_ORDER × Z5_ORDER
  └─→ TRANSIENT_STD          = Z13_ORDER + TRANSIENT_SEQ_CYCLES × SEQ_PERIOD

ATTRACTOR_LSD_OFFSET (=738), STEP_MIN (=52)
  └─→ ATTRACTOR_MSD_HIGH_N = 738 // 52
  └─→ ATTRACTOR_LSD_MIN    = 738 mod 52

ATTRACTOR_PERIOD (=139), ATTRACTOR_MSD_HIGH_N (=14)
  └─→ EARLY_CONVERGENCE_T = 3×139 + 14
```

---

## テスト項目

```python
# 1. transient 測定
eng = FibDriveEngine(mode='A', n_levels=2)
assert eng.attractor_info()['transient'] == TRANSIENT_STD  # 573

# 2. アトラクタLSD格子
eng.reset()
for _ in range(TRANSIENT_STD): eng.step()
lsds = sorted(eng.state_hash())
assert lsds[0] == ATTRACTOR_LSD_MIN       # 10
gaps = [lsds[i+1]-lsds[i] for i in range(59)]
assert all(abs(g - ATTRACTOR_LSD_STEP) <= 2 for g in gaps)  # ±2許容

# 3. 二段階収束
# (実装の必要性に応じて)
eng2 = FibDriveEngine(mode='A', n_levels=2)
# t=EARLY_CONVERGENCE_T(431) で上位20セルがshift>=736 → 実測確認済み
```

---

## 変更ファイル (v0.7.8 → v0.7.9)

- `fib_drive.py`: 上記新定数を追加
- `__init__.py`: エクスポート追加、`__version__ = "0.7.9"`
