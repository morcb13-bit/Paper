# PROOF_CORE.md — 核心証明（最短版）
**対象**: 核62〜64 / 2026-03-19

---

## 1. transient = 573

```
573 = Z13_ORDER + H_ORDER×Z5_ORDER×SEQ_PERIOD
    = 13        + 4        × 5        × 28
```

**証明の骨格（4ステップ）:**

```
Step 1: prefix_A(0..12) mod BASE = BASE-STEP_Z13
  ∵ Sum_A(0..12) = 77, 77 mod 13 = 12, 77×240 mod 3120 = 2880 ✓

Step 2: t=13 での min_LSD_shift = 2880 + 17 = 2897  [=TRANSIENT_INIT_SHIFT]
  17 = floor(carry_in/cell) = floor((6-1/13)×3) = floor(18-3/13) = 17

Step 3: 28ステップごとに min_shift が +48 増加  [=TRANSIENT_STEP_PER_SEQ]
  48 = (BASE-STEP_Z13)/VERTICES

Step 4: アトラクタ入り条件 min_shift ≥ 736
  必要サイクル = ceil((736-2897+BASE) mod BASE / 48)
               = ceil(959/48) = 20 = H_ORDER×Z5_ORDER  ✓
  ∵ 20×28 = 4×(139+1) ≡ 4 = H_ORDER (mod ATTRACTOR_PERIOD)
```

---

## 2. 17 = Z13_ORDER + H_ORDER

```
init_shift = 2897 = (BASE-STEP_Z13) + 17

17 = floor((6 - 1/Z13_ORDER) × fanout)
   = floor(18 - 3/13)
   = 17

等しさ: 18-3/13 ≥ 17 ↔ Z13_ORDER ≥ fanout (13 ≥ 3 ✓)
        18-3/13 < 18  ↔ 3/13 > 0 ✓
```

---

## 3. ATTRACTOR_LSD_OFFSET = 738

```
total_macro_sum = net×VERTICES + CARRY_EXTRA
                = 325×60      + 52         = 19552

net         = total_add // BASE = 1016880 // 3120 = 325
CARRY_EXTRA = STEP_MIN = 52
  ∵ floor((BASE-STEP_Z13)×VERTICES/BASE) - fanout
   = floor(60×12/13) - 3 = 55-3 = 52

shift = net×fanout + CARRY_EXTRA×fanout/VERTICES - STEP_Z13
      = 975        + 52×3/60                    - 240
      = 975        + 2.6                         - 240
      = 737.6

round(737.6) = 738  ∵ 小数部0.6 ≥ 0.5
  0.6 = STEP_MIN×fanout mod VERTICES / VERTICES = 36/60
  36 ≥ 30  ∵ STEP_MIN×fanout mod VERTICES = VERTICES - 2×H_ORDER×fanout = 36 ✓
```

---

## 4. CARRY_EXTRA = STEP_MIN (幾何的等式)

```
STEP_MIN = VERTICES - Z5_ORDER - fanout
52       = 60       - 5        - 3

意味:
  VERTICES = 60:  切端20面体の頂点数
  Z5_ORDER = 5:   ceil(STEP_Z13/STEP_MIN) = ceil(240/52)
  fanout   = 3:   vertex degree
  STEP_MIN = 52:  活性状態数 = Z13_ORDER×H_ORDER = 13×4
```

---

## 脚注（実装不要・参照のみ）

- `17` が carry-geometric と period-algebraic の2経路で一致する詳細: v78追記2参照
- `20×28 ≡ H_ORDER (mod 139)` の証明詳細: 20×28 = 4×140 = 4×(139+1) ≡ 4
- `ATTRACTOR_MSD_HIGH_N=14` が `EARLY_CONVERGENCE_T=431=3×139+14` に再登場する理由: 保留(OPEN_QUESTIONS参照)
