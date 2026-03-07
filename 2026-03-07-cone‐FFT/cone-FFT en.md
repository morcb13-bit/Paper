# cone-FFT

An integer two-channel transform with local cone geometry.

**v29 — 2026-03-07**  
Theory complete. Full invertible pipeline pending (see Status).

---

## What is this?

cone-FFT is a signal transform built on two independent layers:

1. **Parseval core** — A 2-channel integer recursion with exact norm scaling factor 40
2. **Cone geometry** — Local geometric invariants derived from a (3,4) rotation orbit, yielding the golden ratio φ as a *result*, not an assumption

The two layers share exactly one structural link: the prime factor **5**.

---

## Quick Start

```bash
# Run verification (all propositions)
python demo_verify.py

# Run full pipeline demo
python cone_fft_complete.py
```

Requirements: Python 3.10+, standard library only.

---

## Core Results

### The recursion matrix

```python
M = [[6, 2],
     [2,-6]]
```

satisfies **M² = 40I**, giving exact integer norm scaling per step.

```python
from parseval_core import step_rule_40, step_rule_40_inv

# Forward: norm ×40
T1, U1 = step_rule_40(T0, U0)

# Inverse: norm ÷40  (exact integer division)
T0r, U0r = step_rule_40_inv(T1, U1)
assert T0r == T0  # exact recovery
```

### The cone chain

From the 4-point orbit of (3,4):

```
G = ((25,-25),(-25,50))  →  disc = 5⁵  →  ρ_eig = φ²  →  ρ_a = √5
→  R = -φ⁻²  →  C_CONE = 2/φ
```

φ is derived, not assumed.

```python
from cone_local import cone_chain
ch = cone_chain()
# ch['R']      == -0.38196601  (-φ⁻²)
# ch['C_CONE'] ==  1.23606798  (2/φ)
```

---

## Propositions

| # | Statement | Status |
|---|---|---|
| A | ‖T'‖²+‖U'‖² = 40(‖T‖²+‖U‖²) | ✓ proved (M²=40I) |
| B | Gram disc = 5⁵ → φ² → √5 → −φ⁻² → 2/φ | ✓ proved |
| C | Integer invertibility: (6T'+2U')/40 = T | ✓ proved |
| D | 40 determined by 4 independent conditions | ✓ proved |

---

## BASE = 3120

```
3120 = lcm(3, 16, 5, 13)

factor | meaning
-------|--------
3      | balanced ternary base
16     | r_local sign period
5      | √5 geometry (gateway to φ)
13     | B13 lattice
```

Minimum common period for balanced ternary depth 4,
which is the minimum depth at which B13 × √5 coexist.

---

## File Structure

```
parseval_core.py      Integer 2-channel recursion + inverse
cone_local.py         Local Gram chain → C_CONE = 2/φ
cone_fft_complete.py  End-to-end pipeline
demo_verify.py        Full verification suite
PAPER.md              Theory paper (English)
PAPER_ja.md           Theory paper (Japanese)
```

---

## Status

| Component | State |
|---|---|
| Parseval core (40-recursion) | ✅ Complete |
| Cone geometry (φ-chain) | ✅ Complete |
| Integer forward transform | ✅ Complete |
| Integer inverse transform | ✅ Complete |
| cos+sin Parseval identity | ✅ Verified |
| BASE=3120 construction | ✅ Complete |
| **cos+sin full invertible pipeline** | 🔲 Pending |

The single remaining gap: adding the sin channel to close the
end-to-end integer FFT system. All algebraic foundations are in place.

---

## License

MIT

---

*cone-FFT v29*
