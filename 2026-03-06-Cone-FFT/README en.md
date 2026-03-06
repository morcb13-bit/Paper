# cone-FFT / B13 Fractal Phase Library

**Inevitable Emergence of the Golden Ratio φ in BASE=3120 Integer FFT**

---

## Overview

This repository contains the mathematical foundations of cone-FFT, built on the B13 integer fractal phase representation.

Core result:

$$\mu^* = 2\,\mathrm{arcosh}\!\left(\frac{\sqrt{5}+\sqrt{13}}{4}\right) \approx 1.852, \qquad r(\mu^*) = -\varphi^{-2}$$

The golden ratio $\varphi = (1+\sqrt{5})/2$ is **not assumed** — it emerges inevitably from the structure of integer arithmetic.

---

## Philosophy

> Value lies in the remainder. Integer arithmetic only.

- Phase is manipulated using integer addition alone
- Trigonometric evaluation is secondary
- The choice BASE=3120 is not arbitrary

---

## Why BASE=3120?

$$3120 = \mathrm{lcm}(2^4,\; 15,\; 13) = 2 \times 5! \times 13$$

| Factor | Role |
|--------|------|
| $2^4 = 16$ | Sign period of local correlation $r_{\text{local}}$ |
| $15 = 3 \times 5$ | 5-fold symmetry × 3-level fractal hierarchy |
| $13$ | B13 lattice period (cone structure) |

---

## Main Results

### 1. Closed-Form Analytic Expression

$$\mu^* = 2\,\mathrm{arcosh}\!\left(\frac{\sqrt{5}+\sqrt{13}}{4}\right) \approx 1.852$$

### 2. Fixed Point

$$r(\mu^*) = -\varphi^{-2} \approx -0.38197$$

### 3. Causal Chain

```
BASE = lcm(2⁴, 15, 13) = 3120
        ↓
T/U parity split → U·n = T·n + n (mod N)
        ↓
Sign of r_local fully determined by n mod 16
        ↓
n=0, N/2 dominate a(n)  (DC and Nyquist concentration)
        ↓
ρ_a(μ) ≈ cosh(μ) / cosh(μ/2) = √5
        ↓
μ* = 2 arcosh((√5+√13)/4) ≈ 1.852
        ↓
r(μ*) = −φ⁻²   ←  φ emerges as a result
```

### 4. Diagonal Gate Fixed Point Theorem

The center $(7,7)$ of the rice-field problem (diagonal gate) connects to the T/U parity split through the double condition on $m=7$, forcing the $\sqrt{5}$ fixed point.

$$7 \equiv 2 \pmod{5} \;\Rightarrow\; \mathbb{Q}(\sqrt{5})$$
$$7 \equiv -6 \pmod{13} \;\Rightarrow\; \text{B13 boundary}$$

---

## Repository Structure

```
cone_fft_theory.md       # Paper (Japanese)
cone_fft_theory_en.md    # Paper (English)
README.md                # This file (Japanese)
README_en.md             # README (English)
```

---

## Quick Start

```python
import math

PHI = (1 + 5**0.5) / 2
BASE = N = 3120

# Analytic value of μ*
mu_star = 2 * math.acosh((5**0.5 + 13**0.5) / 4)
print(f"mu* = {mu_star:.6f}")  # → 1.852266

# Level-0 integer table
COS_TABLE = [round(BASE * math.cos(2*math.pi*k/BASE)) for k in range(BASE)]
SIN_TABLE = [round(BASE * math.sin(2*math.pi*k/BASE)) for k in range(BASE)]

# T/U split (all-even / all-odd)
T_list = [k for k in range(BASE) if k % 2 == 0]  # 1560 points
U_list = [k for k in range(BASE) if k % 2 == 1]  # 1560 points

# Sign of r_local (determined algebraically)
def r_local_sign(n):
    if n % 16 == 0: return +1   # image sets coincide
    if n % 8  == 0: return -1   # half-period shift
    return None                  # zero point
```

---

## Theory Status

| Item | Status |
|------|--------|
| Closed-form analytic expression for $\mu^*$ | ✅ Complete |
| Algebraic proof of $r = -\varphi^{-2}$ | ✅ Complete |
| Algebraic proof of $r_{\text{local}}$ sign | ✅ Complete |
| Necessity of BASE=3120 | ✅ Complete |
| Connection to the rice-field problem | ✅ Complete |
| Zero difference from fractal integer implementation | ✅ Complete |
| Paper draft | ✅ Complete (pre-review) |
| cone-FFT implementation | 🔲 Future work |
| Peer review submission | 🔲 Future work |

---

## License

TBD (pre-review draft)

---

## Changelog

| Version | Date | Content |
|---------|------|---------|
| v22 | 2026-03-04 | Numerical confirmation of μ* ≈ 1.852 |
| v23 | 2026-03-05 | Causal chain complete; BASE lcm structure |
| v24 | 2026-03-06 | Closed-form μ*; a(n) concentration; algebraic proof of r_local sign |
| v25 | 2026-03-06 | Cone eigenvalue proof; rice-field connection; theory complete |
