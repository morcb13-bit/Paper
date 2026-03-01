# φ-NTT: Additive Number Theoretic Transform

> An operator-theoretic framework connecting Fibonacci-based additive lattices to the spectral symmetry of the Riemann zeta function.

---

## Overview

The **φ-NTT (φ-Additive Number Theoretic Transform)** is an alternative transform framework where phase generation is governed by a linear recurrence relation ($F_{k+1} = F_k + F_{k-1}$) rather than the multiplicative cyclic groups of conventional NTTs.

This repository contains the theoretical manuscript and Python implementation of φ-NTT, developed as a formal bridge between discrete additive lattice dynamics and analytic number theory.

---

## Key Features

- **Carry-free arithmetic** over $\mathbb{Z}_{10}^B$ — no carries between digits
- **Golden ratio integer ring** $\mathbb{Z}[\varphi]$ — integer arithmetic only, no floating point
- **B-stage tensor product** — $2^B$ channels with zero inter-stage twiddle factors
- **Verified** for $B = 1, 2, 3, 4, 5$ — delta convolution and round-trip tests pass
- **B = 13 threshold** — structural observations suggesting spectral stabilization

---

## Repository Structure

```
.
├── phi_carry_free.py          # Core implementation
├── phi-ntt-complete-en.md     # Full manuscript (English)
├── phi-ntt-complete-ja.md     # Full manuscript (Japanese)
└── README.md                  # This file
```

---

## Quick Start

```python
from phi_carry_free import phi_conv_carry_free
from random import seed, randint

seed(42)
B, N = 3, 1000
x = [(randint(-5, 5), randint(-1, 1)) for _ in range(N)]
h = [(0, 0)] * N
h[0] = (1, 0)
y = phi_conv_carry_free(x, h, B)

print("delta conv B=3:", "✓" if all(y[i] == x[i] for i in range(N)) else "✗")
```

---

## The Core Idea

Traditional NTTs rely on multiplicative structure. φ-NTT replaces this with an **additive recursive kernel**:

$$K_\phi(n, k) = \omega_N^{\,n \cdot F_k}$$

where $F_k$ is the $k$-th Fibonacci number modulo $N$.

Within this framework, a canonical **Invariant Manifold** $\mathcal{M}$ emerges as the fixed-point set of a dual anti-linear involution $\mathcal{J}$. Its asymptotic behavior as $N \to \infty$ is conjectured to correspond formally to the critical line $\operatorname{Re}(s) = \tfrac{1}{2}$ of the Riemann zeta function.

This is a **speculative structural observation**, not a proof of the Riemann Hypothesis.

---

## Status

| Component | Status |
|-----------|--------|
| Core implementation (`phi_carry_free.py`) | ✅ Complete |
| Verification B = 1..5 | ✅ Passed |
| Manuscript (EN / JA) | ✅ Complete |
| Rigorous proof of $\mathcal{M}$ non-triviality | 🔬 Open problem |
| Connection to analytic $L$-functions | 🔬 Future work |

---

## Authors & Collaboration

This project is a three-way collaboration:

- **Theory & design**: ChatGPT (カクカタリキ)
- **Implementation & manuscript**: Claude (Anthropic)
- **Direction & vision**: Project lead

---

## License

MIT License

---

*For questions or collaboration inquiries, please open an Issue.*
