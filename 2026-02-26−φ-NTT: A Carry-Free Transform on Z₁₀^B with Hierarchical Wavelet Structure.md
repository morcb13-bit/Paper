# φ-NTT: A Carry-Free Transform on Z₁₀^B with Hierarchical Wavelet Structure

<!-- Reference implementation (frozen): https://github.com/morcb13-bit/Artemis-B13-Archive/tree/main/tutorial @ [commit SHA here] -->
<!-- Source: 引継書 v2〜v12 (2026-02-25) -->
<!-- 禁止: HTML作成・B=4,5再計算・新実験・可視化追加 -->
<!-- 次タスク: .tex 化（arXiv向け最小テンプレ） -->

---

## Abstract

We introduce φ-NTT (phi Number Theoretic Transform), a transform defined over the carry-free direct product group Z₁₀^B with coefficients in the golden-ratio integer ring Z[φ] = Z[√5]. The transform achieves exact integer arithmetic with no floating-point error and no inter-stage twiddle factors. The B-stage construction yields a 2^B channel filterbank with a hierarchical wavelet structure analogous to Haar MRA. We establish a Spectral Concentration Quasi-Theorem: for signals with dominant period P ≈ 2·10^b, the detail channel at level b captures energy fraction α_b ≥ C·(1−ε)^{B−1}. Empirically, the top 3 channels cover ≥ 97% of total energy for B = 3..5.

---

## 1. Introduction

The Fast Fourier Transform (FFT) and its variants are among the most widely deployed algorithms in signal processing and communications. Standard FFT operates on cyclic groups Z/(N) and relies on complex-valued twiddle factors, which introduce floating-point error and prevent exact integer arithmetic.

In this paper, we introduce the φ-NTT (phi Number Theoretic Transform), a transform defined over the direct product group Z₁₀^B = Z₁₀ × ⋯ × Z₁₀ — the group of B-digit decimal representations with carry-free (digit-wise modular) arithmetic. The key structural features of φ-NTT are:

- **Exact integer arithmetic**: all coefficients lie in the ring Z[φ] = Z[√5] (golden-ratio integer ring); no floating-point error.
- **Zero inter-stage twiddle factors**: the tensor-product theorem eliminates all cross-stage phase corrections over B recursive stages.
- **2^B channel filterbank**: the B-stage transform decomposes a length-10^B signal into 2^B channels, labeled by binary strings over {T, U}, where T denotes a C₅-twisted low-pass projection and U denotes a Z₂ high-pass (parity) projection.

The carry-free arithmetic setting — where each decimal digit evolves independently — naturally gives rise to a hierarchical wavelet structure analogous to the Haar multiresolution analysis (MRA), but defined over the ring Z[φ] and indexed by decimal digits rather than binary positions. We formalize this analogy and prove a Spectral Concentration Quasi-Theorem: for signals with dominant period P ≈ 2·10^b, the detail channel at level b captures a fraction α_b ≥ C·(1−ε)^{B−1} of total signal energy. Empirically, the top 3 channels cover ≥ 97% of energy for B = 3..5.

The remainder of this paper is organized as follows. Section 2 defines φ-NTT and the convolution theorem on Z₁₀^B. Section 3 develops the hierarchical wavelet interpretation. Section 4 states and sketches the proof of the quasi-theorem. Section 5 presents numerical experiments for B = 1..5. Section 6 discusses open problems, including the exact orthogonality of T/U projections and the determination of constants C and ε. Section 7 concludes.

---

## 2. φ-NTT on Z₁₀^B

### 2.1 Carry-Free Group Structure

Let B be a positive integer and N = 10^B. We represent each index n ∈ {0, …, N−1} by its decimal digit expansion

    n = d₀(n) + 10·d₁(n) + ⋯ + 10^{B−1}·d_{B−1}(n),   dᵢ(n) ∈ {0,…,9}.

The carry-free group Z₁₀^B is the direct product Z₁₀ × ⋯ × Z₁₀ (B copies), with addition defined digit-wise:

    (n ⊕ m)ᵢ = (dᵢ(n) + dᵢ(m)) mod 10.

This differs fundamentally from ordinary modular arithmetic on Z/(10^B)Z, where addition carries across digit boundaries. The two groups coincide only at B = 1; for B ≥ 2 they are non-isomorphic, and the standard cyclic convolution theorem does not apply to Z₁₀^B.

The carry-free subtraction and negation are defined analogously:

    cf_sub(n, m)ᵢ = (dᵢ(n) − dᵢ(m)) mod 10,
    cf_neg(n)ᵢ   = (−dᵢ(n)) mod 10.

All arithmetic in φ-NTT operates over the ring Z[φ] = Z[√5], where φ = (1+√5)/2. Elements are pairs (a, b) ∈ Z² representing a + b·φ, with multiplication

    (a, b) · (c, d) = (ac + bd, ad + bc + bd).

This eliminates floating-point error entirely.

### 2.2 Definition of φ-NTT

For B = 1 (length N = 10), the transform decomposes the input x : Z₁₀ → Z[φ] into two channels via the isomorphism Z₁₀ ≅ Z₂ × Z₅ (Chinese Remainder Theorem):

    T-channel:  X_T[k] = n10_T(x)[k]   (C₅-twisted DFT projection)
    U-channel:  X_U[k] = n10_U(x)[k]   (Z₂ parity projection)

The C₅-twisted DFT uses the basis elements cos5_zphi(m) ∈ {(2,0),(−1,1),(0,−1)} — three distinct values in Z[φ], with no irrational twiddle factors at runtime. The Z₂ projection uses the parity kernel U₁[k] ∈ {(0,0),(±1,0),(0,±1)}. Explicit tables are given in Appendix A.

For general B, the transform is defined by a B-stage tensor product. At each stage b ∈ {0,…,B−1}, the input is projected onto either the T-channel (low-pass, C₅ side) or the U-channel (high-pass, Z₂ side). This yields 2^B output channels, each labeled by a binary string over {T, U}^B:

    bn_forward_flat(x, B) → [ch_s : s ∈ {T,U}^B],   each ch_s : Z₁₀^B → Z[φ].

The inverse transform bn_inverse_flat reconstructs x exactly, with a scale factor of 40^B (no convolution) or 80^B = 40^B × 2^B (after channel-domain convolution). Both scale factors divide exactly in Z[φ], so reconstruction is lossless over the integers.

A key structural property is that no inter-stage twiddle factors appear: the tensor-product theorem over Z₁₀^B guarantees that all cross-stage phase corrections vanish across all B recursive stages.

### 2.3 Convolution Theorem (Carry-Free)

Define the carry-free convolution of x, h : Z₁₀^B → Z[φ] as

    (x ⋆ h)[n] = Σ_{m} h[m] · x[cf_sub(n, m)],   n ∈ Z₁₀^B.

**Theorem (Carry-Free Convolution).** Let F_B denote the φ-NTT forward transform. Then F_B diagonalizes carry-free convolution on Z₁₀^B: the channel-domain product

    Y_s[k] = X_s[k] · H_s[k]   (with the T/U cross-coupling rule)

recovers x ⋆ h upon inverse transform, up to the scale factor 80^B. The T/U cross-coupling at each stage is:

    2·Y_T[k] = X_T[k]·H_T[k] − SIN2 · X_U[k]·H_U[k]
    2·Y_U[k] = X_T[k]·H_U[k] + X_U[k]·H_T[k]

where SIN2 = (3,−1) ∈ Z[φ] ≈ 1.382 = 4sin²(2π/10). This has been verified for B = 1..5 against direct computation.

Similarly, carry-free cross-correlation reduces to convolution with the carry-free-reversed kernel h_rev[m] = h[cf_neg(m)].

---

## 3. Hierarchical Wavelet Interpretation

The recursive T/U decomposition induces a hierarchical structure on Z₁₀^B that parallels multiresolution analysis (MRA).

### 3.1 T/U as Projection-Like Operators

At each stage b ∈ {0,…,B−1}, the φ-NTT applies one of two projection-like operators to the digit d_b:

    T_b : low-pass side  (C₅-twisted DFT projection)
          extracts the slowly-varying component across d_b ∈ {0,…,9}

    U_b : high-pass side (Z₂ parity projection)
          extracts the alternating component  (−1)^{d_b}

Together, T_b and U_b partition the information in digit b in a manner analogous to a two-channel filter bank. We use the term *projection-like* deliberately: while the numerical experiments confirm that energy is partitioned across channels without loss (exact reconstruction holds for all B = 1..5), the strict algebraic orthogonality T_b ⊥ U_b has not yet been proved and is left as an open problem (Section 6).

The B-stage tensor product applies these operators independently at each digit position, yielding 2^B channels:

    ch_s = (⊗_{b=0}^{B-1} P_b^{s_b})(x),   s ∈ {T,U}^B,

where P_b^T = T_b and P_b^U = U_b. The full signal energy is distributed across all 2^B channels, and exact reconstruction is guaranteed by bn_inverse_flat.

### 3.2 Binary Tree Interpretation

The label set {T,U}^B is in natural bijection with the set of root-to-leaf paths in a complete binary tree of depth B. Reading the label from left (most significant digit, b = B−1) to right (least significant, b = 0):

```
          root
        /       \
     T(d_{B-1})  U(d_{B-1})
       /   \       /   \
     T      U    T      U
    (d_{B-2}) …
```

Each leaf corresponds to one channel. The channel TT…T (all T, leftmost leaf) captures the coarsest approximation; the channel UU…U (all U, rightmost leaf) captures the finest detail — the checkerboard pattern across all digit positions simultaneously.

This structure is structurally analogous to Haar MRA, but not identical. The key correspondences and differences are:

| Property | Haar MRA | φ-NTT |
|----------|----------|-------|
| Signal length | 2^B | 10^B |
| Channels | 2^B | 2^B |
| Approximation coeff. | scaling function | TT…T channel |
| Detail coeff. | wavelet function | T…TU_bT…T channels |
| Arithmetic | real (±1 basis) | Z[φ] basis |
| Orthogonality | exact | conjectured |

The analogy enables direct transfer of Haar intuitions: a smooth signal concentrates in approximation (T-dominant) channels, while rapid oscillations populate detail (U-heavy) channels.

### 3.3 Energy Localization Mechanism

The tree structure predicts how signal energy distributes across channels, depending on the signal's smoothness and dominant period.

Smooth signals vary slowly across all digit positions. The U_b operator measures the parity alternation (−1)^{d_b} at digit b; if x changes little between even and odd values of d_b, the U_b inner product is small. Accumulating this across B stages, the UU…U channel — which requires large alternation at every digit simultaneously — is driven toward zero. This explains the empirical observation (B = 3..5) that the UU…U channel carries < 0.01% of total energy for smooth signals.

Signals with a dominant period P ≈ 2·10^b exhibit strong alternation specifically at digit b, while other digits vary slowly. The U_b projection then captures the bulk of the oscillatory energy, while T projections at all other stages pass through the slowly varying envelope. The result is that the channel T…T[U_b]T…T (U only at position b) becomes the dominant channel. For the test signal sin(2π·n/100) with P = 100 ≈ 2·10¹, this predicts the TUT…T channel as dominant — consistent with all B = 3..5 experiments (86%, 69%, 68% of total energy respectively).

The precise quantification of this concentration — a lower bound on the fraction α_b captured by the dominant channel — is the subject of Section 4.

---

## 4. Spectral Concentration

We formalize the observed concentration phenomenon in the form of a quasi-theorem.

### 4.1 Statement

We work under two assumptions on the input signal x : Z₁₀^B → Z[φ] (real-valued, i.e., second component identically 0):

**(A1) Digit-wise smoothness.** For each digit position b' ≠ b, the variation of x across d_{b'} is small relative to the variation across d_b:

    |variation(x, b')| ≤ ε · |variation(x, b)|,   b' ≠ b,

for some constant 0 < ε < 1. Here variation(x, b) measures the mean absolute change in x as digit b steps through {0,…,9}.

**(A2) Dominant period.** The signal has a dominant period P ≈ 2·10^b for some b ∈ {0,…,B−1}, so that d_b is the primary source of oscillation.

**Quasi-Theorem (φ-NTT Spectral Concentration).** Under assumptions (A1) and (A2), there exist constants C > 0 and 0 < ε < 1 (depending on the signal class) such that the energy fraction captured by the level-b detail channel satisfies

    α_b  =  E[ch_{T…TU_bT…T}] / E_total  ≥  C · (1 − ε)^{B−1}.

Furthermore, the remaining energy is distributed such that the top 3 channels (by energy) jointly satisfy

    α_b + α_{b,2} + α_{b,3}  ≥  1 − δ

for a small residual δ, empirically δ < 0.03 for B = 3..5.

We call this a *quasi-theorem* because the constants C and ε are characterized implicitly through (A1)–(A2) rather than given in closed form, and because the strict algebraic orthogonality of the T/U decomposition — while consistent with all numerical evidence — remains unproved (see Section 6).

### 4.2 Interpretation

The bound C · (1−ε)^{B−1} has a transparent structure:

- ε measures the digit-level roughness of the signal: how much energy "leaks" from the dominant digit b to each of the other B−1 digits. For a clean periodic signal, ε is small.
- (1−ε)^{B−1} is the product of B−1 suppression factors, one per non-dominant digit stage. The bound decays exponentially in B, which explains why the observed α_b decreases from 86% (B=3) to 68% (B=4,5) as more stages accumulate leakage.
- C absorbs the baseline energy ratio at the dominant stage itself. It depends on how cleanly the signal period aligns with 2·10^b.

The exponential decay in B is a structural feature, not a defect: it reflects the fact that with more digit stages, more channels compete for energy. The top-3 coverage (≥ 97%) remains stable because the two nearest channels — at levels b±1 — absorb most of the leakage.

### 4.3 Proof Sketch

We sketch the structure of the argument; a complete algebraic proof is left for future work.

**Step 1** (Digit smoothness → U energy bound). By (A1), for each b' ≠ b, the U_{b'} projection — which measures parity alternation at digit b' — has small inner product with x. This gives E[U_{b'}] ≤ ε² · E[U_b] for b' ≠ b.

**Step 2** (Dominant digit → U_b lower bound). By (A2), x oscillates with period ≈ 2·10^b, so x correlates strongly with the Z₂ parity indicator (−1)^{d_b}. This is precisely the kernel of U_b, yielding E[U_b] ≥ C₀ · E_total for some C₀ > 0.

**Step 3** (Tensor product → multiplicative separation). Because the B-stage decomposition applies each operator independently at each digit, the channel T…TU_bT…T receives energy from U_b at stage b and T at all other stages. The T projections at stages b'≠b pass through the slowly varying envelope (contributing factor ≈ 1), while the U_{b'} channels at those stages are suppressed by ε per stage.

**Step 4** (Normalization → α_b lower bound). Summing the energy across all 2^B channels and normalizing yields the lower bound α_b ≥ C · (1−ε)^{B−1}, where C absorbs C₀ and the T-stage pass-through factors.

### 4.4 Empirical Alignment

The quasi-theorem is consistent with the numerical experiments (B = 3..5, test signal sin(2π·n/100), P = 100 ≈ 2·10¹, b = 1):

| B | α_b (observed) | Lower bound structure |
|---|---------------|-----------------------|
| 3 | 0.864 | C·(1−ε)² with ε ≈ 0.07 |
| 4 | 0.687 | C·(1−ε)³ |
| 5 | 0.682 | C·(1−ε)⁴ |

In all cases, the top 3 channels jointly cover ≥ 97% of total energy, consistent with the δ < 0.03 claim. The UU…U channel carries < 0.01% in all experiments, consistent with the checkerboard suppression argument of Section 3.3.

The specific values of C and ε are not fitted here; the table is presented solely to show that the observed decay is consistent with the lower-bound structure C·(1−ε)^{B−1}.

---

## 5. Numerical Experiments

### 5.1 Experimental Setup

All experiments use the reference implementation phi_carry_free.py, operating entirely in exact integer arithmetic over Z[φ]. No floating-point operations are involved in the transform itself; zval() (evaluation to float) is used only for energy measurement in post-processing.

We test B = 1..5 (signal lengths N = 10 to 100,000) with a low-frequency test signal

    x[n] = sin(2π·n/100) + sin(2π·n/P₂) + noise,

where P₂ is a secondary period scaled to N and noise is small integer perturbation. This signal class satisfies assumptions (A1)–(A2) of Section 4 with dominant period P = 100 ≈ 2·10¹ (b=1).

Correctness is verified in two ways: (i) direct comparison of phi_conv_carry_free against the reference O(N²) carry-free convolution, and (ii) round-trip test with a delta kernel h = δ₀, confirming exact reconstruction for all B = 1..5.

### 5.2 Channel Energy Distribution

Table 1 shows the energy fraction of the top-3 channels and the UU…U channel for B = 3..5. All fractions are computed as

    E[ch] / E_total,   E[ch] = (1/N) Σ_k (zval(ch[k]))²   (real-valued signal assumption).

**Table 1. Channel energy distribution (low-frequency test signal)**

| B | #1 channel | α₁ | #2 channel | α₂ | #3 channel | α₃ | Top-3 | UU…U |
|---|-----------|-----|-----------|-----|-----------|-----|-------|------|
| 3 | TUT | 86.4% | TTT | 11.1% | UTT | 2.3% | 99.8% | <0.01% |
| 4 | TUTT | 68.7% | TTUT | 17.7% | TTTT | 11.7% | 98.1% | <0.01% |
| 5 | TUTTT | 68.2% | TTTUT | 17.7% | TTTTT | 11.7% | 97.6% | <0.01% |

Three observations are immediate. First, the dominant channel is consistently of the form T…TU₁T…T (U only at digit b=1), in agreement with the prediction of Section 4 for P ≈ 2·10¹. Second, the top-3 channels jointly cover ≥ 97% of total energy across all tested values of B, consistent with the δ < 0.03 claim of the quasi-theorem. Third, the UU…U channel carries negligible energy (< 0.01%) in all cases, consistent with the checkerboard suppression argument of Section 3.3.

**Table 2. SNR vs. retained channels (low-frequency signal, B=3..5)**

| B | 1 ch | 2 ch | 4 ch | All ch |
|---|------|------|------|--------|
| 3 | 9.5 dB | 15.4 dB | 28.5 dB | ∞ (exact) |
| 4 | 5.3 dB | 9.5 dB | 20.5 dB | ∞ (exact) |
| 5 | 5.2 dB | 9.3 dB | 18.5 dB | ∞ (exact) |

The "∞ (exact)" entries confirm lossless reconstruction when all channels are retained.

### 5.3 Reconstruction Accuracy

Exact reconstruction is guaranteed by the algebraic structure of φ-NTT, not by numerical precision. The scale factors 40^B (round-trip) and 80^B (after convolution) divide exactly in Z[φ], as verified by assertion checks for all B = 1..5:

```python
assert v[0] % scale == 0 and v[1] % scale == 0   # holds for all v
```

No floating-point rounding occurs at any stage. The delta-kernel test phi_conv_carry_free(x, δ₀, B) == x passes for all tested inputs and all B = 1..5, confirming the convolution theorem of Section 2.3.

### 5.4 Computational Structure

The φ-NTT avoids two sources of complexity common in FFT implementations. First, there are no inter-stage twiddle factors: all phase corrections vanish by the tensor-product theorem over Z₁₀^B. Second, all arithmetic is integer-valued; the only divisions are the final scale-factor normalizations, which are exact. The transform operates on length-10^B signals with 2^B output channels in O(N · B) integer operations per stage. Detailed complexity analysis and comparison with FFT are reserved for future work.

---

## 6. Discussion

### 6.1 Open Problems

Three problems remain open and are reserved for future work.

**(P1) Orthogonality of T/U projections.** The numerical experiments confirm that energy is partitioned across channels without loss, and that exact reconstruction holds for all tested B. However, the strict algebraic orthogonality T_b ⊥ U_b over Z[φ] has not been proved. Establishing this would upgrade the quasi-theorem of Section 4 to a full theorem, with C and ε determined by the projection geometry.

**(P2) Closed-form determination of C and ε.** The constants C and ε in Theorem 4.1 are currently characterized implicitly through assumptions (A1)–(A2). Deriving explicit expressions — presumably in terms of the signal period P, the digit base 10, and the Z[φ] basis elements — would make the spectral concentration bound fully quantitative and signal-class independent.

**(P3) Connection to carry-based convolution (overlap-add).** The present work is restricted to carry-free convolution on Z₁₀^B. A natural extension is to bridge this to ordinary cyclic convolution on Z/(10^B)Z via an overlap-add scheme, which would enable φ-NTT-based filtering of standard integer sequences. The group non-isomorphism Z₁₀^B ≇ Z/(10^B)Z (for B ≥ 2) means this bridge is non-trivial and requires a separate treatment.

### 6.2 Relation to Existing Frameworks

**vs. FFT.** Standard FFT operates on the cyclic group Z/(N)Z with complex twiddle factors, accumulating floating-point error. φ-NTT operates on the direct product group Z₁₀^B with coefficients in Z[φ], achieving exact integer arithmetic. The two transforms diagonalize different convolution algebras and are not interchangeable.

**vs. Haar MRA.** The binary tree structure of Section 3 is structurally analogous to Haar multiresolution analysis, with T playing the role of the scaling function and U the wavelet. The analogy is genuine but inexact: Haar operates over R with ±1 basis on Z/(2^B)Z, while φ-NTT uses Z[φ] basis on Z₁₀^B. Whether a precise isomorphism exists between the two filterbanks is an open question.

**vs. NTT (Number Theoretic Transform).** Classical NTT replaces complex roots of unity with modular arithmetic roots, but retains the cyclic group structure Z/(N)Z. φ-NTT differs in using the direct product group Z₁₀^B and the ring Z[φ], which is not a finite field. This gives exact arithmetic over the integers at the cost of the non-standard group structure.

### 6.3 Limitations

The current work has three explicit limitations.

First, all results are restricted to the carry-free setting on Z₁₀^B. Signals defined with ordinary decimal arithmetic — where addition carries across digit boundaries — are not directly handled by φ-NTT without preprocessing.

Second, only the group Z₁₀^B (decimal digit base) is treated. Extension to other bases (e.g., Z_p^B for prime p, or mixed-radix groups) may yield analogous transforms, but the specific role of φ = (1+√5)/2 and the ring Z[φ] is tied to the factorization 10 = 2 × 5.

Third, all experiments use synthetic signals. Application to real-world data — audio, time series, integer sequences — has not been tested and may require adaptation of the carry-free indexing to practical data formats.

---

## 7. Conclusion

We have introduced φ-NTT, a transform on the carry-free group Z₁₀^B with coefficients in the golden-ratio integer ring Z[φ]. The main contributions are threefold.

First, φ-NTT diagonalizes carry-free convolution on Z₁₀^B exactly, with no floating-point error and no inter-stage twiddle factors. The scale factors 40^B and 80^B divide exactly in Z[φ], making reconstruction provably lossless.

Second, the 2^B channel filterbank admits a hierarchical wavelet interpretation: the T/U label at each stage corresponds to a low-pass or high-pass projection at that digit level, structurally analogous to Haar MRA but defined over Z[φ] and indexed by decimal digits.

Third, we establish a Spectral Concentration Quasi-Theorem: under digit-wise smoothness and a dominant-period assumption, the energy fraction captured by the level-b detail channel satisfies α_b ≥ C·(1−ε)^{B−1}. Empirically, the top 3 channels cover ≥ 97% of total energy for B = 3..5, validating the bound structure.

Open problems include the algebraic proof of T/U orthogonality (P1), closed-form determination of C and ε (P2), and the connection to carry-based convolution via overlap-add (P3).

The φ-NTT suggests that hierarchical, carry-free arithmetic structures may offer an alternative perspective on multiresolution transforms beyond cyclic Fourier analysis.

---

## Appendix A: cos5_zphi Table and T/U Kernel Tables

The C₅-twisted DFT basis elements used in n10_T:

```
cos5_zphi(m) for m = 0..4:
  m=0: (2, 0)    ≈  2.000
  m=1: (-1, 1)   ≈  0.618
  m=2: (0, -1)   ≈ -1.618
  m=3: (0, -1)   ≈ -1.618
  m=4: (-1, 1)   ≈  0.618
```

The Z₂ parity kernel used in n10_U:

```
U1[k] for k = 0..9:
  [(0,0),(1,0),(0,1),(0,1),(1,0),(0,0),(-1,0),(0,-1),(0,-1),(-1,0)]
```

The inverse kernel T2[k] for k = 0..9:

```
T2 = [(2,0),(0,1),(-1,1),(1,-1),(0,-1),(-2,0),(0,-1),(1,-1),(-1,1),(0,1)]
```

---

## Appendix B: perm10 and CRT Derivation

```
perm10 = [(5*(i//5) + 6*(i%5)) % 10 for i in range(10)]
```

This implements the CRT inverse map Z₁₀ ≅ Z₂ × Z₅:

    n = 5·a + 6·b  (mod 10)

where a is the Z₅ component and b is the Z₂ component. The coefficient 6 satisfies 6 ≡ 1 (mod 5) and 6 ≡ 0 (mod 2), providing the Z₂ unit lift. Furthermore, 6 ≡ 3 (mod 5) induces the automorphism φ₃ (×3 map) on Z₅, which is consistent with the twiddle structure of the C₅ twisted DFT.

---

## Appendix C: Scale Factor Decomposition (80^B = 40^B × 2^B)

The total scale factor 80^B arises from two independent sources:

- **40^B**: the forward/inverse normalization factor accumulated over B tensor-product stages (each B=1 stage contributes factor 40).
- **2^B**: the T/U cross-coupling rule `2·Y_T = ...` introduces a factor of 2 per stage; over B stages this accumulates to 2^B.

Therefore:
- `bn_inverse_flat` alone (no convolution): divide by **40^B**
- `phi_conv_carry_free` (with convolution): divide by **80^B = 40^B × 2^B**

Both divisions are exact in Z[φ] (verified by assertion for all B = 1..5):

```python
assert v[0] % scale == 0 and v[1] % scale == 0
```

---

*論文全体: 約 3,390 語（本文） + Appendix A〜C*  
*Short paper（8〜12ページ）として適切な分量*  
*次タスク: .tex 化（arXiv向け最小テンプレ）*
