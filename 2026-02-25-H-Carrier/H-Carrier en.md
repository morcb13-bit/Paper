# Specification Section: H-Carrier (Useful Residual Reinjection) Propagation and Stability Control

---

## Overview

This processor does not discard internally generated residuals (heat) as mere losses.  
Instead, residual signal \( H \) is reused as an **adaptive bias** applied to the next stage (or next time step), affecting thresholds, coupling coefficients, and decay rates.

However, unconditional reinjection of \( H \) may cause feedback oscillation, drift, or saturation.  
Therefore, the following three guardrails (stability conditions) are **mandatory I/O specifications**.

---

## 1. Residual Usefulness Guarantee (Threshold Gate)

### Purpose

Prevent circulation of small stationary noise (e.g., quantization error or thermal-equivalent noise) and allow only structurally meaningful residuals (such as edge transitions) to act as carriers.

---

### Decomposition

- \( H_s \): signed residual  
- \( H_m \): residual magnitude (energy)

Example definitions (per sample or per tile):

\[
H_s := \mathrm{mean}(\text{heat}) \quad \text{or} \quad \sum(\text{heat})
\]

\[
H_m := \mathrm{mean}(|\text{heat}|) \quad \text{or} \quad \sum(|\text{heat}|) \quad \text{or} \quad \sqrt{\sum(\text{heat}^2)}
\]

---

### Threshold Gate

Let the unit-specific noise floor be threshold \( \tau \).  
Reinjection is permitted only when:

\[
H_m > \tau
\]

Gate coefficient:

\[
g := \mathrm{clip}\left(\frac{H_m - \tau}{\kappa}, 0, 1\right)
\]

- \( \kappa \): gate transition width  
  - \( \kappa \to 0 \): hard gate  
  - \( \kappa > 0 \): soft gate  

---

### Reinjection Bias (Basic Form)

\[
b := \alpha \cdot g \cdot H_s
\]

---

### Recommended Notes

- \( \tau \) should be set to:  
  quantization-noise-equivalent level + safety margin under environmental noise.
- For L/R cross-injection, construct \( H \) as a differential residual to suppress common-mode components.

---

## 2. Drift Prevention via Leaky Integration

### Purpose

Prevent residual reinjection from accumulating into DC offset, saturation, bias drift, or slow adaptation under changing environments.

---

### Leaky Accumulation

\[
\hat{H}_s[t] := (1 - \lambda)\hat{H}_s[t-1] + \lambda H_s[t]
\]

where:

\[
0 < \lambda \le 1
\]

Similarly for magnitude if required:

\[
\hat{H}_m[t] := (1 - \lambda_m)\hat{H}_m[t-1] + \lambda_m H_m[t]
\]

---

### Reinjection with Leak

\[
b[t] := \alpha \cdot g[t] \cdot \hat{H}_s[t]
\]

---

### Specification Parameters

- \( \lambda \) (leak coefficient): adaptation time constant  
  - Small \( \lambda \): slow adaptation  
  - Large \( \lambda \): fast adaptation  

Design recommendation:  
Select \( \lambda \) based on the use case and verify absence of drift under rapidly changing noise conditions.

---

## 3. Stability Loop-Gain Constraint

### Purpose

Prevent H-Carrier reinjection from behaving as uncontrolled positive feedback and ensure practical stability under all input conditions.

---

### Fundamental Constraint (Sufficient Condition)

The effective reinjection loop gain must remain below unity:

\[
G_{eff} := |\alpha| \cdot g_{max} \cdot C
\]

Requirement:

\[
G_{eff} < 1
\]

Where:

- \( \alpha \): reinjection gain (design parameter)
- \( g_{max} \): maximum gate coefficient (typically 1)
- \( C \): upper bound of the unit coupling amplification factor  
  (depends on topology and update rule)

---

### Correlation with Inverse Fibonacci Convergence (Recommended)

If the system employs inverse-Fibonacci convergence with:

\[
\text{INV\_PHI} = \phi^{-1}
\]

a conservative design recommendation is:

\[
|\alpha| \le (1 - \text{INV\_PHI}) \cdot \beta
\]

where:

\[
0 < \beta < 1
\]

Typical range: \( \beta = 0.5 \sim 0.9 \)

---

### Intuition

- Larger INV_PHI (weaker convergence) → smaller \( \alpha \) required  
- Stronger convergence → greater H-Carrier tolerance  

---

### Implementation Notes

- Saturation (clipping) is mandatory:  
  \[
  |b| \le b_{max}
  \]
- For inter-unit cross-injection (L ↔ R), combine with:
  - band-limited application, or  
  - differential residual construction  

---

## I/O Specification (Minimum Set)

Each unit must expose:

- \( y_I, y_Q \): primary signal outputs (I/Q or equivalent two-phase components)
- \( H_s, H_m \): residual outputs (signed / magnitude)

If H-Carrier reinjection is enabled, the conditions in Sections 1–3 must be satisfied.

If not satisfied, residuals must be treated as **observation-only (diagnostic port)**.

---

## Terminology

- **H (heat)**: residual / dissipative component; not merely loss, but a carrier of un-settled structural information.
- **H-Carrier**: mechanism that transports \( H \) as adaptive bias to subsequent stages.
- **g (gate)**: coefficient allowing only useful residuals to pass.
- **λ (leak)**: forgetting coefficient preventing drift.
- **α (gain)**: reinjection strength determining loop gain.
