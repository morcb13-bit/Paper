# H-Carrier: Stability-Constrained Residual Reinjection

H-Carrier is a lightweight residual reinjection mechanism designed for adaptive signal processing systems such as hearing aids, spatial audio pipelines, and mobile speech transcription.

Instead of discarding residual signals as noise, H-Carrier conditionally reinjects structured residual components as adaptive bias signals — while enforcing strict stability constraints.

---

## Core Idea

For each processing unit:

- `H_s` — signed residual  
- `H_m` — residual magnitude  

Only meaningful residuals are reinjected using a gated mechanism:

g = clip((H_m - τ) / κ, 0, 1)

b = α · g · H_s

Three guardrails ensure stability:

1. **Threshold Gate** — blocks stationary noise  
2. **Leaky Integration** — prevents drift  
3. **Loop-Gain Constraint** — guarantees bounded stability  

Effective loop gain must satisfy:

G_eff = |α| · g_max · C < 1

---

## Why It Matters

In multi-speaker (cocktail party) environments:

- Speaker transitions are transient-driven  
- Residual channels often contain onset cues  
- Standard pipelines discard these cues  

H-Carrier preserves useful transient structure without causing feedback instability.

---

## Example Use Cases

### Hearing Aids

- Adaptive mask threshold biasing  
- Binaural differential residual enhancement  
- Faster response to speaker switching  

### Mobile / Edge ASR

- Improved Voice Activity Detection recovery  
- Stabilized streaming thresholds  
- Reduced mis-trigger events  

---

## Design Goals

- O(N) per frame
- No additional FFT passes required
- Scalar state variables only
- Hardware-friendly implementation
- Explicit small-gain stability guarantee

---

## Philosophy

Residuals are not always noise.
When handled safely, they can carry structural information.

H-Carrier formalizes how to reuse them without sacrificing stability.

---

## Status

Concept specification and architectural integration examples are provided.
Empirical validation and clinical evaluation are future work.

---

## License

Specify your license here.
