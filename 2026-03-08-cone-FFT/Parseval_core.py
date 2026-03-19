"""
Parseval core for cone-FFT v28.
Implements the integer 2-channel recursion
    (T, U) -> (6T+2U, 2T-6U),
with one-step norm scaling exactly equal to 40.
Proposition A:
    |T'|² + |U'|² = 40(|T|² + |U|²)
Proposition C (integer invertibility):
    The recursion is exactly invertible over the integers:
        T = (6T' + 2U') / 40
        U = (2T' - 6U') / 40
    Proof: 6T'+2U' = 6(6T+2U)+2(2T-6U) = 40T
           2T'-6U' = 2(6T+2U)-6(2T-6U) = 40U  □
"""

def step_rule_40(T, U):
    return [6*t + 2*u for t, u in zip(T, U)], \
           [2*t - 6*u for t, u in zip(T, U)]

def step_rule_40_inv(T, U):
    return [(6*t + 2*u) // 40 for t, u in zip(T, U)], \
           [(2*t - 6*u) // 40 for t, u in zip(T, U)]

def build_fractal_integer_signal(B, T0, U0):
    T, U = T0[:], U0[:]
    for _ in range(B): T, U = step_rule_40(T, U)
    return T, U

def build_fractal_integer_signal_inv(B, T, U):
    for _ in range(B): T, U = step_rule_40_inv(T, U)
    return T, U
