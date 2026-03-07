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


def step_rule_40(T: list[int], U: list[int]) -> tuple[list[int], list[int]]:
    """1段整数再帰。ノルム係数 40。"""
    T_next = [6*t + 2*u for t, u in zip(T, U)]
    U_next = [2*t - 6*u for t, u in zip(T, U)]
    return T_next, U_next


def step_rule_40_inv(T: list[int], U: list[int]) -> tuple[list[int], list[int]]:
    """
    1段整数逆再帰。Proposition C より厳密に整数。
        T_prev = (6T + 2U) / 40
        U_prev = (2T - 6U) / 40
    """
    T_prev = [(6*t + 2*u) // 40 for t, u in zip(T, U)]
    U_prev = [(2*t - 6*u) // 40 for t, u in zip(T, U)]
    return T_prev, U_prev


def build_fractal_integer_signal(
    B: int,
    T0: list[int],
    U0: list[int],
) -> tuple[list[int], list[int]]:
    """B段の整数再帰。ノルム係数 40^B。"""
    T, U = T0[:], U0[:]
    for _ in range(B):
        T, U = step_rule_40(T, U)
    return T, U


def build_fractal_integer_signal_inv(
    B: int,
    T: list[int],
    U: list[int],
) -> tuple[list[int], list[int]]:
    """
    B段の整数逆再帰。build_fractal_integer_signal の完全逆写像。
    Proposition C より整数上で閉じる。
    """
    for _ in range(B):
        T, U = step_rule_40_inv(T, U)
    return T, U
  
