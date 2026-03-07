"""
demo_verify.py — cone-FFT 全命題検証 (21項目)

命題A: ‖T'‖²+‖U'‖² = 40(‖T‖²+‖U'‖²)
命題B: Gram disc=3125 → rho_eig=φ² → rho_a=√5 → R=-φ⁻² → C_CONE=2/φ
命題C: 整数可逆性 (6T'+2U')/40=T, (2T'-6U')/40=U
命題D: 40の4重必然（代数/行列/数論/共有因子5）
"""
import math
import numpy as np
import sys
sys.path.insert(0, '/home/claude')

from parseval_core import (
    step_rule_40, step_rule_40_inv,
    build_fractal_integer_signal,
    build_fractal_integer_signal_inv,
)
from cone_local import cone_chain, gram2, eigvals_sym2
from cone_fft_invertible_v29 import cone_fft_roundtrip, BASE

PHI  = (1 + 5**0.5) / 2
PASS = "✓"
FAIL = "✗"

results = []

def check(name, ok, detail=""):
    mark = PASS if ok else FAIL
    results.append((name, ok))
    print(f"  {mark} {name}" + (f"  [{detail}]" if detail else ""))
    return ok

def section(title):
    print(f"\n{'─'*56}")
    print(f"  {title}")
    print(f"{'─'*56}")

# ── 命題A: 40倍則 ────────────────────────────────────────
section("命題A: ‖T'‖²+‖U'‖² = 40(‖T‖²+‖U‖²)")

rng = np.random.default_rng(0)
for trial, (n, label) in enumerate([(10,"小"), (100,"中"), (BASE//2,"大")]):
    T = list(rng.integers(-1000, 1001, n))
    U = list(rng.integers(-1000, 1001, n))
    T1, U1 = step_rule_40(T, U)
    n0 = sum(x*x for x in T)  + sum(x*x for x in U)
    n1 = sum(x*x for x in T1) + sum(x*x for x in U1)
    ratio = n1 / n0 if n0 else 0
    check(f"A-{trial+1}: n={n} ({label})", abs(ratio-40)<1e-9,
          f"比={ratio:.6f}")

# 1段後の各成分の型確認（整数のまま）
T = [1, 2, 3]; U = [4, 5, 6]
T1, U1 = step_rule_40(T, U)
check("A-4: 出力が整数リスト", all(isinstance(x, int) for x in T1+U1))

# ── 命題B: Gram → C_CONE 変換鎖 ─────────────────────────
section("命題B: disc=3125 → rho_eig=φ² → rho_a=√5 → R=-φ⁻² → C_CONE=2/φ")

CENTER   = (7, 7)
ROTATIONS = [(3,4),(4,-3),(-3,-4),(-4,3)]
V = ROTATIONS
W = [(V[(i+1)%4][0]-V[i][0], V[(i+1)%4][1]-V[i][1]) for i in range(4)]
ch = cone_chain(V, W)

check("B-1: disc = 3125",      ch["disc"] == 3125,
      f"disc={ch['disc']}")
check("B-2: rho_eig = φ²",     abs(ch["rho_eig"] - PHI**2) < 1e-9,
      f"rho_eig={ch['rho_eig']:.8f}")
check("B-3: rho_a = √5",       abs(ch["rho_a"] - math.sqrt(5)) < 1e-9,
      f"rho_a={ch['rho_a']:.8f}")
check("B-4: R = -φ⁻²",         abs(ch["R"] - (-1/PHI**2)) < 1e-9,
      f"R={ch['R']:.8f}")
check("B-5: C_CONE = 2/φ",     abs(ch["C_CONE"] - 2/PHI) < 1e-9,
      f"C_CONE={ch['C_CONE']:.8f}")

# ── 命題C: 整数可逆性 ────────────────────────────────────
section("命題C: 整数可逆性")

rng2 = np.random.default_rng(1)
for trial, n in enumerate([10, 100, 1000]):
    T0 = list(rng2.integers(-500, 501, n))
    U0 = list(rng2.integers(-500, 501, n))
    T1, U1 = step_rule_40(T0, U0)
    Tr, Ur = step_rule_40_inv(T1, U1)
    check(f"C-{trial+1}: 1段逆変換 n={n}", Tr==T0 and Ur==U0)

# 多段
T0 = list(range(100))
U0 = list(range(100, 200))
for B in [1, 2, 3]:
    TB, UB = build_fractal_integer_signal(B, T0, U0)
    Tr, Ur = build_fractal_integer_signal_inv(B, TB, UB)
    check(f"C-{B+3}: {B}段往復", Tr==T0 and Ur==U0, f"B={B}")

# 40整除性
T0 = list(rng2.integers(-1000, 1001, 200))
U0 = list(rng2.integers(-1000, 1001, 200))
T1, U1 = step_rule_40(T0, U0)
ok_div = all((6*t+2*u)%40==0 and (2*t-6*u)%40==0 for t,u in zip(T1,U1))
check("C-7: 40整除性", ok_div)

# ── 命題D: 40の4重必然 ──────────────────────────────────
section("命題D: 40の4重必然")

# 代数的: M=[[6,2],[2,-6]], M²=40I
M = [[6,2],[2,-6]]
M2 = [[M[0][0]*M[0][0]+M[0][1]*M[1][0], M[0][0]*M[0][1]+M[0][1]*M[1][1]],
      [M[1][0]*M[0][0]+M[1][1]*M[1][0], M[1][0]*M[0][1]+M[1][1]*M[1][1]]]
check("D-1: M² = 40I (代数)",
      M2 == [[40,0],[0,40]], f"M²={M2}")

# 行列: det=40, tr=0
det_M = M[0][0]*M[1][1] - M[0][1]*M[1][0]
tr_M  = M[0][0] + M[1][1]
check("D-2: det(M)=−40, tr(M)=0 (行列)",
      det_M == -40 and tr_M == 0, f"det={det_M}, tr={tr_M}")

# 数論: 40 = lcm(8,5) = lcm(5,8)
import math as _math
check("D-3: 40 = lcm(8,5) (数論)",
      _math.lcm(8,5) == 40, f"lcm(8,5)={_math.lcm(8,5)}")

# 共有因子: gcd(6²+2², 2²+6²) = 40
g = _math.gcd(6**2+2**2, 2**2+6**2)
check("D-4: gcd(6²+2², 2²+6²) = 40 (共有因子)",
      g == 40, f"gcd={g}")

# ── 完全可逆パイプライン ─────────────────────────────────
section("完全可逆パイプライン (cos+sin, B=1〜3)")

rng3 = np.random.default_rng(42)
f = list(rng3.integers(-100, 101, BASE).astype(int))
for B in [1, 2, 3]:
    f_rec = cone_fft_roundtrip(f, B)
    err   = max(abs(a-b) for a,b in zip(f,f_rec))
    check(f"Pipeline B={B}: 完全復元", err==0, f"max_err={err}")

# ── 結果サマリー ─────────────────────────────────────────
print(f"\n{'='*56}")
total = len(results)
passed = sum(1 for _, ok in results if ok)
print(f"  結果: {passed}/{total} 通過")
for name, ok in results:
    if not ok:
        print(f"  {FAIL} FAILED: {name}")
print(f"  {'✓ 全テスト通過' if passed==total else '✗ 要確認'}")
print(f"{'='*56}")
