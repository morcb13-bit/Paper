"""整数のまま（剰余を取らずに）一歩の最小多項式を出し、Z 上で因数分解する。"""
from fractions import Fraction as F
import sympy
from sympy import symbols, Poly
from cage_min import NB, DARTS, DI, E, C
x = symbols('x')

def step_int(psi):
    S = {v: sum(psi[DI[(u, v)]] for u in NB[v]) for v in NB}
    out = [0]*E
    for i, (u, v) in enumerate(DARTS):
        out[DI[(v, u)]] = C[len(NB[v])]*S[v] - 12*psi[i]
    return out

psi0 = [0]*E; psi0[0] = 1
rows, piv, combos = [], [], []
psi = [F(a) for a in psi0]
k = 0
while True:
    v = psi[:]; comb = [F(0)]*(k+1); comb[k] = F(1)
    for r, p, cb in zip(rows, piv, combos):
        f = v[p]
        if f:
            v = [a - f*b for a, b in zip(v, r)]
            comb = [a - f*b for a, b in zip(comb, cb + [F(0)]*(k+1-len(cb)))]
    nz = next((j for j, a in enumerate(v) if a), None)
    if nz is None:
        g = comb; break
    inv = 1/v[nz]
    rows.append([a*inv for a in v]); piv.append(nz)
    combos.append([a*inv for a in comb])
    psi = [F(a) for a in step_int([int(a) for a in psi])]; k += 1

P = Poly(list(reversed(g)), x).monic()
den = sympy.lcm([sympy.denom(c) for c in P.all_coeffs()])
P = Poly([sympy.Integer(c*den) for c in P.all_coeffs()], x)
print("整数上の最小多項式 deg =", P.degree())
for f, e in sympy.factor_list(P)[1]:
    print(f"  次数{f.degree():>3} 重複{e}   {f.as_expr()}")
