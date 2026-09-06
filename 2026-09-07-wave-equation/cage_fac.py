import sys, sympy
from sympy import Poly, symbols, GF, factor_list, gcd
from cage_min import minpoly, order_of_x, E, step
x = symbols('x')

def analyse(m):
    psi0 = [0]*E; psi0[0] = 1
    g = minpoly(psi0, m)
    P = Poly(list(reversed(g)), x, domain=GF(m))
    fl = sympy.factor_list(P)
    facs = []
    for f, e in fl[1]:
        cl = [int(c) % m for c in reversed(f.all_coeffs())]
        d = len(cl)-1
        o = order_of_x(cl, m) if d >= 1 else None
        facs.append((d, e, o, f.as_expr()))
    facs.sort(key=lambda t: (t[0], t[1], str(t[3])))
    N = 1
    for d, e, o, _ in facs:
        if o: N = sympy.ilcm(N, o)
    emax = max(e for _, e, _, _ in facs)
    k = 1
    while m**k < emax: k += 1
    N *= m**(0 if emax == 1 else k)
    return g, facs, N

for m in [int(a) for a in sys.argv[1:]]:
    g, facs, N = analyse(m)
    print(f"=== m={m}  deg g={len(g)-1} ===")
    for d, e, o, ex in facs:
        mark = "  ← 13" if o and o % 13 == 0 else ""
        print(f"  次数{d:>3} 重複{e}  x の位数 {o}  {sympy.factorint(o) if o else ''}{mark}")
        print(f"        {ex}")
    print(f"  合成した周期 = {N} = {sympy.factorint(N)}")
