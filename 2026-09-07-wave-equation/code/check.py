import sys, sympy
from cage_min import minpoly, polypow, E

def check(m, N):
    psi0 = [0]*E; psi0[0] = 1
    g = minpoly(psi0, m); d = len(g)-1
    one = [1]+[0]*(d-1); xx = [0,1]+[0]*(d-2)
    ok = polypow(xx, N, g, m) == one
    minimal = all(polypow(xx, N//p, g, m) != one for p in sympy.primefactors(N))
    print(f"m={m}  N={N}  x^N=1 : {ok}   どの素因数を落としても≠1 : {minimal}")

check(5, 1560); check(5, 780)
check(13, 5859189829080)
check(23, 39155492640)
