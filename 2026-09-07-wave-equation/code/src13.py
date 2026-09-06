import sympy
from sympy import symbols, Poly, GF
from cage_min import order_of_x
x = symbols('x')
INT = ["x - 12", "x + 12", "x**2 + 144",
 "x**4 + 8*x**3 + 32*x**2 + 1152*x + 20736",
 "x**4 + 16*x**3 + 224*x**2 + 2304*x + 20736",
 "x**12 - 8*x**11 + 192*x**10 - 1408*x**9 + 35584*x**8 - 393216*x**7 + 5308416*x**6 - 56623104*x**5 + 737869824*x**4 - 4204265472*x**3 + 82556485632*x**2 - 495338913792*x + 8916100448256",
 "x**12 + 8*x**11 + 384*x**10 + 2944*x**9 + 84736*x**8 + 466944*x**7 + 13860864*x**6 + 67239936*x**5 + 1757085696*x**4 + 8790736896*x**3 + 165112971264*x**2 + 495338913792*x + 8916100448256"]
for m in [5,7,11,13,17,19,23]:
    o13 = sympy.n_order(m, 13) if m != 13 else 0
    hits = []
    for s in INT:
        P = Poly(sympy.sympify(s), x)
        for f, e in sympy.factor_list(Poly(P, x, domain=GF(m)))[1]:
            cl = [int(c) % m for c in reversed(f.all_coeffs())]
            o = order_of_x(cl, m)
            if o % 13 == 0:
                hits.append((P.degree(), f.degree(), o))
    print(f"m={m:>2}  ord_13(m)={o13 if o13 else "-":>2}  13を生む因数: "
          + ("なし" if not hits else
             ", ".join(f"整数側{a}次→mod{m}で{b}次(位数{c})" for a,b,c in hits)))
