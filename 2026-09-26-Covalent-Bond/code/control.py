import random, sympy as sp, helix_check as h
from run_exact import step, measure
def chain_signs(signs):
    T=h.V[:]; cs=[h.cen(T)]
    for s in signs: T=step(T,s); cs.append(h.cen(T))
    return cs
random.seed(13)
C={"対照A：途中でねじの向きを一度だけ反転（1,1,1,1 の後に交互）":[1,1,1,1]+[-1,1]*4,
   "対照B：無作為の符号（seed 13）":[random.choice((1,-1)) for _ in range(12)]}
for k,s in C.items():
    ok,_=measure(chain_signs(s),k); print("符号",s,"→ 崩れたか:",not ok)
