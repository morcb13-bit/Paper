_PT={}
def ptr_s(g):
    ptr,sh=state[g]
    m=_PT.get(g)
    if m is None:
        base=U.zsub(ring5(g)[0],g); m={U.zmul(base,U.zt(e)):e for e in range(10)}; _PT[g]=m
    return m[ptr],sh
_R5={}
_ring5=ring5
def ring5(g):
    r=_R5.get(g)
    if r is None: r=_ring5(g); _R5[g]=r
    return r
