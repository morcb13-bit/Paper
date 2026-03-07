"""
Local cone invariant chain for cone-FFT v27.
Proposition B: V=(3,4) rotations -> C_CONE = 2/phi.
"""
import math
PHI = (1 + 5**0.5) / 2

def dot2(a, b):   return a[0]*b[0] + a[1]*b[1]
def norm2_2(a):   return dot2(a, a)
def sub2(a, b):   return (a[0]-b[0], a[1]-b[1])
def corr2(a, b):
    den = math.sqrt(norm2_2(a) * norm2_2(b))
    return dot2(a, b) / den if den else 0.0

def gram2(a, b):
    return ((dot2(a,a), dot2(a,b)), (dot2(a,b), dot2(b,b)))

def eigvals_sym2(G):
    a, b = G[0]; _, d = G[1]
    tr = a+d; det = a*d - b*b
    s  = math.sqrt(tr*tr - 4*det)
    return ((tr+s)/2, (tr-s)/2)

def rho_from_eigs(lm, ln):
    return float("inf") if ln==0 else math.sqrt(lm/ln)

def rho_a_from_eigs(lm, ln):
    if ln == 0: return float("inf")
    r = rho_from_eigs(lm, ln)
    return r - 1.0/r

def rho_to_R(rho): return (1-rho)/(1+rho)

def cone_chain(V, W):
    G = gram2(V[0], W[0])
    l1, l2 = eigvals_sym2(G)
    lm, ln = max(l1,l2), min(l1,l2)
    a, b = G[0]; _, d = G[1]
    tr = a+d; det = a*d - b*b
    rho_eig = rho_from_eigs(lm, ln)
    rho_a   = rho_a_from_eigs(lm, ln)
    R       = rho_to_R(rho_a)
    return {"G": G, "disc": tr*tr-4*det, "lam": (lm, ln),
            "rho_eig": rho_eig, "rho_a": rho_a,
            "R": R, "C_CONE": 2+2*R}
