import itertools
L={(2*(a+d),2*(b+d),2*(c+d)) for a in range(4) for b in range(4) for c in range(4) for d in range(4) if a+b+c+d==3}
print(len(L))
cent=[c for c in itertools.product(range(0,7),repeat=3) if all(tuple(c[i]+s*(i==j)*2 for i in range(3)) in L for j in range(3) for s in (1,-1))]
print("oct centres",cent)
def tet(c,o): return frozenset(tuple(c[i]*2+2*n[i] for i in range(3)) for n in itertools.product([0,1],repeat=3) if sum(n)%2==o)
up=[(0,0,2),(0,2,0),(2,0,0),(2,2,2),(0,1,1),(1,0,1),(1,1,0),(1,1,2),(1,2,1),(2,1,1)]
B=tet((1,1,1),0)
print("B verts",sorted(B))
for c in up:
    s=tet(c,1)&B
    if len(s)>=2: print(c,sorted(s))
print("all tet verts in L:",all(v in L for c in up for v in tet(c,1)), all(v in L for v in B))
