import itertools
def tet(c,s,h):  # cube center c, sign s (+1 A,-1 B), half edge h
    return {tuple(c[i]+e[i]*h for i in range(3)) for e in itertools.product([1,-1],repeat=3) if e[0]*e[1]*e[2]==s}
# (1) Sierpinski level1 in big cube [0,12]^3, big A = product-positive corners about center (6,6,6)
big=tet((6,6,6),1,6); print("bigA",sorted(big))
mids={tuple((a[i]+b[i])//2 for i in range(3)) for a in big for b in big}
kids=[]
for v in big:
    kid={tuple((v[i]+w[i])//2 for i in range(3)) for w in big}
    kids.append(kid)
cubes=[]
for n in itertools.product([0,1],repeat=3):
    c=tuple(3+6*x for x in n); s=(-1)**sum(n)
    cubes.append((n,s,tet(c,s,3),tet(c,-s,3)))
for kid in kids:
    m=[(n,s,'own' if kid==t else 'mirror') for n,s,t,tm in cubes if kid in (t,tm)]
    print(sorted(kid),m)
print("pair shares",[len(a&b) for a,b in itertools.combinations(kids,2)])
# (2) 3x3x3, corners + center, alternating sign
keep=[n for n in itertools.product([0,1,2],repeat=3) if all(x in(0,2) for x in n) or n==(1,1,1)]
T={n:tet(tuple(3+6*x for x in n),(-1)**sum(n),3) for n in keep}
for a,b in itertools.combinations(keep,2):
    k=T[a]&T[b]
    if k: print(a,b,k)
