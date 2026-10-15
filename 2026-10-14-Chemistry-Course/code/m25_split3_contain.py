import itertools
from collections import Counter
E=2
cells=list(itertools.product(range(3),repeat=3))
A={(2,0,0),(0,2,0),(0,0,2),(2,2,2)}
bigA={(6,0,0),(0,6,0),(0,0,6),(6,6,6)}
def tet(c): return frozenset(tuple(c[i]*E+E*n[i] for i in range(3)) for n in itertools.product([0,1],repeat=3) if sum(n)%2==1)
face=lambda a,b:sum(abs(a[i]-b[i]) for i in range(3))==1
ev=[c for c in cells if sum(c)%2==0]          # 14 cells: 8 corners + 6 face centres
def ok(S):
    T=[tet(c) for c in S]
    if any(len(a&b)>=2 for a,b in itertools.combinations(T,2)):return False
    seen={0};st=[0]
    while st:
        i=st.pop()
        for j in range(len(T)):
            if j not in seen and T[i]&T[j]:seen.add(j);st.append(j)
    return len(seen)==len(T) and bigA<=frozenset().union(*T)
# all subsets of every size from the 14 even cells (all faces non-adjacent automatically)
bysize=Counter()
good9=[]
for k in range(1,15):
    for S in itertools.combinations(ev,k):
        if ok(S):
            bysize[k]+=1
            if k==9: good9.append(S)
print("connected, contains big A, by size:",dict(sorted(bysize.items())))
# symmetry classes of the 150
G=[]
for p in itertools.permutations(range(3)):
    for s in itertools.product([1,-1],repeat=3):
        g=(p,s)
        f=lambda c,g=g:tuple(g[1][i]*(c[g[0][i]]-1)+1 for i in range(3))
        if {f(c) for c in A}==A: G.append(f)
cls={}
for S in good9:
    key=min(tuple(sorted(f(c) for c in S)) for f in G)
    cls.setdefault(key,[]).append(S)
print("classes",len(cls))
cornerA=A; cornerB={(0,0,0),(2,2,0),(2,0,2),(0,2,2)}
for key,v in cls.items():
    S=key; T=[tet(c) for c in S]; cnt=Counter(x for t in T for x in t)
    print(" size",len(v),"A-corners",len(set(S)&cornerA),"B-corners",len(set(S)&cornerB),"faces",sum(1 for c in S if sum(1 for x in c if x==1)==2),
          "points",len(cnt),"shared",Counter(cnt.values()))

# full parity class (14)
S14=tuple(ev); T=[tet(c) for c in S14]; cnt=Counter(x for t in T for x in t)
print("14-set ok",ok(S14),"points",len(cnt),"mult",Counter(cnt.values()))
# inside big A? integer half-space test via determinants
V=[(6,0,0),(0,6,0),(0,0,6),(6,6,6)]
def det(a,b,c):return a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])
def sub(a,b):return tuple(a[i]-b[i] for i in range(3))
def inside(p):
    for i in range(4):
        o=[V[j] for j in range(4) if j!=i]
        s1=det(sub(o[1],o[0]),sub(o[2],o[0]),sub(V[i],o[0]))
        s2=det(sub(o[1],o[0]),sub(o[2],o[0]),sub(p,o[0]))
        if s1*s2<0:return False
    return True
for orient in (1,0):
    inn=[c for c in cells if all(inside(p) for p in frozenset(tuple(c[i]*E+E*n[i] for i in range(3)) for n in itertools.product([0,1],repeat=3) if sum(n)%2==orient))]
    print("orient",'A' if orient else 'B',"cells whose tetra lie inside big A:",inn)
# 2-split check of same containment

def tetO(c,o): return frozenset(tuple(c[i]*E+E*n[i] for i in range(3)) for n in itertools.product([0,1],repeat=3) if sum(n)%2==o)
up=[c for c in [(0,0,2),(0,2,0),(2,0,0),(2,2,2),(0,1,1),(1,0,1),(1,1,0),(1,1,2),(1,2,1),(2,1,1)]]
TU=[tetO(c,1) for c in up]; cA=tetO((1,1,1),1); cB=tetO((1,1,1),0)
print("centre A shares edge with an upright:",any(len(cA&t)>=2 for t in TU))
print("centre B shares edge:",any(len(cB&t)>=2 for t in TU),"touch points:",sum(1 for t in TU if cB&t))
print("uprights pairwise edge-share:",any(len(a&b)>=2 for a,b in itertools.combinations(TU,2)))
# volume x3 : small tetra 8, octahedron 32, big A 216
print("10*8 + 1*8 + 4*32 =",10*8+8+4*32)
# 2-split same containment
E2=3
def t2(c,o): return [tuple(c[i]*E2+E2*n[i] for i in range(3)) for n in itertools.product([0,1],repeat=3) if sum(n)%2==o]
for o in (1,0):
    print("2-split orient",'A' if o else 'B',[c for c in itertools.product(range(2),repeat=3) if all(inside(p) for p in t2(c,o))])
