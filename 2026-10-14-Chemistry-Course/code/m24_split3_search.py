import itertools
def run(m, K, label):
    E=6//m  # sub-edge (m=2 ->3, m=3 ->2); big cube edge 6
    cells=list(itertools.product(range(m),repeat=3))
    bigA={(6,0,0),(0,6,0),(0,0,6),(6,6,6)}
    def tet(c,useOdd):
        return frozenset(tuple(c[i]*E+E*n[i] for i in range(3)) for n in itertools.product([0,1],repeat=3) if (sum(n)%2==1)==useOdd)
    rules={"alt0":lambda c:sum(c)%2==1,"alt1":lambda c:sum(c)%2==0,"allA":lambda c:True,"allB":lambda c:False}
    face=lambda a,b:sum(abs(a[i]-b[i]) for i in range(3))==1
    # enumerate K-subsets with no face-adjacent pair (backtracking)
    sets=[]
    def bt(start,cur):
        if len(cur)==K: sets.append(tuple(cur)); return
        for i in range(start,len(cells)):
            c=cells[i]
            if any(face(c,d) for d in cur): continue
            cur.append(c); bt(i+1,cur); cur.pop()
    bt(0,[])
    res={}
    for rn,rule in rules.items():
        ok=[]
        for S in sets:
            T=[tet(c,rule(c)) for c in S]
            if any(len(a&b)>=2 for a,b in itertools.combinations(T,2)): continue
            # connectivity via shared vertices
            seen={0}; st=[0]
            while st:
                i=st.pop()
                for j in range(K):
                    if j not in seen and T[i]&T[j]: seen.add(j); st.append(j)
            if len(seen)!=K: continue
            if not bigA <= frozenset().union(*T): continue
            ok.append(S)
        res[rn]=ok
    print(label,"m=",m,"keep",K,"no-face sets",len(sets),{k:len(v) for k,v in res.items()})
    return res,tet,rules
# positive control
r2,_,_=run(2,4,"POS")
print(" POS sets:",{k:v for k,v in r2.items() if v})
# main
r3,tet,rules=run(3,9,"MAIN")
for k,v in r3.items():
    for S in v[:3]: print(" ",k,S)
# negative control: corners + center
NEG=tuple(c for c in itertools.product(range(3),repeat=3) if all(x in(0,2) for x in c) or c==(1,1,1))
for k,v in r3.items(): print("NEG in",k,":",NEG in v)

# post-hoc: parity used, symmetry under rotations preserving big A (cell coords 0..2)
ok=r3["allA"]
print("all even-sum cells:",all(sum(c)%2==0 for S in ok for c in S))
def mats():
    out=[]
    for p in itertools.permutations(range(3)):
        for s in itertools.product([1,-1],repeat=3):
            out.append((p,s))
    return out
def ap(g,c):
    p,s=g; v=[c[p[i]]-1 for i in range(3)]; return tuple(s[i]*v[i]+1 for i in range(3))
A={(2,0,0),(0,2,0),(0,0,2),(2,2,2)}
G=[g for g in mats() if {ap(g,c) for c in A}==A]
print("group size",len(G))
sym=[S for S in ok if all({ap(g,c) for c in S}==set(S) for g in G)]
print("symmetric sets",len(sym))
ev=[c for c in itertools.product(range(3),repeat=3) if sum(c)%2==0]
orb=set(frozenset(ap(g,c) for g in G) for c in ev); print("orbit sizes",sorted(len(o) for o in orb))
