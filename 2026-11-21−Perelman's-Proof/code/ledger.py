# firing ledger: Fibonacci by addition only, subtract m once when sum >= m
def run(m, a=1, b=1):
    seq=[(a,0),(b,0)]; r0,r1=a,b
    while True:
        s=r0+r1; f=0
        if s>=m: s-=m; f=1
        assert s<m
        seq.append((s,f)); r0,r1=r1,s
        if (r0,r1)==(a,b): break
    # period: steps until (a,b) repeats -> residues list excluding last two (they're the restart)
    res=seq[:-2]; return res
def stats(m):
    res=run(m); P=len(res)
    # firings belong to steps producing res[2..P+1]; count firings over one full period of produced terms
    full=run(m)
    # recompute firing per term over one period cyclically
    r=[x for x,_ in res]
    F=sum(1 for n in range(P) if r[n]+r[(n+1)%P]>=m)
    S=sum(r)
    return P,F,S,r
for m in [7,11,13,19,23,29,41,61]:
    P,F,S,r=stats(m); print(m,P,F,S,m*F==S)
ok=all((lambda P,F,S,r:m*F==S)(*stats(m)) for m in range(2,400)); print("mF=S 2..399",ok)
P,F,S,r=stats(13); print(r)
print([1 if r[n]+r[(n+1)%P]>=13 else 0 for n in range(P)])
print("half-period sign", all((r[n+14]+r[n])%13==0 for n in range(14)))
print("distinct F over 30 moduli 2..31:", len({stats(m)[1] for m in range(2,32)}))
