import itertools
from collections import Counter
def kids(o,e):
    h=e//2
    return [tuple(o[i]+h*n[i] for i in range(3)) for n in itertools.product([0,1],repeat=3) if sum(n)%2==1]
def tetA(o,e): return [tuple(o[i]+e*n[i] for i in range(3)) for n in itertools.product([0,1],repeat=3) if sum(n)%2==1]
cubes=[((0,0,0),96)]; rows={}
for k in range(1,5):
    cubes=[(c,e//2) for o,e in cubes for c in kids(o,e)]
    cnt=Counter(v for o,e in cubes for v in tetA(o,e))
    T=len(cubes); V=len(cnt); S=sum(1 for x in cnt.values() if x==2)
    rows[k]=dict(T=T,V=V,S=S,T2=T//2,V2=V//2,S2=S//2)
shell=[2*n*n for n in range(1,6)]  # 2,8,18,32,50 (足し算でも同じ)
print("shell",shell)
for key in ["T","V","S","T2","V2","S2"]:
    seq=[rows[k][key] for k in range(1,5)]
    hit=[x for x in seq if x in shell]
    print(key,seq,"hit",hit)
