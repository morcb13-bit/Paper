import math, itertools, os
from collections import Counter, deque
exec(open('wind_core.py').read())
exec(open('pentagram_neuron.py').read().split('def run_tests')[0])
F, faces, SC = carrier()
n=len(SC); P4=PHI**4
L=[(i,j) for i,j in itertools.combinations(range(n),2) if abs(math.dist(SC[i],SC[j])-P4)<1e-6]
adj={i:[] for i in range(n)}
for i,j in L: adj[i].append(j); adj[j].append(i)
xs=[x for x,_ in SC]; cx=(min(xs)+max(xs))/2
axis=[i for i in range(n) if abs(SC[i][0]-cx)<1e-6]
# 二色
col={}; bip=True; comps=0
for s in range(n):
    if s in col: continue
    comps+=1; col[s]=0; dq=deque([s])
    while dq:
        a=dq.popleft()
        for b in adj[a]:
            if b not in col: col[b]=1-col[a]; dq.append(b)
            elif col[b]==col[a]: bip=False
c=Counter(col.values())
print(f"行 {os.environ.get('B13_ROWS','13')}  星 {n}  φ⁴の結合 {len(L)}本  次数 {dict(sorted(Counter(len(v) for v in adj.values()).items()))}"
      f"  連結成分 {comps}  二部 {bip} {c[0]}/{c[1]}  軸の上の星 {len(axis)}  軸x={cx:.6f}")
