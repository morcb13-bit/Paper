import json, math, random
exec(open('testA_edge_chain.py').read().split("# 最初の胞を一つ")[0])
PHI=(1+5**0.5)/2
f=lambda x:(x[0]+x[1]*PHI)/2
a=0; b=next(j for j in range(N) if ADJ[a][j]); ring=cells_on_edge(a,b); c=ring[0]; d=next(w for w in ring if ADJ[c][w])
s=(a,b,c,d); cells=[list(s)]
for _ in range(15):
    s=label(step(s,1),2); cells.append(list(s))
assert set(cells[15])==set(cells[0])
used=sorted({v for cc in cells[:15] for v in cc})
coords={i:[round(f(x),9) for x in V[i]] for i in range(N)}
json.dump({'V':[coords[i] for i in range(N)],'cells':cells[:15],'pivot':[cc[:2] for cc in cells[:15]]},open('edge_ring.json','w'))
print(len(used), 'points; cells', len(cells)-1)
