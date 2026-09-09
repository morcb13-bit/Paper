#  gen_data.py ── 動く図に埋め込む担体データを厳密整数から作る
import json, math
exec(open('rt_walk.py').read().split('# ================================================================ 実行')[0])
from collections import Counter

verts,n_ico,edges,elen2 = rt_graph()
faces,vadj = rt_faces(verts,n_ico,edges)
PHI = (1+5**0.5)/2
def tofloat(v): return [c[0]+c[1]*PHI for c in v]
XYZ = [tofloat(v) for v in verts]

fedges={}
for i,f in enumerate(faces):
    for k in range(4): fedges.setdefault(frozenset((f[k],f[(k+1)%4])),[]).append(i)
nbr={i:[] for i in range(30)}
for e,fs in fedges.items():
    if len(fs)==2: nbr[fs[0]].append(fs[1]); nbr[fs[1]].append(fs[0])
for i in nbr: nbr[i].sort()

def direction(u,v):
    d=vsub(verts[u],verts[v])
    for c in d:
        s=psign(c)
        if s:
            if s<0: d=tuple(pneg(x) for x in d)
            break
    return d
dirs={}
for a,b in edges: dirs.setdefault(direction(a,b),[]).append(frozenset((a,b)))
zones=[]
for dv,es in sorted(dirs.items()):
    fs=set()
    for e in es: fs.update(fedges[e])
    zones.append(sorted(fs))
belt=set(zones[0])
restf=[f for f in range(30) if f not in belt]
seen=set(); caps=[]
for f in restf:
    if f in seen: continue
    comp=[f]; seen.add(f); st=[f]
    while st:
        x=st.pop()
        for y in nbr[x]:
            if y in restf and y not in seen: seen.add(y); comp.append(y); st.append(y)
    caps.append(sorted(comp))
cap=caps[0]

def carrier(sub):
    S=set(sub)
    nb={f:[g for g in nbr[f] if g in S] for f in sorted(sub)}
    degs=sorted(set(len(nb[f]) for f in nb))
    L=1
    for d in degs: L=L*d//math.gcd(L,d)
    t=L if L%2 else L//2                      # 2t が全次数で割れる最小の t
    while any((2*t)%d for d in degs): t+=1
    arcs=[]
    for f in sorted(nb):
        for g in nb[f]:
            k=[kk for kk in range(4)
               if g in fedges[frozenset((faces[f][kk],faces[f][(kk+1)%4]))]][0]
            arcs.append([f,g,k])
    idx={(a[0],a[1]):i for i,a in enumerate(arcs)}
    rows=[]
    for (u,v,_k) in arcs:
        d=len(nb[u]); c=(2*t)//d
        row={}
        for w in nb[u]: row[idx[(w,u)]]=row.get(idx[(w,u)],0)+c
        row[idx[(v,u)]]=row.get(idx[(v,u)],0)-t
        rows.append(sorted([k,val] for k,val in row.items() if val))
    walls=[f for f in sorted(nb) if len(nb[f])==min(degs)] if len(degs)>1 else []
    # 既定の起点は、既定の視角で手前に来る面
    AX=[0.5257,0.0,-0.8507]
    def cendir(f):
        c=[sum(XYZ[i][k] for i in faces[f])/4 for k in range(3)]
        L=sum(x*x for x in c)**0.5
        return sum(c[k]*AX[k] for k in range(3))/L
    seed=max(sorted(sub), key=cendir)
    return {"faces":sorted(sub),"arcs":arcs,"rows":rows,"t":t,"seed":seed,
            "deg":{str(k):v for k,v in sorted(Counter(len(nb[f]) for f in nb).items())},
            "walls":walls,"edges":sum(len(nb[f]) for f in nb)//2}

data={"xyz":XYZ,"faces":faces,
      "carriers":{"cap":carrier(cap),"belt":carrier(sorted(belt)),"rt":carrier(list(range(30)))}}
json.dump(data,open("rt_data.json","w"),separators=(",",":"))
for k,v in data["carriers"].items():
    print(k, "面",len(v["faces"]), "次数",v["deg"], "辺",v["edges"], "t",v["t"], "有向辺",len(v["arcs"]))
