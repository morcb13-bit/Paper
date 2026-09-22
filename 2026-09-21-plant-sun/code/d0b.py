exec(open('maze.py').read().split("for name,wall in")[0])
import collections
for T in (10000,40000):
  open_=Wset-set(WALL_L); d=bfs(FOOD,open_); C={i:0 for i in W}; first={}
  for t in range(T):
    mol_step(C,open_,FOOD)
    for i in open_:
      if C[i]>0 and i not in first: first[i]=t
  ks=[i for i in open_ if i in first]
  inv=sum(1 for a in ks for b in ks if d[a]<d[b] and first[a]>first[b])
  byd=collections.defaultdict(list)
  for i in ks: byd[d[i]].append(first[i])
  print(T,'届いた',len(ks),'/',len(open_),'逆転',inv,'歩数ごとの最初と最後',[(k,min(v),max(v)) for k,v in sorted(byd.items())][:22])
