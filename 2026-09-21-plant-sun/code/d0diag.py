exec(open('d0rw.py').read().split("for name,wall in")[0])
import collections
open_=Wset-set(WALL_L); d=bfs(FOOD,open_); S=walk(WALL_L,20000,1)
by=collections.defaultdict(list)
for i in open_: by[d[i]].append(S[i])
print('歩数ごとのS（最小〜最大）',[(k,min(v),max(v)) for k,v in sorted(by.items())])
good=0;tot=0;far=[]
for i in open_:
  if i==FOOD: continue
  p,ok=climb(S,open_,i); tot+=1
  if ok and len(p)-1==d[i]: good+=1
  else: far.append(d[i])
print('最短歩数で着く出発点',good,'/',tot,'着かない出発点の歩数',sorted(far))
