exec(open('pi1.py').read().split("def climb")[0])
def climb_all(D,open_,start,d):
  cur={start}; steps=0; used=set()
  while True:
    if cur=={FOOD}: return '着く',steps,used
    nxt=set()
    for c in cur:
      if c==FOOD: nxt.add(c); continue
      nb=[j for j in ADJ[c] if j in open_]; m=max(D[j] for j in nb)
      if m<=D[c]: return f'止まる（段{R[c][0]}）',steps,used
      nxt|={j for j in nb if D[j]==m}
    cur=nxt; steps+=1; used|=cur
    if steps>100: return '回る',steps,used
for name,wall,ignore in (('本番 左の隙間',WALL_L,False),('対照1 壁を無視して広げる',WALL_L,True),('対照2 右の隙間',WALL_R,False)):
  open_=Wset-set(wall); d=bfs(FOOD,open_)
  for T in (100,300):
    D=field(Wset if ignore else open_,T); out=[]
    for s in SEED:
      st,n,used=climb_all(D,open_,s,d)
      out.append(f'種{s}:{st} 歩数{n}/最短{d[s]} 隙間{sorted(x for x in used if x in (28,29,34,35))} 枝の環{len(used)} {"OK" if st=="着く" and n==d[s] else "NG"}')
    print(name,f'T={T}',' ｜ '.join(out))
