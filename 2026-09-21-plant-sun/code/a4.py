# A3 と同じ判定。餌が毎コマ出す分子を Q 個に増やす
exec(open('smell3.py').read())
for Qv in (10,100):
  Q=Qv
  for wall,gname in ((WALL_L,'左'),(WALL_R,'右')):
    open_=Wset-set(wall); d=bfs(FOOD,open_); best=min(d[s] for s in SEED)
    for M in (10,16):
      for sd in (1,2):
        r=run(wall,L=1,T=6000,seed=sd,M=M); p=r['path']
        steps=None if p is None else len(p)-1; gap=[] if p is None else sorted(set(x for x in p if x in (28,29,34,35)))
        near=sum(1 for i in W if r['E'][i] and d.get(i,99)<=2)
        print(f'Q={Q} 隙間{gname} M={M} 乱数{sd}  届いたコマ{r["reach"]} 先端の歩数{steps}（最短{best}） 隙間{gap} 最後の生{r["hist"][-1]}（餌から2歩以内{near}）',flush=True)
