# 検定A3（走らせる前に決めた）  体の量を M 環に限る。誕生のたびに M を超えるなら、いちばん古い環（しっぽ）を消す
#  OK  餌に届く ∧ 先端の路の歩数＝最短歩数（左の隙間15、右16）
#  NG  届かない／最短より長い
#  対照  匂いなし（+1だけ）→ 向きがないので届かない、または最短より長いはず ／ 隙間を右へ → 右の隙間を通る
exec(open('smell2.py').read())
cases=[('本番 左',WALL_L,True,M) for M in (6,10,16)]+[('対照 匂いなし 左',WALL_L,False,10),('対照 右の隙間',WALL_R,True,10)]
for name,wall,sense,M in cases:
  open_=Wset-set(wall); d=bfs(FOOD,open_); best=min(d[s] for s in SEED)
  for sd in (1,2):
    r=run(wall,L=1,T=12000,seed=sd,sense=sense,M=M); p=r['path']
    steps=None if p is None else len(p)-1; gap=[] if p is None else [x for x in p if x in (28,29,34,35)]
    near=sum(1 for i in W if r['E'][i] and d.get(i,99)<=2)
    print(f'{name} M={M} 乱数{sd}  届いたコマ{r["reach"]} 先端の歩数{steps}（最短{best}） 隙間{gap} 最後の生{r["hist"][-1]}（餌から2歩以内{near}） 生{r["births"]}消{r["deaths"]}',flush=True)
