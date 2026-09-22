# 検定A2（登録済み）  匂い＝地面だけを漂う分子（層1枚）、体＝core10の規則、感覚器＝隣と自分の分子数の大小
#  OK  餌に届く ∧ 先端の路の歩数＝壁を避けた最短歩数 ∧ 伸びて縮む（生の最大＞最後、消去あり）
#  NG  壁の前で止まる／届かない／最短より長い
#  対照1 +1なし（匂いだけ） 対照2 匂いなし（+1だけ） 対照3 隙間を右へ
exec(open('smell.py').read())
import sys
cases=[('本番 左の隙間',WALL_L,1,True),('対照1 匂いだけ',WALL_L,0,True),('対照2 +1だけ',WALL_L,1,False),('対照3 右の隙間',WALL_R,1,True)]
for name,wall,drive,sense in cases:
  open_=Wset-set(wall); d=bfs(FOOD,open_); best=min(d[s] for s in SEED)
  for sd in (1,2):
    r=run(wall,L=1,T=12000,seed=sd,drive=drive,sense=sense); h=r['hist']; p=r['path']
    steps=None if p is None else len(p)-1
    gap=[] if p is None else [x for x in p if x in (28,29,34,35)]
    front=max((R[i][0] for i in W if r['E'][i]),default=-1)
    print(f'{name} 乱数{sd}  届いたコマ{r["reach"]} 先端の歩数{steps}（最短{best}） 隙間{gap} 生 最大{max(h)} 最後{h[-1]} 生{r["births"]}消{r["deaths"]} 最も深い段{front}',flush=True)
