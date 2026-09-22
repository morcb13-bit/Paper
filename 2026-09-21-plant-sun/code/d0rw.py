# 検定D0'  餌の分子をランダムウォークにする（監督の指定：ブラウン運動、道を選ばない）
#  分子  餌の環が毎コマ1個出す。各分子は毎コマ、自分の環の4つの口から1つを等しく選ぶ。
#        口の先が通れる環なら移る／壁なら留まる／担体の外（隣が4に足りない分）なら出て行く
#  判定の替え（走らせる前）：乱歩では最初に届くコマの順は判定に使えない。アメーバが読むのは数なので、
#        各環に居た分子の数をコマごとに足した積算 S（整数）で判定する
#  OK   種（頂点の環）から、S が最も大きい通れる隣へ移ることを繰り返すと、最短歩数ちょうどで餌に着く
#  NG   途中で止まる（周りより大きい所）／最短歩数より長い
#  負の対照  S をでたらめに並べ替える → 着かないはず
import numpy as np
exec(open('maze.py').read().split("def mol_step")[0])
def walk(wall,T,seed):
  rng=np.random.default_rng(seed); open_=Wset-set(wall)
  idx={i:n for n,i in enumerate(W)}; C=np.zeros(len(W),dtype=np.int64); S=np.zeros(len(W),dtype=np.int64)
  slots={i:sorted(ADJ[i])+[None]*(4-len(ADJ[i])) for i in W}
  for t in range(T):
    C[idx[FOOD]]+=1; N=C.copy(); 
    for i in open_:
      c=C[idx[i]]
      if c==0: continue
      m=rng.multinomial(c,[0.25]*4); N[idx[i]]-=c
      for k,j in enumerate(slots[i]):
        if m[k]==0: continue
        if j is None: continue                      # 外へ
        if j in open_: N[idx[j]]+=m[k]
        else: N[idx[i]]+=m[k]                       # 壁で留まる
    C=N; S+=C
  return {i:int(S[idx[i]]) for i in W}
def climb(S,open_,start):
  cur=start;path=[cur]
  while cur!=FOOD:
    nb=[j for j in ADJ[cur] if j in open_]; b=max(nb,key=lambda j:S[j])
    if S[b]<=S[cur]: return path,False
    cur=b;path.append(cur)
    if len(path)>100: return path,False
  return path,True
WALL_L=list(range(30,36));WALL_R=list(range(28,34));FOOD=97
for name,wall in (('左の隙間',WALL_L),('右の隙間',WALL_R)):
  open_=Wset-set(wall); d=bfs(FOOD,open_)
  for sd in (1,2,3):
    S=walk(wall,20000,sd); p,ok=climb(S,open_,0)
    mono=all(d[p[k+1]]==d[p[k]]-1 for k in range(len(p)-1))
    gap=[x for x in p if x in (28,29,34,35)]
    rng=np.random.default_rng(100+sd); vals=list(S.values()); rng.shuffle(vals); Sr=dict(zip(S.keys(),vals))
    pr,okr=climb(Sr,open_,0)
    print(f'{name} 乱数{sd}  着いた{ok} 歩数{len(p)-1}（最短{d[0]}） 一歩ごとに近づく{mono} 通った隙間{gap}  頂点のS {S[0]}  ／負の対照 着いた{okr} 歩数{len(pr)-1}')
