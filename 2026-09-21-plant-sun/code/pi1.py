# 検定PI1（登録済み）  経路に限定して均等に広がる餌を、整数と加算だけで書く
#   D_t(i) = Σ(通れる隣 j の D_{t-1}(j)) + （餌なら 4^t）      ＝ 長さ n の道の寄与 (1/4)^n を 4^t 倍したもの
#   壁と担体の外へ出た分は戻らない
#  OK  種から「最も濃い通れる隣へ移る」を繰り返すと、最短歩数ちょうどで餌に着く
#  対照1 D を壁を無視して広げる → 壁の前で止まるはず   対照2 隙間を右へ → 右の隙間を通る最短路
exec(open('maze.py').read().split("def mol_step")[0])
FOOD=71; SEED=(0,1,2)
WALL_L=list(range(30,36)); WALL_R=list(range(28,34))
def field(spread,T):
  D={i:0 for i in W}; p4=1
  for t in range(T):
    p4*=4
    D={i:(sum(D[j] for j in ADJ[i] if j in spread) if i in spread else 0) for i in W}
    D[FOOD]+=p4
  return D
def climb(D,open_,start):
  cur=start;p=[cur]
  while cur!=FOOD:
    nb=[j for j in ADJ[cur] if j in open_]; m=max(D[j] for j in nb); B=[j for j in nb if D[j]==m]
    if m<=D[cur]: return p,'止まる'
    if len(B)>1: return p,'同点%d'%len(B)
    cur=B[0];p.append(cur)
  return p,'着く'
for name,wall,ignore in (('本番 左の隙間',WALL_L,False),('対照1 壁を無視して広げる',WALL_L,True),('対照2 右の隙間',WALL_R,False)):
  open_=Wset-set(wall); d=bfs(FOOD,open_)
  for T in (100,300):
    D=field(Wset if ignore else open_,T)
    out=[]
    for s in SEED:
      p,st=climb(D,open_,s)
      ok=st=='着く' and len(p)-1==d[s]
      out.append(f'種{s}:{st} 歩数{len(p)-1}/最短{d[s]} 隙間{[x for x in p if x in (28,29,34,35)]} 止まった段{R[p[-1]][0]} {"OK" if ok else "NG"}')
    print(name,f'T={T}',' ｜ '.join(out))
