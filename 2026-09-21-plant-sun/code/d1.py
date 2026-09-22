# 検定D1（走らせる前に決めた）  匂いが壁を越えて届くか
#  OK  乱数3通りすべてで、地面の通れる環すべてに S>0、かつ 壁の向こう（段0〜6）の S の和が 層1枚より大きい
#  参考（判定に使わない）：地面の S を種から登ったときに止まる場所
exec(open('smell.py').read())
def climb(S,open_,start):
  cur=start;p=[cur]
  while cur!=FOOD:
    nb=[j for j in ADJ[cur] if j in open_]; b=max(nb,key=lambda j:S[j])
    if S[b]<=S[cur]: break
    cur=b;p.append(cur)
  return p
open_=Wset-set(WALL_L); far=[i for i in open_ if R[i][0]<=6]
for L in (1,3):
  for sd in (1,2,3):
    r=run(WALL_L,L=L,T=20000,seed=sd,body=False); S=r['S']
    z=sum(1 for i in open_ if S[i]==0); p=climb(S,open_,0)
    print(f'層{L} 乱数{sd}  S=0の環{z}  壁の向こうのSの和{sum(S[i] for i in far)}  餌のS{S[FOOD]}  頂点から登って止まる環{p[-1]}（段{R[p[-1]][0]}）')
