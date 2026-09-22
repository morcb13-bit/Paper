# 検定PI0  一歩の向きの ζ^k を掛けて、種（翼0の頂点）から各環まで全部の親子の鎖を足す（Z[ζ10] の4整数）
#  OK（打ち消しあり）振幅が ζ^k×経路数 より小さくなる環がある ／ NG（打ち消しなし）全環でちょうど ζ^k×経路数
import json,b13_chain_units as U
G=json.load(open('geo10.json'));R=G['R']
w0=[i for i,r in enumerate(R) if r[5]==0]
R14=[[tuple(c) for c in r] for r in json.load(open('R14.json'))]
cen=[c for row in R14 for c in sorted(row,key=lambda c:U.xy(c)[0])]
V=[(-3,2,-1,4),(-5,3,-1,6),(-1,2,1,3),(-2,4,1,5)]; K={V[0]:3,V[1]:3,V[2]:2,V[3]:2}   # 左下 ζ^3・右下 ζ^2（短長とも）
A={w0[0]:(1,0,0,0)}; C={w0[0]:1}; bad=0; zero=0
for i in w0[1:]:
  a=(0,0,0,0); c=0
  for p in R[i][3]:
    d=U.zsub(cen[i],cen[p]); a=U.zadd(a,U.zrot(A[p],K[d])); c+=C[p]
  A[i]=a; C[i]=c
  if not any(U.zrot((c,0,0,0),k)==a for k in range(10)): bad+=1
  if a==(0,0,0,0): zero+=1
print('環',len(w0),'振幅＝ζ^k×経路数でないもの',bad,'打ち消して0',zero,'最大経路数',max(C.values()))
