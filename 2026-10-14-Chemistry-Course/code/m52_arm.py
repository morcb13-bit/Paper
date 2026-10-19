# 炭素＝シェルピンスキー 2 段目（大きな正四面体の中心に核）。中心から各点までの長さ²（整数だけ、3 次元のまま）
from itertools import combinations
d2=lambda a:sum(x*x for x in a)
BIG=[(6,6,6),(6,-6,-6),(-6,6,-6),(-6,-6,6)]                      # 2 層目の頂点（大きな正四面体の角、辺² 288）
MID=sorted({tuple((a[i]+b[i])//2 for i in range(3)) for a,b in combinations(BIG,2)})   # 辺の中点＝小さな正四面体どうしの共有点
SUB=[tuple((v[i]+sum(m[i] for m in MID if d2(tuple(m[j]-v[j] for j in range(3)))==72))//4 for i in range(3)) for v in BIG]   # 角の小さな正四面体の中心
print("中心→2 層目の頂点（角）",[d2(v) for v in BIG])
print("中心→辺の中点（小さな正四面体の共有点）",sorted({d2(m) for m in MID}),"個数",len(MID))
print("中心→角の小さな正四面体の中心",[d2(s) for s in SUB],SUB)
print("角→その小さな正四面体の中心",[d2(tuple(v[i]-s[i] for i in range(3))) for v,s in zip(BIG,SUB)])
print("小さな正四面体の中心どうし",sorted({d2(tuple(a[i]-b[i] for i in range(3))) for a,b in combinations(SUB,2)}))
