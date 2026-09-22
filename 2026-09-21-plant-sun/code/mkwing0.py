# 翼0（105環・五角形721枚）の図形データ wing0.json を geo10.json から作る（動く図が読む）
import json,math
exec(open('maze.py').read().split("def mol_step")[0])
wid={i:n for n,i in enumerate(W)}
pent=sorted({k for i in W for k in R[i][4]}); pid={k:n for n,k in enumerate(pent)}
polys=[[[round(P[k][0]+math.cos(math.radians((P[k][2]+2*m)*36-18)),3),round(P[k][1]+math.sin(math.radians((P[k][2]+2*m)*36-18)),3)] for m in range(5)] for k in pent]
data=dict(poly=polys,rings=[dict(x=round(R[i][1],3),y=round(R[i][2],3),row=R[i][0],p=[pid[k] for k in R[i][4]],n=sorted(wid[j] for j in ADJ[i])) for i in W],
  food=wid[71],seed=[wid[s] for s in (0,1,2)],wallL=[wid[w] for w in range(30,36)],wallR=[wid[w] for w in range(28,34)])
json.dump(data,open('wing0.json','w'),separators=(',',':'))
print('五角形',len(polys),'環',len(data['rings']))
