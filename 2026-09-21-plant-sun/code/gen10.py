# 10枚の担体：72°ずつの5枚＋その間に36°回して外へずらした5枚（ずらす量はただ一通り t）
import json,cmath,b13_chain_units as U
R=[[tuple(c) for c in r] for r in json.load(open("R14.json"))]
z0=(2,-2,0,-3); t=(-11,4,-4,11)
V=[(-3,2,-1,4),(-5,3,-1,6),(-1,2,1,3),(-2,4,1,5)]
X=lambda c:U.xy(c)[0]
def place(c,k): 
    p=U.zadd(U.zrot(U.zsub(c,z0),k),z0)
    return U.zadd(p,U.zrot(t,k-1)) if k%2 else p
rings=[];idx={}
for w,k in enumerate(range(10)):
    for r,row in enumerate(R):
        for c in sorted(row,key=X):
            idx[(w,c)]=len(rings); rings.append((w,r,c,k))
pid={};P=[];out=[]
for i,(w,r,c,k) in enumerate(rings):
    par=[idx[(w,U.zsub(c,v))] for v in V if (w,U.zsub(c,v)) in idx and rings[idx[(w,U.zsub(c,v))]][1]==r-1] if r>0 else []
    cr=place(c,k); pl=[]
    for q,a in U.ring_cells(cr):
        if q not in pid: pid[q]=len(P); x,y=U.xy(q); P.append([round(x,4),round(y,4),a%10])
        pl.append(pid[q])
    x,y=U.xy(cr); out.append([r,round(x,4),round(y,4),par,pl,w])
assert U.fits([place(c,k) for (w,r,c,k) in rings]) is not None
json.dump(dict(rot=round(cmath.phase(U.ROT),6),P=P,R=out,center=list(U.xy(z0))),open('geo10.json','w'),separators=(',',':'))
print("環",len(out),"五角形",len(P))
