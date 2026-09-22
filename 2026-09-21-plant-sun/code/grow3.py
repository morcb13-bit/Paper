import json,b13_chain_units as U
R=[[tuple(c) for c in r] for r in json.load(open("R14.json"))]
V=[(-3,2,-1,4),(-5,3,-1,6),(-1,2,1,3),(-2,4,1,5)]
G=[(1,2,2,1),(2,0,2,-1),(3,1,2,-1),(4,-1,2,-3)]; GG=G+[tuple(-t for t in g) for g in G]
X=lambda c:U.xy(c)[0]
def step(rows):
    have=sum(rows,[]); bot=sorted(rows[-1],key=X)
    kids=[]
    for a,b in zip(bot,bot[1:]):                      # 内側：二つの親から一歩で届く点
        c=[U.zadd(a,v) for v in V if U.zsub(U.zadd(a,v),b) in [tuple(x) for x in V] or any(U.zadd(b,w)==U.zadd(a,v) for w in V)]
        c=[x for x in c if U.fits(have+[x]) is not None]
        if len(c)!=1: return None,"内側の候補 %d"%len(c)
        kids.append(c[0])
    for p,side in ((bot[0],-1),(bot[-1],1)):          # 端：一つの親から一歩、兄弟と段の間隔で接し、置ける
        pg=(U.zsub(bot[1],p) if side<0 else U.zsub(p,bot[-2])) if len(bot)>1 else None
        c=[U.zadd(p,v) for v in V if side*(X(U.zadd(p,v))-X(p))>0]
        c=[x for x in c if U.fits(have+kids+[x]) is not None]
        if kids:
            s=kids[0] if side<0 else kids[-1]
            c=[x for x in c if U.zsub(s,x) in GG and (pg is None or (U.zsub(s,x) if side<0 else U.zsub(x,s))!=pg)]
        if len(c)!=1: return None,"端の候補 %d"%len(c)
        kids=[c[0]]+kids if side<0 else kids+[c[0]]
    return rows+[kids],"ok"
rows=[R[0],R[1]]
for n in range(3,15):
    rows,msg=step(rows)
    if rows is None: print(n,"止まった",msg); break
    print(n,"段の円環",len(rows[-1]),"担体と一致",set(rows[-1])==set(R[n-1]))
