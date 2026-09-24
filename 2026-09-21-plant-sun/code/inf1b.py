import collections, inf1 as I, b13_chain_units as U
def edges(T,swap=False):
    c=collections.Counter()
    for _,A,B,C in T:
        for x,y in ((A,B),(A,C),(B,C)): c[frozenset((x,y))]+=1
    return c
for name,sw in (("正",False),("入れ替え",True)):
    T=I.wheel
    for n in range(7): T=I.step(T,sw)
    c=edges(T); over=sum(1 for v in c.values() if v>2)
    tri=collections.Counter(frozenset(t[1:]) for t in T)
    print(name,"7段",len(T),"枚  3枚以上が共有する辺",over,"  重複した三角形",sum(1 for v in tri.values() if v>1))
