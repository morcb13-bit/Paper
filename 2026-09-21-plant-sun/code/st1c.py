import random
exec(open('st1.py').read().split('random.seed(6)\nfor')[0])
print("元の表の (1,1,1):",tabs[3][(1,1,1)])
lim=(5**10-1)//2
for key,val in (((1,1,1),(0,0)),((0,0,0),(1,0))):
    bad=[dict(t) for t in tabs[:10]]; bad[3][key]=val; random.seed(1)
    P=[(random.randint(-lim,lim),random.randint(-lim,lim)) for _ in range(300)]
    hit=sum(1 for X,Y in P if (todig_n(X,10)[3],todig_n(Y,10)[3])==key[:2])
    print("壊した所",key,"→",val," 和の一致",sum(run(todig_n(X,10),todig_n(Y,10),bad)[1]==X+Y for X,Y in P),"/300  桁3で (x,y) がそこに当たった組",hit)
