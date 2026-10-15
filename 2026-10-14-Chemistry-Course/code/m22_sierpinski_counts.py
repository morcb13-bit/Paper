import itertools
from collections import Counter
def kids(o,e):  # kept subcubes: odd index sum
    h=e//2
    return [tuple(o[i]+h*n[i] for i in range(3)) for n in itertools.product([0,1],repeat=3) if sum(n)%2==1]
def tetA(o,e): return [tuple(o[i]+e*n[i] for i in range(3)) for n in itertools.product([0,1],repeat=3) if sum(n)%2==1]
cubes=[((0,0,0),48)]
for k in range(1,4):
    cubes=[(c,e//2) for o,e in cubes for c in kids(o,e)]
    cnt=Counter(v for o,e in cubes for v in tetA(o,e))
    print(k,len(cubes),"verts",len(cnt),"shared",sum(1 for v in cnt.values() if v==2),"single",sum(1 for v in cnt.values() if v==1),"max",max(cnt.values()))
# Check level1 tetra union equals Sierpinski: kids of big A
