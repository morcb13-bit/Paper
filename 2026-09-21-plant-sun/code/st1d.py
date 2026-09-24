exec(open('st1.py').read().split('random.seed(6)\nfor')[0])
from collections import defaultdict
g=defaultdict(list)
for x in range(-2,3):
    for y in range(-2,3):
        f=tuple(ref[(x,y,c)][1] for c in (-1,0,1)); g[f].append(x+y)
for f,v in sorted(g.items()): print("c_in −1,0,+1 → c_out",f,"  x+y",sorted(set(v)))
