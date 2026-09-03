import json, math, sys
from collections import defaultdict, Counter
import b13_chain_units as U
N=12; K=int(sys.argv[1])
exec(open("or3_exact.py").read().split('starts=[st for R,sts')[0].split('N=12')[1].replace('K=int(sys.argv[1]) if len(sys.argv)>1 else 21',''))
starts=[st for R,sts in SS[:K] for st in sts]
print("五芒星 %d 個／出発 %d 本" % (K,len(starts)))
sizes=Counter(); sigs=set(); acc=0
def dfs(al,d):
    global acc
    if d==N:
        acc+=1; sizes[len(al)]+=1
        sigs.add(hash(frozenset(i for i,_ in al))); return
    for x in (-1,0,1):
        n2=[(i,(e[1],O[e][x])) for i,e in al if x in O[e]]
        if n2: dfs(n2,d+1)
sys.setrecursionlimit(100)
dfs([(i,st) for i,st in enumerate(starts)],0)
n=len(starts); tot=sum(sizes.values())
mean=sum(k*v for k,v in sizes.items())/tot
srt=[]
for k in sorted(sizes): srt += [k]*0
med=None; c=0
for k in sorted(sizes):
    c+=sizes[k]
    if med is None and c>=tot//2: med=k
print("受理語 %d / %d (%.3f%%)｜|W| 最小%d 中央%d 最大%d 平均%.1f｜|W|/出発 平均%.4f｜相異なるW %d 種"
      % (tot,3**N,100*tot/3**N,min(sizes),med,max(sizes),mean,mean/n,len(sigs)))
h=Counter()
for k,v in sizes.items(): h[int(20*k/n)]+=v
print("|W|/出発 の分布 %s" % {"%d%%"%(5*b):h[b] for b in sorted(h)})
