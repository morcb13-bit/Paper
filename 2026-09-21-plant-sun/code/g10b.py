import sys; sys.path.insert(0,'.')
src=open('gen10.py').read(); exec(src[:src.index('pid={}')])
from collections import Counter,defaultdict
cell_w=defaultdict(set); cell_a={}
for (w,r,c,k) in rings:
    for q,a in U.ring_cells(place(c,k)):
        cell_w[q].add(w); cell_a.setdefault(q,set()).add(a%10)
Q=list(cell_w)
by=defaultdict(list)
for q in Q:
    p=U.xy(q); by[(int(p[0]//2),int(p[1]//2))].append(q)
adj=defaultdict(set)
for q in Q:
    p=U.xy(q); gx,gy=int(p[0]//2),int(p[1]//2)
    for dx in (-1,0,1):
        for dy in (-1,0,1):
            for r in by.get((gx+dx,gy+dy),()):
                if r!=q and U.norm2(U.zsub(q,r))==U.NCELL: adj[q].add(r)
E=Counter()
for q in Q:
    for r in adj[q]:
        for a in cell_w[q]:
            for b in cell_w[r]:
                if a<b: E[(a,b)]+=1
print('五角形',len(Q),'辺',sum(len(v) for v in adj.values())//2)
print('扇どうしをつなぐ辺（向きつき数え）',sorted(E.items()))
# 連結か
seen={Q[0]}; st=[Q[0]]
while st:
    u=st.pop()
    for w in adj[u]:
        if w not in seen: seen.add(w); st.append(w)
print('連結',len(seen)==len(Q))
print('隣の番地の差',Counter(min((b-a)%10 for a in cell_a[q] for b in cell_a[r]) for q in Q for r in adj[q]))
import pickle; pickle.dump(dict(Q=Q,adj={k:list(v) for k,v in adj.items()},cell_w=dict(cell_w),cell_a=cell_a,z0=z0),open('/home/claude/prop/g10.pkl','wb'))
