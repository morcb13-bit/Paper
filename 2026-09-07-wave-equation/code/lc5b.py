import json, collections
G=json.load(open('carrier_1245_graph.json'))
ADJ=[list(a) for a in G['adj']]; N=len(ADJ)
def bfs(s):
    d=[-1]*N; p=[-1]*N; d[s]=0; q=collections.deque([s])
    while q:
        u=q.popleft()
        for v in ADJ[u]:
            if d[v]<0: d[v]=d[u]+1; p[v]=u; q.append(v)
    return d,p
d,_=bfs(0); a=max(range(N),key=lambda i:d[i])
d,_=bfs(a); b=max(range(N),key=lambda i:d[i])
d2,p=bfs(b)
path=[a]
while path[-1]!=b: path.append(p[path[-1]])
c=path[len(path)//2]
print(f"最遠点対 ({a},{b}) 歩数={d2[a]}  経路の中点 id={c}")
for s,lab in ((c,"中点"),(a,"端点a"),(b,"端点b")):
    dd,_=bfs(s); print(f"  {lab} id={s}: 最大τ={max(dd)}")
open('center.txt','w').write(str(c))
