"""検定LC5 τ-等高集合の内部構造
     入力は隣接だけ。座標・角度・距離は読まない。
     出すもの：|E_k| / E_k 内部辺 / E_k→E_{k+1} 辺数 / 分裂 / 合流
       分裂 = E_k を内部辺だけで分けたときの連結成分の個数
       合流 = E_{k+1} の頂点のうち、E_k に隣を2つ以上持つものの個数
     必ず落ちる設定：出発点を変えて表が全く同じなら、この表は τ を見ていない
"""
import json, collections, sys
G=json.load(open('carrier_1245_graph.json'))
ADJ=[list(a) for a in G['adj']]; N=len(ADJ)

def layers(src):
    d=[-1]*N; d[src]=0; q=collections.deque([src])
    while q:
        u=q.popleft()
        for v in ADJ[u]:
            if d[v]<0: d[v]=d[u]+1; q.append(v)
    return d

def table(src):
    d=layers(src)
    mx=max(d)
    E=[[] for _ in range(mx+1)]
    for v in range(N): E[d[v]].append(v)
    rows=[]
    for k in range(mx+1):
        S=set(E[k])
        inner=sum(1 for u in E[k] for w in ADJ[u] if w in S)//2
        fwd=sum(1 for u in E[k] for w in ADJ[u] if d[w]==k+1)
        # 連結成分（内部辺のみ）
        seen=set(); comp=0
        for u in E[k]:
            if u in seen: continue
            comp+=1; st=[u]; seen.add(u)
            while st:
                x=st.pop()
                for w in ADJ[x]:
                    if w in S and w not in seen: seen.add(w); st.append(w)
            # 連結成分
        merge=0
        if k+1<=mx:
            for w in E[k+1]:
                if sum(1 for x in ADJ[w] if d[x]==k)>=2: merge+=1
        rows.append((k,len(E[k]),inner,fwd,comp,merge))
    return rows

src=int(sys.argv[1]) if len(sys.argv)>1 else None
if src is None:
    # 出発点は隣接だけで決める：次数2の頂点のうち id 最小 …ではなく、単に 0
    src=0
rows=table(src)
print(f"出発点 id={src}   最大τ={rows[-1][0]}")
print(f"{'k':>4} {'|E_k|':>7} {'内部辺':>7} {'→E_k+1':>8} {'分裂':>5} {'合流':>6}")
for r in rows:
    print(f"{r[0]:>4} {r[1]:>7} {r[2]:>7} {r[3]:>8} {r[4]:>5} {r[5]:>6}")
