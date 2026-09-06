"""LC5-1 τ だけから出る4本の整数系列と一次差分。判定はしない。生データのみ。
   S_k = E_k の連結成分数
   M_k = E_{k+1} の頂点のうち E_k に2本以上接続するものの個数
   P_k = E_k と E_{k+1} の間の辺数（無向）
   I_k = E_k 内部の辺数（無向）
"""
import json, collections
G=json.load(open('carrier_1245_graph.json'))
ADJ=[list(a) for a in G['adj']]; N=len(ADJ)
c=int(open('center.txt').read())
d=[-1]*N; d[c]=0; q=collections.deque([c])
while q:
    u=q.popleft()
    for v in ADJ[u]:
        if d[v]<0: d[v]=d[u]+1; q.append(v)
mx=max(d)
E=[[] for _ in range(mx+1)]
for v in range(N): E[d[v]].append(v)
S=[];M=[];P=[];I=[]
for k in range(mx+1):
    Sk=set(E[k])
    I.append(sum(1 for u in E[k] for w in ADJ[u] if w in Sk)//2)
    P.append(sum(1 for u in E[k] for w in ADJ[u] if d[w]==k+1))
    seen=set(); comp=0
    for u in E[k]:
        if u in seen: continue
        comp+=1; st=[u]; seen.add(u)
        while st:
            x=st.pop()
            for w in ADJ[x]:
                if w in Sk and w not in seen: seen.add(w); st.append(w)
    S.append(comp)
    M.append(sum(1 for w in E[k+1] if sum(1 for x in ADJ[w] if d[x]==k)>=2) if k<mx else 0)
def dump(name,a):
    print(f"\n{name}  (k=0..{len(a)-1}、{len(a)}値)")
    for i in range(0,len(a),16):
        print("  " + " ".join(f"{v:>4}" for v in a[i:i+16]))
    df=[a[i+1]-a[i] for i in range(len(a)-1)]
    print(f"  一次差分 Δ{name}")
    for i in range(0,len(df),16):
        print("  " + " ".join(f"{v:>4}" for v in df[i:i+16]))
for n,a in (("S",S),("M",M),("P",P),("I",I)): dump(n,a)
json.dump({'S':S,'M':M,'P':P,'I':I,'E':[len(e) for e in E]},open('lc5_1.json','w'))
