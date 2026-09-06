"""LC5-2 局所状態の重複検定。意味付けはしない。一致を数えるだけ。"""
import json, collections, random, math

def seqs(adj, src, N):
    d=[-1]*N; d[src]=0; q=collections.deque([src])
    while q:
        u=q.popleft()
        for v in adj[u]:
            if d[v]<0: d[v]=d[u]+1; q.append(v)
    mx=max(d); E=[[] for _ in range(mx+1)]
    for v in range(N): E[d[v]].append(v)
    S=[];M=[];P=[];I=[]
    for k in range(mx+1):
        Sk=set(E[k])
        I.append(sum(1 for u in E[k] for w in adj[u] if w in Sk)//2)
        P.append(sum(1 for u in E[k] for w in adj[u] if d[w]==k+1))
        seen=set(); c=0
        for u in E[k]:
            if u in seen: continue
            c+=1; st=[u]; seen.add(u)
            while st:
                x=st.pop()
                for w in adj[x]:
                    if w in Sk and w not in seen: seen.add(w); st.append(w)
        S.append(c)
        M.append(sum(1 for w in E[k+1] if sum(1 for x in adj[w] if d[x]==k)>=2) if k<mx else 0)
    return S,M,P,I

def dups(tuples):
    c=collections.Counter(tuples)
    rep={t:n for t,n in c.items() if n>1}
    pairs=sum(n*(n-1)//2 for n in c.values())
    return rep, pairs

def null(cols, trials=1000, seed=1):
    rnd=random.Random(seed); out=[]
    for _ in range(trials):
        sh=[rnd.sample(c,len(c)) for c in cols]
        out.append(dups(list(zip(*sh)))[1])
    return out

def report(name, S,M,P,I, lo, hi):
    lv=[(S[k],M[k],P[k],I[k]) for k in range(lo,hi)]
    dv=[(S[k+1]-S[k],M[k+1]-M[k],P[k+1]-P[k],I[k+1]-I[k]) for k in range(lo,hi-1)]
    print(f"\n=== {name}  レコード {len(lv)} 本（k={lo}..{hi-1}）===")
    for lab,tp in (("水準 (S,M,P,I)",lv),("差分 (ΔS,ΔM,ΔP,ΔI)",dv)):
        rep,pairs=dups(tp)
        cols=list(map(list,zip(*tp)))
        nl=null(cols)
        nm=sum(nl)/len(nl)
        ge=sum(1 for x in nl if x>=pairs)/len(nl)
        print(f"  {lab}: 一致する組 {pairs} 組 / 重複した値 {len(rep)} 種")
        print(f"      並べ替えの帰無：平均 {nm:.2f} 組、実測以上が出た割合 {ge:.3f}")
        if rep:
            for t,n in sorted(rep.items(), key=lambda x:-x[1])[:12]:
                pos=[k+lo for k in range(len(tp)) if tp[k]==t]
                print(f"      {t} ×{n}  位置 k={pos}")
    # 部分一致（3つ組・2つ組）も数える
    import itertools
    for r in (3,2):
        for cb in itertools.combinations(range(4),r):
            nm_=["S","M","P","I"]
            sub=[tuple(t[i] for i in cb) for t in dv]
            _,p=dups(sub)
            print(f"      差分 部分 {''.join('Δ'+nm_[i] for i in cb)}: {p} 組")

D=json.load(open('lc5_1.json'))
report("ペンローズ担体", D['S'],D['M'],D['P'],D['I'], 0, 126)   # 最終層は除外

# 対照：正方格子
S_=60
def sid(x,y): return (x+S_)*(2*S_+1)+(y+S_)
adj=[[] for _ in range((2*S_+1)**2)]
for x in range(-S_,S_+1):
    for y in range(-S_,S_+1):
        for dx,dy in ((1,0),(-1,0),(0,1),(0,-1)):
            if -S_<=x+dx<=S_ and -S_<=y+dy<=S_: adj[sid(x,y)].append(sid(x+dx,y+dy))
s,m,p,i=seqs(adj,sid(0,0),len(adj))
report("対照 正方格子", s,m,p,i, 0, min(len(s)-1,S_))
