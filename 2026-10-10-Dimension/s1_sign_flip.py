# 検定S1：時間の座標だけ符号を入れ替えた長さ²（…＋x_{n-1}²−x_n²）で、中の形の辺を測り直す
# 基準（走らせる前に決めた）：
#   5次元の中の正方形120：時刻の中の24は辺が 8・8 のまま、時間をまたぐ96はすべて一辺が 0 になる → OK
#   対照（3次元の中の線分16）：時間をまたぐ12のうち、0 になるのは一部だけのはず（全部 0 なら手順を疑う）
from itertools import product
def dot(a,b): return sum(x*y for x,y in zip(a,b))
def lor(a): return sum(x*x for x in a[:-1])-a[-1]*a[-1]
def inner(n,k):
    V=set(product((1,-1),repeat=n)); out={}
    for v0 in V:
        D=[tuple(w[i]-v0[i] for i in range(n)) for w in V if w!=v0]
        def grow(ch,st):
            if len(ch)==k:
                pts=set()
                for m in range(1<<k):
                    p=list(v0)
                    for j in range(k):
                        if m>>j&1: p=[x+y for x,y in zip(p,ch[j])]
                    pts.add(tuple(p))
                if pts<=V and dot(ch[0],ch[0])>4: out.setdefault(frozenset(pts),[tuple(c) for c in ch])
                return
            for i in range(st,len(D)):
                a=D[i]
                if ch and dot(a,a)!=dot(ch[0],ch[0]): continue
                if all(dot(a,c)==0 for c in ch): grow(ch+[a],i+1)
        grow([],0)
    return out
from collections import Counter
for n,k in ((3,1),(5,2)):
    S=inner(n,k); tally=Counter()
    for pts,ch in S.items():
        still=len({p[-1] for p in pts})==1
        tally[("時刻の中" if still else "時間をまたぐ", tuple(sorted(lor(c) for c in ch)))]+=1
    print(f"n={n}（中の正{k}次元立方体 {len(S)}）")
    for (g,ls),c in sorted(tally.items()): print(f"   {g}  符号を入れ替えた辺の長さ² {ls} : {c}")
