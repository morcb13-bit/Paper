# ダイヤモンドを炭素 1 個から育てる：結合を 1 段ずつたどる（世代）
from itertools import product
B=list(product([0,1],repeat=3)); odd=lambda n:sum(n)%2==1
SEED=(9,3,3)                                  # A の小立方体 (1,0,0) の中心
DA=[(3,3,3),(3,-3,-3),(-3,3,-3),(-3,-3,3)]    # A の席の炭素から出る 4 本
DB=[tuple(-x for x in d) for d in DA]         # B の席の炭素から出る 4 本（向きが反転）
isA=lambda p:((p[0]//6)+(p[1]//6)+(p[2]//6))%2==1 and all(x%6==3 for x in p)
gen={SEED:0}; front=[SEED]; shells=[[SEED]]
for k in range(1,9):
    nxt=[]
    for p in front:
        for d in (DA if isA(p) else DB):
            q=tuple(p[i]+d[i] for i in range(3))
            if q not in gen: gen[q]=k; nxt.append(q)
    shells.append(nxt); front=nxt
print("世代ごとの炭素",[len(s) for s in shells])
print("A の席か",[set(isA(p) for p in s) for s in shells])
cA=0;cB=0;outA=[];outB=[]
for k,s in enumerate(shells):
    if k%2==0: cA+=len(s); outA.append(cA)
    else: cB+=len(s); outB.append(cB)
print("A の正四面体の数（累計、偶数世代）",outA)
print("B の席の炭素（累計、奇数世代）",outB)
fib=[1,1];[fib.append(fib[-1]+fib[-2]) for _ in range(20)]; print("フィボナッチ数か",[(x, x in fib) for x in outA])
# 結合 4 本・どの炭素も同じか（内部）
allp=set(gen); print("内側の炭素の結合数",set(sum(1 for d in (DA if isA(p) else DB) if tuple(p[i]+d[i] for i in range(3)) in allp) for p in allp if gen[p]<=6))
# 立方体 24 の中で：ダイヤモンドの A の正四面体と、シェルピンスキー 3 段目の正四面体
def sier(o,e,lev):
    if lev==1: return [o]
    h=e//2; return [t for n in B if odd(n) for t in sier(tuple(o[i]+n[i]*h for i in range(3)),h,lev-1)]
S3=set(sier((0,0,0),24,3))        # 3 段目：辺 6 の小立方体の原点 16 個
D=set(tuple(6*n[i] for i in range(3)) for n in product(range(4),repeat=3) if odd(n))  # ダイヤモンドの A：奇の小立方体すべて
print("立方体24：シェルピンスキー3段目",len(S3),"ダイヤモンドのA",len(D),"3段目はダイヤモンドに含まれる",S3<=D)
rest=D-S3; big=set(tuple((x//12) for x in o) for o in rest); print("残り",len(rest),"が入る大きな小立方体（辺12）",sorted(big),"偶か",[not odd(b) for b in sorted(big)])
