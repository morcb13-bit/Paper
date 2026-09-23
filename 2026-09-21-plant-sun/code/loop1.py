# 測定LP1 中心のまわりを三角形（扇）を跨いで一周する最短の道の歩数 L
#  網：辺の繋がり＋φ²の繋がり（扇10枚の輪）。巻き数は扇の番号の差（±1,±2、10で一周）を積んで数える
#  見るもの（判定なしの測定。先に何が言えるかを決めておく）
#   L mod 10 が 0 → 一周で番地の位相（一歩 312）が戻る。0 でない → 一周するたびに位相がずれる＝周回が桁上がりを持てる
#   L が中心からの距離に一次で伸びる → 周回も一本の整数の展開
#   巻き +1 と −1 で L が同じか（向きで違えば、回る向きが構造として区別される）
import sys,pickle; sys.path.insert(0,'.')
import b13_chain_units as U
from collections import defaultdict,deque
d=pickle.load(open('/home/claude/prop/g10.pkl','rb')); Q=d['Q']; cw=d['cell_w']; adj={k:set(v) for k,v in d['adj'].items()}
for q in Q: adj.setdefault(q,set())
by=defaultdict(list)
for q in Q:
    p=U.xy(q); by[(int(p[0]//3),int(p[1]//3))].append(q)
for q in Q:
    p=U.xy(q); gx,gy=int(p[0]//3),int(p[1]//3)
    for dx in (-1,0,1):
        for dy in (-1,0,1):
            for r in by.get((gx+dx,gy+dy),()):
                if (min(cw[q])%2)!=(min(cw[r])%2) and U.norm2(U.zsub(q,r))==(2,3): adj[q].add(r); adj[r].add(q)
W={q:min(cw[q]) for q in Q}
def dw(a,b):
    x=(W[b]-W[a])%10; return x-10 if x>5 else x
# 中心からの層
z0=d['z0']; S=set(Q)
c0=[q for q in Q if U.norm2(U.zsub(q,z0))==U.NCELL]
lay={q:0 for q in c0}; dq=deque(c0)
while dq:
    u=dq.popleft()
    for w in adj[u]:
        if w not in lay: lay[w]=lay[u]+1; dq.append(w)
def loop(s,target):
    seen={(s,0):0}; dq=deque([(s,0)])
    while dq:
        u,k=dq.popleft()
        if abs(k)>12: continue
        for w in adj[u]:
            st=(w,k+dw(u,w))
            if st not in seen:
                seen[st]=seen[(u,k)]+1
                if st==(s,target): return seen[st]
                dq.append(st)
rows=[]
for L0 in range(0,40,3):
    cand=[q for q in Q if lay.get(q)==L0 and W[q]==0]
    if not cand: continue
    s=min(cand,key=lambda q:U.xy(q)[0]**2)
    a=loop(s,10); b=loop(s,-10)
    rows.append((L0,a,b))
    print(f'層{L0:2d}  一周+ {a}  一周− {b}   mod4 {a%4}  mod5 {a%5}  mod10 {a%10}')
print('--- 二周（巻き20）と、一周のうち渡った扇の内訳')
def loop2(s,target,lim=25):
    seen={(s,0):0}; par={}; dq=deque([(s,0)])
    while dq:
        u,k=dq.popleft()
        if abs(k)>lim: continue
        for w in adj[u]:
            st=(w,k+dw(u,w))
            if st not in seen:
                seen[st]=seen[(u,k)]+1; par[st]=(u,k)
                if st==(s,target):
                    path=[st]
                    while path[-1] in par: path.append(par[path[-1]])
                    return seen[st],path
                dq.append(st)
    return None,None
from collections import Counter
for L0 in (9,12,15,24):
    s=min([q for q in Q if lay.get(q)==L0 and W[q]==0],key=lambda q:U.xy(q)[0]**2)
    one,p=loop2(s,10); two,_=loop2(s,20)
    jumps=Counter(abs(dw(p[i+1][0],p[i][0])) for i in range(len(p)-1))
    phi2=sum(1 for i in range(len(p)-1) if p[i][0] not in d['adj'].get(p[i+1][0],[]) and p[i+1][0]!=p[i][0] and U.norm2(U.zsub(p[i][0],p[i+1][0]))==(2,3))
    print(f'層{L0}  一周 {one}  二周 {two}  二周−2×一周 {two-2*one}   扇の番号の飛び {dict(jumps)}  φ²の繋がりを使った回数 {phi2}')
print('--- 帯の中だけで一周（層 k〜k+1 から出ない）。最短の一周は中心へ潜る投げ縄だったため')
def loopband(s,target,lo,hi):
    seen={(s,0):0}; dq=deque([(s,0)])
    while dq:
        u,k=dq.popleft()
        if abs(k)>12: continue
        for w in adj[u]:
            if not (lo<=lay.get(w,-1)<=hi): continue
            st=(w,k+dw(u,w))
            if st not in seen:
                seen[st]=seen[(u,k)]+1
                if st==(s,target): return seen[st]
                dq.append(st)
for L0 in range(3,40,3):
    c=[q for q in Q if lay.get(q)==L0 and W[q]==0]
    s=min(c,key=lambda q:U.xy(q)[0]**2)
    a=loopband(s,10,L0,L0+1)
    print(f'層{L0:2d}〜{L0+1}  一周 {a}  {"" if a is None else f"mod4 {a%4} mod5 {a%5} mod10 {a%10}  一周/層 {a/L0:.2f}"}')
print('--- 中心へ潜らせない（層 k 未満を通らない）一周。ずらしの扇を跨ぐことになる')
for L0 in range(3,40,3):
    c=[q for q in Q if lay.get(q)==L0 and W[q]==0]
    s=min(c,key=lambda q:U.xy(q)[0]**2)
    seen={(s,0):0}; par={}; dq=deque([(s,0)]); a=None
    while dq and a is None:
        u,k=dq.popleft()
        if abs(k)>12: continue
        for w in adj[u]:
            if lay.get(w,-1)<L0: continue
            st=(w,k+dw(u,w))
            if st not in seen:
                seen[st]=seen[(u,k)]+1; par[st]=(u,k)
                if st==(s,10): a=seen[st]; break
                dq.append(st)
    if a is None: print(f'層{L0:2d}以上  一周できない'); continue
    path=[(s,10)]
    while path[-1] in par: path.append(par[path[-1]])
    ph=sum(1 for i in range(len(path)-1) if U.norm2(U.zsub(path[i][0],path[i+1][0]))==(2,3))
    mx=max(lay[p[0]] for p in path)
    print(f'層{L0:2d}以上  一周 {a}  mod4 {a%4} mod5 {a%5} mod10 {a%10}  φ²を渡った回数 {ph}  通った最も外の層 {mx}')
