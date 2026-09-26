# H2 共有結合模型：面で繋いだ正12面体2つの上を歩く電子（整数と加算のみ）
# 数は a+bφ を (a,b) で持つ。φ²=φ+1。
from itertools import product
def add(x,y): return (x[0]+y[0], x[1]+y[1])
def sub(x,y): return (x[0]-y[0], x[1]-y[1])
def mul(x,y):
    a,b=x; c,d=y   # (a+bφ)(c+dφ)=ac+(ad+bc)φ+bdφ² = (ac+bd)+(ad+bc+bd)φ
    return (a*c+b*d, a*d+b*c+b*d)
def dot(u,v):
    s=(0,0)
    for i in range(3): s=add(s,mul(u[i],v[i]))
    return s
Z=(0,0); O=(1,0); P=(0,1); IP=(-1,1)  # 1/φ = φ-1
def neg(x): return (-x[0],-x[1])
V=[]
for s in product([O,neg(O)],repeat=3): V.append(tuple(s))
for sa in [IP,neg(IP)]:
    for sb in [P,neg(P)]:
        base=[Z,sa,sb]
        for k in range(3): V.append(tuple(base[(i-k)%3] for i in range(3)))
V=list(dict.fromkeys(V))
V=[tuple(mul((5,0),c) for c in v) for v in V]  # 5倍に拡大して÷5を整数で閉じさせる
EDGE2=(200,-100)  # 25×(8-4φ)
def d2(u,v):
    w=tuple(sub(u[i],v[i]) for i in range(3)); return dot(w,w)
def edges(Vs): return {(i,j) for i in range(len(Vs)) for j in range(i+1,len(Vs)) if d2(Vs[i],Vs[j])==EDGE2}
EA=edges(V)
# 鏡映：面法線 n=(0,φ,1)、面は x·n = φ²、n·n = 2+φ。÷(2+φ) は ×(3-φ) して ÷5（割り切れを検査）
n=(Z,P,O); PHI2=(5,5); NN=(2,1); CONJ=(3,-1)
def div_nn(x):
    y=mul(x,CONJ); assert y[0]%5==0 and y[1]%5==0; return (y[0]//5,y[1]//5)
def refl(x):
    t=sub(dot(x,n),PHI2); c=div_nn(mul((2,0),t))
    return tuple(sub(x[i],mul(c,n[i])) for i in range(3))
VB=[refl(v) for v in V]
face=[i for i,v in enumerate(V) if dot(v,n)==PHI2]
print("単体 頂点",len(V),"辺",len(EA),"面の頂点",len(face))
# 照合：五角形（長さ5の単純閉路）の数
adjA={i:set() for i in range(20)}
for i,j in EA: adjA[i].add(j); adjA[j].add(i)
def count_cycles(adj,L):
    N=len(adj); c=0
    def dfs(start,cur,path,depth):
        nonlocal c
        for w in adj[cur]:
            if depth==L-1:
                if w==start: c+=1
            elif w>start and w not in path:
                path.add(w); dfs(start,w,path,depth+1); path.remove(w)
    for s in range(N): dfs(s,s,{s},0)
    return c//2
print("単体 長さ5の輪",count_cycles(adjA,5))

def build(shared_idx):
    # A:0..19、B の頂点は A の頂点と座標が一致すれば同一視（shared_idx に入るものだけ）
    VAll=list(V); idB={}
    for k,v in enumerate(VB):
        if v in V and V.index(v) in shared_idx: idB[k]=V.index(v)
        else: idB[k]=len(VAll); VAll.append(("B",k))
    E=set(EA)
    for i,j in edges(VB):
        a,b=sorted((idB[i],idB[j])); E.add((a,b))
    adj={i:set() for i in range(len(VAll))}
    for a,b in E: adj[a].add(b); adj[b].add(a)
    region={i:("S" if i in shared_idx else "A") for i in range(20)}
    for k in idB:
        if idB[k]>=20: region[idB[k]]="B"
    return adj,region,len(E)

def walk(adj,region,start,Lmax):
    cnt={v:0 for v in adj}; cnt[start]=1
    # 閉じた道のうち B を通ったもの：状態に「B を踏んだか」を持つ
    st={(start,0):1}
    rows=[]
    for L in range(1,Lmax+1):
        new={v:0 for v in adj}
        for v in adj:
            for u in adj[v]: new[v]+=cnt[u]
        cnt=new
        ns={}
        for (v,f),c in st.items():
            for u in adj[v]:
                g=f or (region[u]=="B")
                ns[(u,g)]=ns.get((u,g),0)+c
        st=ns
        tot={"A":0,"S":0,"B":0}
        for v,c in cnt.items(): tot[region[v]]+=c
        rows.append((L,tot["A"],tot["S"],tot["B"],st.get((start,0),0),st.get((start,1),0)))
    return rows

# 出発点：共有面から最も遠い A の頂点（対蹠の面の頂点）
far=[i for i,v in enumerate(V) if dot(v,n)==neg(PHI2)]
start=far[0]
# 接触の段階：面一枚／辺一本／頂点一つ／離す
f=face
e=next((i,j) for i,j in EA if i in f and j in f)
cases={"面(5頂点)":set(f),"辺(2頂点)":set(e),"頂点(1)":{f[0]},"離す(0)":set()}
for name,sh in cases.items():
    adj,region,ne=build(sh)
    print(f"\n== {name}  頂点{len(adj)} 辺{ne}")
    print(" L   A側     共有     B側    閉(B不通過) 閉(B通過)")
    for r in walk(adj,region,start,12):
        print(f"{r[0]:2d} {r[1]:7d} {r[2]:7d} {r[3]:7d} {r[4]:9d} {r[5]:9d}")
