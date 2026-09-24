import pickle,random,b13_chain_units as U
d=pickle.load(open("/home/claude/prop/st.pkl","rb")); tabs,tin,C2,D0=d["tabs"],d["tin"],d["C2"],d["D0"]
PHIm1=U.zsub(U.PHI,U.ONE); PHI3=U.zmul(U.zmul(U.PHI,U.PHI),U.PHI); PHI3i=U.zmul(U.zmul(PHIm1,PHIm1),PHIm1)
def upow(n):
    p=U.ONE
    for _ in range(abs(n)): p=U.zmul(p,PHI3 if n>0 else PHI3i)
    return p
def todig(X):
    out=[]
    while X: r=((X+2)%5)-2; out.append(r); X=(X-r)//5
    return out or [0]
def todig_n(X,n):
    d=todig(X); return d+[0]*(n-len(d))
_R={i:U.zadd(C2,U.zmul(upow(i),D0)) for i in range(-40,41)}
def ref2(i): return _R[i]
LO=40
TAB={i:(tabs[i] if i>=0 else tin[-i-1]) for i in range(-LO,41)}
IDX2={ref2(i):i for i in range(-LO,41)}
def fadd(xd,yd,lo,hi,wrong=False):     # 桁 lo..hi−1 に入力、その下 −LO..lo−1 は 0
    rng=list(range(-LO,hi)); X=dict(zip(range(lo,hi),xd)); Y=dict(zip(range(lo,hi),yd))
    f={i:tuple(TAB[i][(X.get(i,0),Y.get(i,0),c)][1] for c in (-1,0,1)) for i in rng}
    A=dict(f); P=dict(f); a,b=1,1; ra=rb=PHI3i; steps=0
    while not all(len(set(v))==1 for v in A.values()) and steps<12:
        s=U.zmul(ra,PHIm1) if wrong else ra; nA={}
        for i in rng:
            j=IDX2.get(U.zadd(C2,U.zmul(s,U.zsub(ref2(i),C2))))
            q=P[j] if j is not None else (-1,0,1)
            nA[i]=tuple(A[i][q[c]+1] for c in range(3))
        P=A; A=nA; a,b=a+b,a; ra,rb=U.zmul(ra,rb),ra; steps+=1
    c={i+1:A[i][1] for i in rng}; c[-LO]=0
    s={i:TAB[i][(X.get(i,0),Y.get(i,0),c[i])][0] for i in range(lo,hi)}
    return steps+1, s, c[hi]
def val(s,lo,hi,cN): return sum(s[i]*5**(i-lo) for i in range(lo,hi))+cN*5**(hi-lo)
def num(d): return sum(v*5**i for i,v in enumerate(d))
random.seed(6)
for nd in (10,20,40):
    xd=[2]+[1]*(nd-1); k,s,cN=fadd(xd,xd,0,nd); kw,sw,cw=fadd(xd,xd,0,nd,wrong=True)
    lim=(5**nd-1)//2; ok=0
    for _ in range(100):
        X,Y=random.randint(-lim,lim),random.randint(-lim,lim); kk,ss,cc=fadd(todig_n(X,nd),todig_n(Y,nd),0,nd); ok+=val(ss,0,nd,cc)==X+Y
    print(f"T2b（止め＝全桁の関数が定数になったら） {nd}桁 貫く入力 和 {val(s,0,nd,cN)==2*num(xd)} 段数 {k}   無作為 {ok}/100   対照（比ずらし）和 {val(sw,0,nd,cw)==2*num(xd)}")
lim=(5**15-1)//2; ok=0; ks=set()
for _ in range(100):
    X,Y=random.randint(-lim,lim),random.randint(-lim,lim)          # 5^10 倍した整数：下10桁が小数部
    kk,ss,cc=fadd(todig_n(X,15),todig_n(Y,15),-10,5); ok+=val(ss,-10,5,cc)==X+Y; ks.add(kk)
print(f"T3b（同じ止め） 整数部5桁＋小数部10桁 和 {ok}/100 段数 {sorted(ks)}")
