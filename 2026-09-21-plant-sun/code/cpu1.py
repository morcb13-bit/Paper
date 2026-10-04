# 仮想のプロセッサ一式（担体なし）：プログラム内蔵・フェッチ・解読・PC・分岐・NEG・SHL/SHR
#  事前の基準（観察者側の整数の模型と照合）
#   T-mul  掛け算 A×B（繰り返しの足し算、JNZ で回る）            40通り 停止して積が一致
#   T-fib  フィボナッチ n 回（LOAD/STORE/ADD/JMP/JNZ）             30通り 停止して F(n) が一致
#   T-abs  絶対値（JNEG と NEG）                                    50通り
#   T-sh   ×25（SHL 2回）と ÷5（SHR）                              50通り
#   対照1 解読で命令の位相を鏡に通す  対照2 PC を進めない（400手で止まらない）  対照3 符号の読みを逆
#  構造で持つもの：番地＝位相の組でベクトルを回した点／解読＝命令の位相でベクトルを回した点に置いた装置
#   PC の＋1＝同じ加算器／分岐先＝命令の番地の位相をそのまま PC へ／NEG＝鏡／SHL・SHR＝φ^(∓3) の相似で隣の層から移す
#   符号＝上の桁から「0なら素通し、0でなければ自分の符号」を、φ^(+3a) の相似で合成（相手がどこにも無くなったら止める）
#  規則として残すもの：位相の符号の読み（指しが鏡のどちら側か）、繰り上がりの規則、置いた装置の種類そのもの
import random
exec(open('mem1.py').read().split('# ---- メモリ')[0])
_PC={}
def partner(i,r,up=False):
    key=(i,r,up)
    if key not in _PC:
        _PC[key]=IDX.get(U.zadd(C2,U.zmul(r,U.zsub(REF[i],C2))))
    return _PC[key]
def add_step(acc,opB):                    # 相手の引きを覚えておくだけ速くした版（中身は同じ）
    rng=range(-LO,HI); R=dict(acc); B=dict(opB)
    for i in range(-LO,0): R[i]=(0,0); B[i]=(0,0)
    f={i:tuple(D(i,*R[i],*B[i],c)[0] for c in (-1,0,1)) for i in rng}
    A=dict(f); P=dict(f); rs=[PHI3i,PHI3i]; st=0
    while not all(len(set(v))==1 for v in A.values()):
        ra=RA[st]; nA={}
        for i in rng:
            j=partner(i,ra); q=P[j] if j is not None else (-1,0,1)
            nA[i]=tuple(A[i][q[c]+1] for c in range(3))
        P=A; A=nA; st+=1
    cin={i:(A[i-1][1] if i>0 else A[-1][1]) for i in range(HI)}
    return {i:D(i,*R[i],*B[i],cin[i])[1] for i in range(HI)}
RA=[]; a_,b_=PHI3i,PHI3i
for _ in range(12): RA.append(a_); a_,b_=U.zmul(a_,b_),a_
UP=[]; a_,b_=PHI3,PHI3
for _ in range(10): UP.append(a_); a_,b_=U.zmul(a_,b_),a_
PH={((-2*r)%10,(-2*r)%4):r for r in range(-2,3)}
def phase(r): return ((-2*r)%10,(-2*r)%4)
def mirror(p): return ((-p[0])%10,(-p[1])%4)
Z=(0,0); ONE=phase(1)
M=5**HI
def wrap(X): return ((X+(M-1)//2)%M)-(M-1)//2
def todig(X,n):
    out=[]
    while X: r=((X+2)%5)-2; out.append(r); X=(X-r)//5
    return (out+[0]*n)[:n]
def val(reg): return sum(PH[reg[i]]*5**i for i in range(HI))
def reg_of(X): return {i:phase(r) for i,r in enumerate(todig(wrap(X),HI))}
# ---- 番地：語の二つの番地の桁の位相 (k0,k1) で E0・E1 を回して足した点 ----
E0=U.zmul(upw(20),D0); E1=U.zmul(upw(22),U.zmul(U.zt(1),D0))
def wpts(p0,p1): 
    T=U.zadd(U.zrot(E0,p0[0]),U.zrot(E1,p1[0])); return [U.zadd(REF[i],T) for i in range(HI)]
ALLW=[(phase(a),phase(b)) for a in range(-2,3) for b in range(-2,3)]
assert len({p for w in ALLW for p in wpts(*w)})==25*HI
def addr_of(n): d=todig(n,2); return (phase(d[0]),phase(d[1]))
# ---- 解読：命令の上二桁の位相 (k8,k9) で V0・V1 を回して足した点に、装置を置く ----
V0=U.zmul(upw(5),D0); V1=U.zmul(upw(6),U.zmul(U.zt(1),D0))
def opt(p8,p9): return U.zadd(U.zrot(V0,p8[0]),U.zrot(V1,p9[0]))
CODES={'LOAD':(1,0),'STORE':(2,0),'ADD':(-1,0),'NEG':(-2,0),'JMP':(0,1),'JNZ':(1,1),'JNEG':(2,1),'SHL':(1,-1),'SHR':(-1,-1),'HALT':(0,0)}
UNIT={opt(phase(r8),phase(r9)):name for name,(r8,r9) in CODES.items()}
assert len(UNIT)==len(CODES)
def instr(name,n=0):
    d=todig(n,2); r8,r9=CODES[name]; w=[0]*HI; w[0],w[1],w[8],w[9]=d[0],d[1],r8,r9
    return {i:phase(w[i]) for i in range(HI)}
# ---- 符号：上から合成（0 の桁＝素通し、0 でない桁＝自分の符号）。φ^(+3a) の相似で上の相手を引く ----
def sgn_of_phase(p,flip=False):
    s=0 if p==Z else (1 if p[0] in (6,8) else -1)          # 指しが鏡のどちら側か（規則として残す）
    return -s if flip else s
def sign(reg,flip=False):
    f={i:((-1,0,1) if reg[i]==Z else (sgn_of_phase(reg[i],flip),)*3) for i in range(HI)}
    A=dict(f); P=dict(f); st=0
    while True:
        r=UP[st]; nA={}; found=False
        for i in range(HI):
            j=partner(i,r,True)
            if j is not None and j<HI: found=True; q=P[j]
            else: q=(-1,0,1)
            nA[i]=tuple(q[A[i][c]+1] for c in range(3))
        if not found: break
        P=A; A=nA; st+=1
    return A[0][1]
# ---- 一手 ----
def shift(reg,up):                         # SHL：一つ下の層から、SHR：一つ上の層から、相似で移す
    r=PHI3i if up else PHI3; out={}
    for i in range(HI):
        j=partner(i,r); out[i]=reg[j] if (j is not None and 0<=j<HI) else Z
    return out
def run(mem,ctl='ok',limit=400):
    acc={i:Z for i in range(HI)}; pc={i:Z for i in range(HI)}
    for t in range(limit):
        w={i:mem[p] for i,p in enumerate(wpts(pc[0],pc[1]))}              # フェッチ
        p8,p9=(mirror(w[8]),mirror(w[9])) if ctl=='dec' else (w[8],w[9])
        op=UNIT.get(opt(p8,p9))                                           # 解読
        if op is None: return None
        ap=wpts(w[0],w[1])
        if op=='HALT': return acc
        nxt=pc if ctl=='pc' else add_step(pc,{i:(mirror(ONE) if i==0 else Z) for i in range(HI)})
        if op=='LOAD': acc={i:mem[p] for i,p in enumerate(ap)}
        elif op=='STORE':
            for i,p in enumerate(ap): mem[p]=acc[i]
        elif op=='ADD': acc=add_step(acc,{i:mirror(mem[p]) for i,p in enumerate(ap)})
        elif op=='NEG': acc={i:mirror(acc[i]) for i in range(HI)}
        elif op=='SHL': acc=shift(acc,True)
        elif op=='SHR': acc=shift(acc,False)
        elif op=='JMP': nxt={i:(w[i] if i<2 else Z) for i in range(HI)}
        elif op=='JNZ' and sign(acc)!=0: nxt={i:(w[i] if i<2 else Z) for i in range(HI)}
        elif op=='JNEG' and sign(acc,ctl=='sgn')<0: nxt={i:(w[i] if i<2 else Z) for i in range(HI)}
        pc=nxt
    return None
def load(prog,data):
    mem={p:Z for w in ALLW for p in wpts(*w)}
    for n,(name,a) in enumerate(prog):
        for i,p in enumerate(wpts(*addr_of(n))): mem[p]=instr(name,a)[i]
    for n,X in data.items():
        for i,p in enumerate(wpts(*addr_of(n))): mem[p]=reg_of(X)[i]
    return mem
def rd(mem,n): return val({i:mem[p] for i,p in enumerate(wpts(*addr_of(n)))})
MUL=[('LOAD',-2),('JNZ',3),('HALT',0),('ADD',-4),('STORE',-2),('LOAD',-3),('ADD',-1),('STORE',-3),('JMP',0)]
FIB=[('LOAD',-1),('JNZ',3),('HALT',0),('ADD',-5),('STORE',-1),('LOAD',-2),('ADD',-3),('STORE',-4),('LOAD',-3),('STORE',-2),('LOAD',-4),('STORE',-3),('JMP',0)]
ABS=[('LOAD',-1),('JNEG',3),('HALT',0),('NEG',0),('HALT',0)]
SH=[('LOAD',-1),('SHL',0),('SHL',0),('STORE',-2),('LOAD',-1),('SHR',0),('HALT',0)]
def fib(n):
    a,b=0,1
    for _ in range(n): a,b=b,a+b
    return a
def shr(X): d=todig(wrap(X),HI)[0]; return (wrap(X)-d)//5
def trial(ctl):
    random.seed(21); res={}
    ok=0
    for _ in range(40):
        A=random.randint(-10**5,10**5); B=random.randint(0,10)
        m=load(MUL,{-1:A,-2:B,-3:0,-4:-1}); r=run(m,ctl); ok+= r is not None and rd(m,-3)==wrap(A*B)
    res['T-mul']=f"{ok}/40"; ok=0
    for _ in range(30):
        n=random.randint(1,12); m=load(FIB,{-1:n,-2:0,-3:1,-4:0,-5:-1}); r=run(m,ctl); ok+= r is not None and rd(m,-2)==fib(n)
    res['T-fib']=f"{ok}/30"; ok=0
    for _ in range(50):
        X=random.randint(-(M-1)//2,(M-1)//2); m=load(ABS,{-1:X}); r=run(m,ctl); ok+= r is not None and val(r)==abs(X)
    res['T-abs']=f"{ok}/50"; ok=0
    for _ in range(50):
        X=random.randint(-(M-1)//2,(M-1)//2); m=load(SH,{-1:X}); r=run(m,ctl); ok+= r is not None and val(r)==shr(X) and rd(m,-2)==wrap(25*X)
    res['T-sh']=f"{ok}/50"
    return res
for ctl,name in (('ok','本番'),('dec','対照1 解読で鏡を通す'),('pc','対照2 PCを進めない'),('sgn','対照3 符号の読みを逆')):
    print(name, trial(ctl), flush=True)
