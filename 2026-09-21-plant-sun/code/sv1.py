# 検定SV1 同じ機械で、絶対番地の篩（素数）と相対番地の篩（幸運数）を走らせる
#  語（5桁）：主の扇の五芒星36個＋ずらしの扇の五芒星36個（頂点を合わせて回し、ずらし t を加えて置く）＝72語
#   番地 0..65 をプログラムと変数、66..71 を内部レジスタに使う
#  印の置き場：翼0の五角形のうち、どの五芒星も囲んでいない五角形を行の順に並べ、番地 1..N を割り当てる（N=300）
#  命令 0 HALT 1 LDA 2 STA 3 ADD 4 MRK（印[ACC]を立てる）5 TST（ACC←印[ACC]）6 NEG（鏡）7 JNZ 8 JNEG 9 JMP
# 事前の基準
#   V1 素数のプログラム：終了後、2..300 で印の立っていない数が既知の素数と一致
#   V2 幸運数のプログラム：終了後、1..300 で印の立っていない数が既知の幸運数と一致
#   対照 素数のプログラムを、印の番地を「生き残りの中の順位」で読む相対番地のまま走らせる → 素数にならない
#   NG 一つでも食い違う
import sys
exec(open('/home/claude/sp1.py').read().split("MULP=")[0].replace("MEM=[regset(k) for k in range(30)]","MEM=None").replace("for k,nm in zip(range(30,36),['ACC','PC','IR','OFF','ONE','Z']): NAMES[nm]=regset(k)",""))
t=(-11,4,-4,11)
def oddset(k): return [U.zadd(rotw(AX[k],w),U.zrot(t,w-1)) for w in (1,3,5,7,9)]
W=[regset(k) for k in range(36)]+[oddset(k) for k in range(36)]
allst=[g for w in W for g in w]
print('語',len(W),' 五芒星の重なりなし',len(set(allst))==len(allst),' すべて五芒星',all(ring5(g) for g in allst))
exec(open('/home/claude/fast.py').read())
MEM=W[:66]
for k,nm in zip(range(66,72),['ACC','PC','IR','OFF','ONE','Z']): NAMES[nm]=W[k]
used=set(p for g in allst for p in ring5(g))
S0all=[]
for q in wing[0]:
    for m in range(10):
        g=U.zsub(q,U.zmul(U.PHI,U.zt(m))); r=ring5(g)
        if r: used|=set(r)
free=sorted([q for q in wing[0] if q not in used],key=lambda q:(round(U.xy(q)[1],4),round(U.xy(q)[0],4)))
N=300; CELL=free[:N+1]; mark={}
print('印の置き場（五芒星を囲まない五角形）',len(free),'→ 使う',N)
MODE=['abs']
def cell_of(a):
    if MODE[0]=='abs': return CELL[a] if 0<=a<=N else None
    alive=[i for i in range(1,N+1) if not mark.get(CELL[i])]           # 相対番地：生き残りの中の a 番目
    return CELL[alive[a-1]] if 1<=a<=len(alive) else None
def machine2(maxstep=2_000_000):
    store(NAMES['PC'],0); store(NAMES['ONE'],1); store(NAMES['ACC'],0)
    for n in range(maxstep):
        pc=value(NAMES['PC'])
        copy(MEM[pc],NAMES['IR']); ir=NAMES['IR']
        op=dval(ptr_s(ir[4])[0])*5+dval(ptr_s(ir[3])[0]); arg=sum(dval(ptr_s(ir[i])[0])*5**i for i in range(3))
        copy(ir[:3],NAMES['OFF'][:3]); [state.__setitem__(g,(U.zsub(ring5(g)[0],g),1)) for g in NAMES['OFF'][3:]]
        jump=False
        if op==0: return True
        elif op==1: copy(MEM[arg],NAMES['ACC'])
        elif op==2: copy(NAMES['ACC'],MEM[arg])
        elif op==3: addregs(NAMES['ACC'],MEM[arg])
        elif op==4:
            c=cell_of(value(NAMES['ACC']))
            if c is not None: mark[c]=1
        elif op==5: c=cell_of(value(NAMES['ACC'])); store(NAMES['ACC'],1 if (c is not None and mark.get(c)) else 0)
        elif op==6:
            for g in NAMES['ACC']:
                e,sh=ptr_s(g); state[g]=(U.zmul(U.zsub(ring5(g)[0],g),U.zt(-e)),sh)
        elif op==7: jump=sgn(NAMES['ACC'])!=0
        elif op==8: jump=sgn(NAMES['ACC'])<0
        elif op==9: jump=True
        addregs(NAMES['PC'], NAMES['OFF'] if jump else NAMES['ONE'])
    return False
def asm(src,vars_):
    lab={}; code=[]
    for line in src:
        if line.endswith(':'): lab[line[:-1]]=len(code)
        else: code.append(line.split())
    OPS={'HALT':0,'LDA':1,'STA':2,'ADD':3,'MRK':4,'TST':5,'NEG':6,'JNZ':7,'JNEG':8,'JMP':9}
    out=[]
    for i,c in enumerate(code):
        op=OPS[c[0]]; arg=0
        if len(c)>1: arg= (lab[c[1]]-i) if op in (7,8,9) else vars_[c[1]]
        out.append((op,arg))
    return out
PRIME=['TOP:','LDA v','ADD NN','JNEG GO','HALT','GO:','LDA v','TST','JNZ NEXT','LDA v','ADD v','STA j',
       'LOOP:','LDA j','ADD NN','JNEG MK','JMP NEXT','MK:','LDA j','MRK','ADD v','STA j','JMP LOOP',
       'NEXT:','LDA v','ADD ONE','STA v','JMP TOP']
PV=dict(v=40,j=41,NN=42,ONE=43)
LUCKY=['L0:','LDA ONE','STA POS','LDA ZERO','STA C',
       'SCAN:','LDA POS','ADD NN','JNEG S1','JMP FIND','S1:','LDA POS','TST','JNZ SKIP',
       'LDA C','ADD ONE','STA C','ADD NK','JNZ SKIP','LDA POS','MRK','LDA ZERO','STA C',
       'SKIP:','LDA POS','ADD ONE','STA POS','JMP SCAN',
       'FIND:','LDA ONE','STA POS','LDA ZERO','STA C',
       'F2:','LDA POS','ADD NN','JNEG F3','HALT','F3:','LDA POS','TST','JNZ FSKIP',
       'LDA C','ADD ONE','STA C','ADD NRR','JNZ FSKIP',
       'LDA POS','NEG','STA NK','LDA NRR','ADD M1','STA NRR','JMP L0',
       'FSKIP:','LDA POS','ADD ONE','STA POS','JMP F2']
LV=dict(POS=52,C=53,NN=54,ONE=55,ZERO=56,NK=57,NRR=58,M1=59)
def loadp(prog):
    for a in range(66): store(MEM[a],0)
    for a,(op,arg) in enumerate(prog): store(MEM[a],word(op,arg))
def primes_ref(n): return [p for p in range(2,n+1) if all(p%d for d in range(2,int(p**.5)+1))]
def lucky_ref(n):
    L=list(range(1,n+1,2)); i=1
    while i<len(L) and L[i]<=len(L):
        k=L[i]; L=[x for j,x in enumerate(L) if (j+1)%k]; i+=1
    return L
P=asm(PRIME,PV); Lp=asm(LUCKY,LV)
print('プログラムの長さ 素数',len(P),'語  幸運数',len(Lp),'語')
for mode,lab in (('abs','V1 素数（絶対番地）'),('rel','対照 素数のプログラムを相対番地で')):
    MODE[0]=mode; mark.clear(); loadp(P); store(MEM[40],2); store(MEM[42],-(N+1)); store(MEM[43],1)
    f=machine2(); MODE[0]='abs'
    got=[i for i in range(2,N+1) if not mark.get(CELL[i])]
    print(f'{lab}: 停止={f}  一致={got==primes_ref(N)}  個数 {len(got)}（素数 {len(primes_ref(N))}）  先頭 {got[:12]}')
MODE[0]='abs'; mark.clear(); loadp(Lp)
for k,v in dict(NN=-(N+1),ONE=1,ZERO=0,NK=-2,NRR=-2,M1=-1).items(): store(MEM[LV[k]],v)
f=machine2()
got=[i for i in range(1,N+1) if not mark.get(CELL[i])]
print(f'V2 幸運数（相対番地を数えて除く）: 停止={f}  一致={got==lucky_ref(N)}  個数 {len(got)}（幸運数 {len(lucky_ref(N))}）  先頭 {got[:14]}')
