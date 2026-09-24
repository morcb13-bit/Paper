# 検定M1 レジスタどうしの命令と条件分岐で、掛け算・割り算を担体の上のプログラムとして走らせる
#  レジスタ：翼0の軸上の五芒星を下から順に使い、頂点を合わせて主の扇5枚に回して置く（1本＝5桁）
#  命令
#   ADD r,s   r ← r+s（A を r から、B を s から読み出す。足す数も外から与えない）
#   NEG r     鏡で折る：指しを共役（ζ^k→ζ^-k）、表裏はそのまま
#   SHL r / SHR r   扇を一つ上／下へ回す（×5 / ÷5）。SHR で最下位が切れ目を越えたものを余り REM に
#   JNZ r,k / JNEG r,k / JMP k   相対ジャンプ。符号＝上の桁から最初に 0 でない桁の指しが鏡のどちら側か
#  プログラムの列は今は外に持つ（担体に置くのは次の段）
# 事前の基準
#   P1 掛け算 A×B（A∈-40..40, B∈1..30, 積が±1562内）を ADD と JNZ のループで：全組一致
#   P2 割り算 A÷B（A∈0..1500, B∈1..60）を ADD・NEG・JNEG のループで：商と余りが全組一致
#   P3 SHL・SHR：×5 と ÷5（余り -2..2）が 1000 組一致
#   対照 符号判定を殺す（JNEG が跳ばない／JNZ が跳ばない）→ P1・P2 が崩れる
import sys,random
sys.argv=['x']; exec(open('/home/claude/ac1.py').read().split('random.seed(1)')[0])
AX=[]
for q in wing[0]:
    for m in range(10):
        g=U.zsub(q,U.zmul(U.PHI,U.zt(m))); r=ring5(g)
        if r and all(t in wing[0] for t in r) and g not in AX: AX.append(g)
AX.sort(key=lambda g:(U.xy(g)[1],U.xy(g)[0]))
print('翼0の中の五芒星（レジスタの置き場）',len(AX))
def regset(k): return [rotw(AX[k],w) for w in (0,2,4,6,8)]
NAMES={}
def ptr_s(g):
    ptr,sh=state[g]; base=U.zsub(ring5(g)[0],g)
    for e in range(10):
        if U.zmul(base,U.zt(e))==ptr: return e,sh
def setreg(name,X):
    for i,s in enumerate(to_digits(X)):
        g=NAMES[name][i]; state[g]=(U.zmul(U.zsub(ring5(g)[0],g),U.zt(2*s)), -1 if s%2 else 1)
def side(e): return 0 if e==0 else (1 if e in (2,4) else -1)      # 指しが鏡の軸のどちら側か
def sign(name):
    for i in range(4,-1,-1):
        e,_=ptr_s(NAMES[name][i])
        if side(e): return side(e)
    return 0
def getval(name):  # 検算用（判定の外）
    t=0
    for i in range(5):
        e,_=ptr_s(NAMES[name][i]); t+= {0:0,2:1,4:2,6:-2,8:-1}[e]*5**i
    return t
from functools import lru_cache
@lru_cache(maxsize=None)
def digit_rr(i,eA,shA,eB,shB,cin):
    g=D[i]; a=-eA-2*cin; b=eB
    L=[q for q in g['scr'] if rep(ZR,g['PA'][q],a%10)==rep(ZR,g['PB'][q],b%10)]
    lv={(g['dA'][q]-g['dB'][q])//2 for q in L}
    if len(lv)!=1: return None
    s=lv.pop(); sA=shA*(-1 if cin%2 else 1)
    aS=2 if sA==-1 else 0; bS=2 if shB==-1 else 0
    same={rep(TW,g['SA'][q],aS)==rep(TW,g['SB'][q],bS) for q in L}
    if len(same)!=1: return None
    if same.pop(): return s,0
    if s!=0: return s,(-1 if s>0 else 1)
    return s,((1 if cin>0 else -1) if cin else 0)
REM=[0]
def ex(ins):
    op=ins[0]
    if op=='ADD':
        r,s=NAMES[ins[1]],NAMES[ins[2]]; c=0; new=[]
        for i in range(5):
            eA,shA=ptr_s(r[i]); eB,shB=ptr_s(s[i]); v=digit_rr(i,eA,shA,eB,shB,c); s_,c=v; new.append(s_)
        for i,s_ in enumerate(new):
            g=r[i]; state[g]=(U.zmul(U.zsub(ring5(g)[0],g),U.zt(2*s_)), -1 if s_%2 else 1)
    elif op=='NEG':
        for g in NAMES[ins[1]]:
            ptr,sh=state[g]; base=U.zsub(ring5(g)[0],g); e,_=ptr_s(g)
            state[g]=(U.zmul(base,U.zt(-e)),sh)
    elif op in ('SHL','SHR'):
        r=NAMES[ins[1]]; old=[ptr_s(g) for g in r]
        def put(g,e,sh): state[g]=(U.zmul(U.zsub(ring5(g)[0],g),U.zt(e)),sh)
        if op=='SHL':
            for i in range(4,0,-1): put(r[i],*old[i-1])
            put(r[0],0,1)
        else:
            REM[0]={0:0,2:1,4:2,6:-2,8:-1}[old[0][0]]
            for i in range(4): put(r[i],*old[i+1])
            put(r[4],0,1)
def run(prog,maxstep=5000,kill=None):
    pc=0; n=0
    while pc<len(prog) and n<maxstep:
        ins=prog[pc]; n+=1
        if ins[0]=='HALT': return True
        if ins[0]=='JMP': pc+=ins[1]; continue
        if ins[0]=='JNZ':
            pc+= ins[2] if (sign(ins[1])!=0 and kill!='JNZ') else 1; continue
        if ins[0]=='JNEG':
            pc+= ins[2] if (sign(ins[1])<0 and kill!='JNEG') else 1; continue
        ex(ins); pc+=1
    return False
for k,nm in enumerate(['acc','a','b','cnt','m1','one','q','r','nb']): NAMES[nm]=regset(k)
MUL=[('ADD','acc','a'),('ADD','cnt','m1'),('JNZ','cnt',-2),('HALT',)]
DIV=[('ADD','r','nb'),('JNEG','r',3),('ADD','q','one'),('JMP',-3),('ADD','r','b'),('HALT',)]
random.seed(2)
for kill in (None,'JNZ','JNEG'):
    okm=totm=0
    for _ in range(80):
        A=random.randint(-40,40); B=random.randint(1,30)
        if abs(A*B)>1562: continue
        setreg('acc',0); setreg('a',A); setreg('cnt',B); setreg('m1',-1)
        run(MUL,kill=kill); totm+=1; okm+= getval('acc')==A*B
    okd=totd=0
    for _ in range(80):
        A=random.randint(0,600); B=random.randint(1,40)
        setreg('r',A); setreg('b',B); setreg('nb',B); ex(('NEG','nb')); setreg('q',0); setreg('one',1)
        fin=run(DIV,kill=kill); totd+=1; okd+= fin and getval('q')==A//B and getval('r')==A%B
    print(f'{"設計どおり" if kill is None else "対照 "+kill+" を殺す"}: P1 掛け算 {okm}/{totm}   P2 割り算 {okd}/{totd}')
ok3=0
for _ in range(1000):
    X=random.randint(-312,312); setreg('acc',X); ex(('SHL','acc')); a=getval('acc')==5*X
    Z=random.randint(-1562,1562); setreg('acc',Z); ex(('SHR','acc'))
    rr=((Z+2)%5)-2; b=getval('acc')==(Z-rr)//5 and REM[0]==rr
    ok3+= a and b
print(f'P3 ×5 と ÷5（余りつき） {ok3}/1000')
