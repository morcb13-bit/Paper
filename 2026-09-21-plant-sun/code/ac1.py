# 検定AC1 アキュムレータ：読み出し → 5桁加算器（PC5）→ 書き戻し を担体の中で N 回まわす
#  レジスタ：翼0の軸上の五芒星一つを選び、頂点を合わせて ζ^w で回して主の扇5枚に置く（相対アドレス）
#   桁 s の持ち方：指し＝そのまわり5枚のうち ζ^(2s) の向きの一枚、表裏＝(-1)^s
#  書き戻し：着地の段 s → 指しを ζ^(2s)、表裏を (-1)^s に置く（私が決めた写しの規則）
#  読み出し：指しを鏡で折って（共役 ζ^k→ζ^-k）A の送り出しの初めの指しにする、表裏もそのまま A の初めに
#  足す数 Y は B の遅れで毎回同じに与える。外から与えるのは最初の値 R0 と Y だけ
# 事前の基準
#   OK  無作為な R0,Y（各200組）で N=1..40 回まわした後のレジスタが R0+N·Y（5⁵=3125 で折り返し、-1562..1562）と一致
#   対照 書き戻しを外す（毎回 R0 を読む）→ N≥2 で一致しない
#   NG  一つでも食い違う
import sys,pickle,random
exec(open('/home/claude/pc5.py').read().split('# D1')[0])
st0=[g for g in stars_all] if False else None
# 翼0の軸上の五芒星（レジスタの置き場）
x0=U.xy(rotw(sa0,0))[0]
S0=[]
for q in wing[0]:
    for m in range(10):
        g=U.zsub(q,U.zmul(U.PHI,U.zt(m)))
        r=ring5(g)
        if r and all(t in wing[0] for t in r) and abs(U.xy(g)[0]-(U.xy(sa0)[0]+U.xy(sb0)[0])/2)<1e-6: S0.append(g)
reg0=sorted(set(S0),key=lambda g:U.xy(g)[1])[0]
REG=[rotw(reg0,w) for w in (0,2,4,6,8)]
print('レジスタの五芒星（翼0の軸上）を5枚の扇に回して置いた:',len(set(REG)),'個、すべて五芒星:',all(ring5(g) for g in REG))
state={}   # 五芒星 → (指し Z[ζ], 表裏 ±1)
def write(i,s): 
    g=REG[i]; state[g]=(U.zmul(U.zsub(ring5(g)[0],g),U.zt(2*s)), -1 if s%2 else 1)
def read(i):
    g=REG[i]; ptr,sh=state[g]; base=U.zsub(ring5(g)[0],g)
    # 指しを鏡で折って向きの差を取り出す：ptr·conj(base)/|base|² の代わりに、ζ の冪を回して一致を探す構造のまま使う
    for e in range(10):
        if U.zmul(base,U.zt(e))==ptr: return e, sh
def digit_reg(i,e,sh,y,cin):
    # A の初めの指し：送り出しを ζ^(-e) と繰り上がり分 ζ^(-2cin) だけ遅らせる（鏡で折った向き）
    g=D[i]; a=-e-2*cin; b=2*y
    L=[q for q in g['scr'] if rep(ZR,g['PA'][q],a%10)==rep(ZR,g['PB'][q],b%10)]
    lv={(g['dA'][q]-g['dB'][q])//2 for q in L}
    if len(lv)!=1: return None
    s=lv.pop()
    # 表裏：A の初めの表裏は sh（(-1)^x）と繰り上がりの半周 (-1)^cin
    shA=sh*(-1 if cin%2 else 1); shB=-1 if y%2 else 1
    aS=2 if shA==-1 else 0
    same={ (rep(TW,g['SA'][q],aS)==rep(TW,g['SB'][q],b%4)) for q in L}
    if len(same)!=1: return None
    if same.pop(): return s,0
    if s!=0: return s,(-1 if s>0 else 1)
    return s,((1 if cin>0 else -1) if cin else 0)
def to_digits(X):
    X=((X+1562)%3125)-1562; out=[]
    for _ in range(5):
        r=((X+2)%5)-2; out.append(r); X=(X-r)//5
    return out
def val(ds): return sum(ds[i]*5**i for i in range(5))
def run(R0,Y,N,writeback=True):
    for i,s in enumerate(to_digits(R0)): write(i,s)
    Yd=to_digits(Y)
    for _ in range(N):
        c=0; new=[]
        for i in range(5):
            e,sh=read(i); v=digit_reg(i,e,sh,Yd[i],c)
            if v is None: return None
            s,c=v; new.append(s)
        if writeback:
            for i,s in enumerate(new): write(i,s)
    return val(new)
random.seed(1); ok=bad=0; okc=0; tot=0
for _ in range(200):
    R0=random.randint(-1562,1562); Y=random.randint(-1562,1562)
    for N in (1,2,3,7,15,40):
        want=((R0+N*Y+1562)%3125)-1562; tot+=1
        ok+= run(R0,Y,N)==want
        okc+= run(R0,Y,N,False)==want
print(f'AC1 書き戻しあり {ok}/{tot}   対照 書き戻しなし {okc}/{tot}（N=1 の {tot//6} 組は一致して当然）')
