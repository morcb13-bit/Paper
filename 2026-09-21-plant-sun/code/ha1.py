# 検定HA1 半加算器（出口二つ）。入力＝スリットから送る(1)/送らない(0)
# 振幅：歩数 d を4で割った余りで (1,0)(0,1)(-1,0)(0,-1)（一歩780、BASE3120）。強度＝re²+im²、整数のみ
# 出口：和＝歩数差が4で割って2余る着地（ずれなしで暗の着地）、強度>0 を1
#       桁上がり＝等距離の着地、強度=4 を1
# 事前の基準
#  OK  入力4通りすべてで、和の出口すべてが XOR、桁上がりの出口すべてが AND
#  NG  一つでも食い違う
#  対照1 鏡の対でない五芒星二つ（出口は同じ規則で選ぶ）
#  対照2 片方の入力の位相を 780（四分の一周期）ずらす → 表が崩れるはず
#  対照3 各入力に 0〜3119 の無作為な位相を足す（20回）→ 大半が崩れるはず
import sys,random,pickle
exec(open('slit14.py').read().split('# 鏡の対になる五芒星を選ぶ')[0].replace("import pickle",""))
AMP=[(1,0),(0,1),(-1,0),(0,-1)]
def amp(d,off):
    ph=(STEP*d+off)%BASE
    if ph%STEP: return None                 # 780 の倍数でない位相は整数の四値に載らない
    return AMP[ph//STEP]
def run(sa,sb,offA=0,offB=0,label='',quiet=False):
    dA=bfs(stars[sa]); dB=bfs(stars[sb])
    sumX=[q for q in screen if (dA[q]-dB[q])%4==2]
    carX=[q for q in screen if dA[q]==dB[q]]
    ok=True; lines=[]
    for a in (0,1):
        for b in (0,1):
            def I(q):
                re=im=0
                for on,d,off in ((a,dA[q],offA),(b,dB[q],offB)):
                    if on:
                        v=amp(d,off)
                        if v is None: return None
                        re+=v[0]; im+=v[1]
                return re*re+im*im
            s=[I(q) for q in sumX]; c=[I(q) for q in carX]
            sbit=[None if v is None else int(v>0) for v in s]
            cbit=[None if v is None else int(v==4) for v in c]
            good=all(x==(a^b) for x in sbit) and all(x==(a&b) for x in cbit)
            ok&=good and len(sumX)>0 and len(carX)>0
            lines.append(f'   入力({a},{b})  和の強度{s} 桁上がりの強度{c}  → 和{set(sbit)} 桁{set(cbit)}  {"OK" if good else "NG"}')
    if not quiet:
        print(f'{label}: 和の出口{len(sumX)}個 桁上がりの出口{len(carX)}個  → {"OK" if ok and len(sumX) and len(carX) else "NG"}')
        print('\n'.join(lines))
    return ok
d=pickle.load(open('/home/claude/prop/slit14.pkl','rb')); sa,sb=d['sa'],d['sb']
y=U.xy(sa)[1]; other=[h for h in stars if h not in (sa,sb) and abs(U.xy(h)[1]-y)>0.5]
run(sa,sb,label='鏡の対')
run(sa,other[0],label='対照1 鏡の対でない')
run(sa,sb,0,780,label='対照2 片方を四分の一周期ずらす')
random.seed(1); n=sum(run(sa,sb,random.randrange(BASE),random.randrange(BASE),quiet=True) for _ in range(20))
print('対照3 無作為な位相 20回のうち表どおり',n)
