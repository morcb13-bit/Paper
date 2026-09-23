# 検定AD2 二桁の足し算（平衡5進、各桁 -2..2 → X,Y ∈ -12..12、625通り）を構造だけのオートマトンで
#  一歩：指し先を ζ¹ 回す（番地一つ）＋表裏を入れ替え、裏→表で符号反転
#  一桁：間隔4の鏡の対。指し先が一致する着地の段 d/2 が桁 s
#  繰り上がり（構造で読む）：その着地で表裏の状態も一致すれば 0、一致しなければ -sign(s)
#     （表裏は d/2 の偶奇を持つので、和と s の偶奇が違う＝5 ずれた、と読める。私の導出）
#  次の桁へ送る：繰り上がり c を、次の対の A の指し先を ζ^(-2c) 余分に回すことで渡す
# 事前の基準
#   OK  625通りすべてで X+Y = 25·c1 + 5·s1 + s0
#   予測 二桁目の和が ±5（例 2+2+1）のとき s=0 で符号が決まらず落ちる。数は 2+2+1 と -2-2-1 の場合の数だけ
#   NG  予測以外で一つでも食い違う
#   対照 表裏のねじれを外す
import pickle
exec(open('ps1.py').read().split('sa,sb=pairs[0] if False')[0])
pairs=pickle.load(open('/home/claude/prop/br2.pkl','rb'))['pairs']
def rot(z,k):
    for _ in range(k%10): z=U.zmul(z,U.zt(1))
    return z
def twn(s,k,op):
    for _ in range(k%4): s=op(s)
    return s
def digit(pr,x,y,twist=True):
    sa,sb=pr; base=U.zsub(stars[sa][0],sa); op=TW if twist else NT
    PA=automaton(sa,rot(base,-2*x),lambda z:U.zmul(z,U.zt(1)))
    PB=automaton(sb,rot(base,2*y),lambda z:U.zmul(z,U.zt(1)))
    SA=automaton(sa,twn((1,0),-2*x,op),op); SB=automaton(sb,twn((1,0),2*y,op),op)
    L=[q for q in screen if PA[q]==PB[q]]
    lv={(dist_star[sa][q]-dist_star[sb][q])//2 for q in L}
    if len(lv)!=1: return None
    s=lv.pop(); same={SA[q]==SB[q] for q in L}
    if len(same)!=1: return None
    if same.pop(): return s,0
    if s==0: return s,None
    return s,(-1 if s>0 else 1)
def add(X0,X1,Y0,Y1,twist=True):
    r0=digit(pairs[0],X0,Y0,twist)
    if r0 is None or r0[1] is None: return None
    s0,c0=r0
    r1=digit(pairs[1],X1+c0,Y1,twist)
    if r1 is None or r1[1] is None: return None
    s1,c1=r1; return 25*c1+5*s1+s0
D=range(-2,3)
for tw,lab in ((True,'ねじれあり'),(False,'対照 ねじれなし')):
    ok=0; amb=0; bad=[]
    for x0 in D:
     for x1 in D:
      for y0 in D:
       for y1 in D:
        r=add(x0,x1,y0,y1,tw); X=5*x1+x0; Y=5*y1+y0
        if r==X+Y: ok+=1
        else:
            s0=(x0+y0); c0=0 if -2<=s0<=2 else (1 if s0>0 else -1)
            if abs(x1+y1+c0)==5: amb+=1
            else: bad.append((X,Y,r))
    print(f'{lab}: 正しい {ok}/625   二桁目の和が±5で落ちたもの {amb}   それ以外で食い違い {len(bad)}  例{bad[:3]}')
