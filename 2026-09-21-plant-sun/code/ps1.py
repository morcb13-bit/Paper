# 検定PS1 位相を定数でも表でもなく、ペンローズの構造だけで運ぶオートマトン
#  二枚（表・裏）：一歩ごとに隣の五角形へ渡り、同時に表裏を入れ替える。裏から表へ戻るとき符号を反転（ねじれ）
#     → 状態 (表,裏) の整数二つ。(u,v)→(-v,u)。4歩で一周（780 を与えない）
#  五芒星のまわり5枚：レジスタの指し先を、一歩ごとに五芒星の中心のまわりで ζ² 回す（Z[ζ] の整数4成分、72°）
#     → 5歩で一周（624 を与えない）
#  入力：半周＝符号反転（二枚の状態に -1）、平衡5進の遅れ＝指し先を ζ^(2·遅れ) 回す
# 事前の基準
#   P1 IF1 の明暗（入力4通り×着地35）がすべて一致   P2 HA1 の表 4/4   P3 RG1 25/25（間隔4の対 0,1,5,7）
#   対照1 ねじれなし（(u,v)→(v,u)）→ P1/P2 が崩れる   対照2 指し先を ζ¹ 回す（10歩で一周）→ P3 が崩れる
#   NG いずれか一つでも一致しない
import pickle
exec(open('br2.py').read().split('# 鏡の対と軸上の五芒星')[0])
pairs=pickle.load(open('/home/claude/prop/br2.pkl','rb'))['pairs']
TW=lambda s:(-s[1],s[0]); NT=lambda s:(s[1],s[0])
def automaton(src,init,op):
    # 源の輪から一歩ずつ広がる：各五角形は親の状態に op を一回かける（最短の層ごと）
    st={q:init for q in stars[src]}; frontier=list(stars[src]); seen=set(frontier)
    while frontier:
        nxt=[]
        for u in frontier:
            for w in adj[u]:
                if w not in seen: seen.add(w); st[w]=op(st[u]); nxt.append(w)
        frontier=nxt
    return st
def sheets(sa,sb,a,b,twist=True):
    op=TW if twist else NT
    A=automaton(sa,(-1 if a else 1,0),op); B=automaton(sb,(-1 if b else 1,0),op)
    out={}
    for q in screen:
        re=A[q][0]+B[q][0]; im=A[q][1]+B[q][1]; I=re*re+im*im
        out[q]='B' if I==4 else ('D' if I==0 else '-')
    return out
sa,sb=pairs[0] if False else pickle.load(open('/home/claude/prop/slit14.pkl','rb'))['sa'],pickle.load(open('/home/claude/prop/slit14.pkl','rb'))['sb']
dA,dB=dist_star[sa],dist_star[sb]
def ref_if1(a,b): return {q:status(dA[q],dB[q],(1560*a-1560*b)%BASE) for q in screen}
for tw,lab in ((True,'ねじれあり'),(False,'対照1 ねじれなし')):
    m=sum(sum(sheets(sa,sb,a,b,tw)[q]==ref_if1(a,b)[q] for q in screen) for a in (0,1) for b in (0,1))
    # HA1：送る/送らない。片方だけのときは片方の状態だけ
    eq=[q for q in screen if dA[q]==dB[q]]; sx=[q for q in screen if (dA[q]-dB[q])%4==2]
    op=TW if tw else NT
    SA=automaton(sa,(1,0),op); SB=automaton(sb,(1,0),op); ha=0
    for a in (0,1):
        for b in (0,1):
            def I(q):
                re=(SA[q][0] if a else 0)+(SB[q][0] if b else 0); im=(SA[q][1] if a else 0)+(SB[q][1] if b else 0); return re*re+im*im
            ha+= all(int(I(q)>0)==(a^b) for q in sx) and all(int(I(q)==4)==(a&b) for q in eq)
    print(f'{lab}: P1 IF1 の明暗一致 {m}/140   P2 HA1 {ha}/4')
# P3：レジスタの指し先を五芒星のまわりで回す
def reg(sa,sb,rot):
    g=sa; base=U.zsub(stars[g][0],g)          # 指し先の最初の向き（五芒星の中心から一枚目へ）
    R=U.zt(rot)
    turn=lambda z:U.zmul(z,R)
    ok=0
    for x in range(-2,3):
        for y in range(-2,3):
            a0=base
            for _ in range((-2*x)%10): a0=turn(a0)      # 遅れ -2x 歩ぶん回しておく
            b0=base
            for _ in range((2*y)%10): b0=turn(b0)
            PA=automaton(sa,a0,turn); PB=automaton(sb,b0,turn)
            lv={(dist_star[sa][q]-dist_star[sb][q])//2 for q in screen if PA[q]==PB[q]}
            exact=any(dist_star[sa][q]-dist_star[sb][q]==2*(x+y) for q in screen)
            c=0 if exact else (1 if x+y>0 else -1)
            if len(lv)==1:
                s=lv.pop()
                if x+y==5*c+s: ok+=1
    return ok
for k in (0,1,5,7):
    s1,s2=pairs[k]
    print(f'対{k} P3 ζ²（5歩で一周）: {reg(s1,s2,2)}/25   対照2 ζ¹（10歩で一周）: {reg(s1,s2,1)}/25')
print('--- 事後に足した対照（対照2が対照になっていなかったため）')
for k in (0,1,5,7):
    s1,s2=pairs[k]
    print(f'対{k}  ζ^5（2歩で一周）: {reg(s1,s2,5)}/25   ζ^4: {reg(s1,s2,4)}/25   ζ^3: {reg(s1,s2,3)}/25')
