# 検定RG1 五芒星のレジスタ一桁（平衡5進）に二つの桁の和を書き、繰り上がりを出す
# 仕組み（整数のみ）
#   一歩の位相 624（5歩で一周、BASE 3120）。入力の桁 x,y ∈ {-2..2} を送り出しの遅れにする：A は -2x 歩、B は +2y 歩
#   着地 q で二つの位相が一致（624·(dA-dB) ≡ 624·(2y+2x) mod 3120）する段の d/2 をレジスタの桁 s とする
#   繰り上がり c：歩数そのものが一致（dA-dB = 2(x+y)）する着地があれば 0、無ければ x+y の符号
# 事前の基準
#   OK   25通りすべてで、一致する着地の d/2 がただ一つの値 s に揃い、x+y = 5c + s、s ∈ {-2..2}
#   NG   一致する段が二つ以上ある／無い／x+y≠5c+s が一つでもある
#   対照1 一歩 780（4歩で一周） 対照2 間隔の広い鏡の対（歩数10）  → どちらも崩れるはず
import pickle
exec(open('br2.py').read().split('# 鏡の対と軸上の五芒星')[0])
pairs=pickle.load(open('/home/claude/prop/br2.pkl','rb'))['pairs']
def reg(sa,sb,step=624):
    dA,dB=dist_star[sa],dist_star[sb]; ok=0; bad=[]
    per=BASE//step
    for x in range(-2,3):
        for y in range(-2,3):
            lv={(dA[q]-dB[q])//2 for q in screen if (step*((dA[q]-dB[q])-2*(x+y)))%BASE==0}
            exact=any(dA[q]-dB[q]==2*(x+y) for q in screen)
            c=0 if exact else (1 if x+y>0 else -1)
            if len(lv)==1:
                s=lv.pop()
                if x+y==5*c+s and -2<=s<=2: ok+=1; continue
            bad.append((x,y,sorted(lv) if isinstance(lv,set) else lv,c))
    return ok,bad
for k,(sa,sb) in enumerate(pairs):
    sep=min(dist_star[sa][q] for q in stars[sb])
    ok,bad=reg(sa,sb)
    line=f'対{k:2d} 間隔{sep:2d}  一歩624: {ok}/25'
    if sep==4:
        ok2,_=reg(sa,sb,780); line+=f'   対照1 一歩780: {ok2}/25'
    print(line+('   例 '+str(bad[:2]) if bad else ''))
