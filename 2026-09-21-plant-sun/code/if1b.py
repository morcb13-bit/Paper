# IF1 残り2行。事前の基準：(1560,0)≡(0,1560)、(1560,1560)≡(0,0) が着地35枚すべてで一致。一枚でも違えばNG
import pickle,sys
sys.argv=['x']; exec(open('slit14.py').read().split('# 鏡の対になる五芒星を選ぶ')[0].replace("import pickle",""))
def rows(sa,sb,label):
    dA=bfs(stars[sa]); dB=bfs(stars[sb])
    tab={}
    for a in (0,1560):
        for b in (0,1560):
            tab[(a,b)]={q:status(dA[q],dB[q],(a-b)%BASE) for q in screen}
    c1=sum(tab[(1560,0)][q]==tab[(0,1560)][q] for q in screen)
    c2=sum(tab[(1560,1560)][q]==tab[(0,0)][q] for q in screen)
    eq=[q for q in screen if dA[q]==dB[q]]
    print(f'{label}: (1560,0)=(0,1560) {c1}/{len(screen)}  (1560,1560)=(0,0) {c2}/{len(screen)}  等距離{len(eq)}個')
    for (a,b),st in tab.items():
        print(f'   入力({a:4},{b:4})  明{sum(v=="B" for v in st.values()):2} 暗{sum(v=="D" for v in st.values()):2}   等距離の着地: '+''.join(st[q] for q in eq))
d=pickle.load(open('/home/claude/prop/slit14.pkl','rb'))
sa,sb=d['sa'],d['sb']
y=U.xy(sa)[1]
other=[h for h in stars if h not in (sa,sb) and abs(U.xy(h)[1]-y)>0.5]
rows(sa,sb,'鏡の対')
rows(sa,other[0],'対照（鏡の対でない）')
