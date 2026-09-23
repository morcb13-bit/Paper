# 検定AD3 10枚の輪（主5・ずらし5）で平衡5進の二桁の足し算
# C0（構造）φ² の繋がり（norm2 = (2,3)）が、ずらした扇 w を主の扇 w±1 の両方とだけつなぐ → 10枚の輪が閉じる
# AD3 桁 s は主の扇 2s（mod 10）にいる指し。一歩＝φ² の繋がりで扇を一枚渡る。y を足す＝符号の向きへ 2|y| 歩
#     繰り上がり＝扇5（主4と主6の間のずらし）を渡った回数を、渡った向きの符号つきで数える
#     （輪が閉じれば正しさは数の性質から出る。検定の中身は、輪が閉じることと向きが構造から読めること）
# 事前の基準  OK 625/625   対照 向きを捨てて渡った回数だけ数える（符号を |c| にする）   NG 一つでも違う
import sys,pickle; sys.path.insert(0,'.')
import b13_chain_units as U
from collections import Counter,defaultdict
d=pickle.load(open('/home/claude/prop/g10.pkl','rb')); Q=d['Q']; cw=d['cell_w']
by=defaultdict(list)
for q in Q:
    p=U.xy(q); by[(int(p[0]//3),int(p[1]//3))].append(q)
link=Counter()
for q in Q:
    if not all(w%2 for w in cw[q]): continue
    p=U.xy(q); gx,gy=int(p[0]//3),int(p[1]//3)
    for dx in (-1,0,1):
        for dy in (-1,0,1):
            for r in by.get((gx+dx,gy+dy),()):
                if any(w%2==0 for w in cw[r]) and U.norm2(U.zsub(q,r))==(2,3):
                    for a in cw[q]:
                        for b in cw[r]: link[(a,b)]+=1
print('φ²の繋がり（ずらしの扇, 主の扇）:',sorted(link.items()))
nb=defaultdict(set)
for (a,b) in link: nb[a].add(b); nb[b].add(a)
ring=all(nb[w]=={(w-1)%10,(w+1)%10} for w in range(1,10,2))
print('C0 ずらしの扇が主の扇 w±1 の両方とだけつながる:',ring)
def move(pos,n,signed=True):
    c=0; step=1 if n>0 else -1
    for _ in range(abs(n)):
        nxt=(pos+step)%10
        assert nxt in nb[pos]
        if {pos,nxt}=={4,5}: c+= (1 if step>0 else -1) if signed else 1
        pos=nxt
    return pos,c
def val(pos): return {0:0,2:1,4:2,6:-2,8:-1}[pos]
def add(x0,x1,y0,y1,signed=True):
    p,c0=move((2*x0)%10,2*y0,signed); s0=val(p)
    p,a=move((2*x1)%10,2*y1,signed); p,b=move(p,2*c0,signed); s1=val(p); c1=a+b
    return 25*c1+5*s1+s0
D=range(-2,3)
for sg,lab in ((True,'向きつき'),(False,'対照 向きを捨てる')):
    ok=sum(add(x0,x1,y0,y1,sg)==(5*x1+x0)+(5*y1+y0) for x0 in D for x1 in D for y0 in D for y1 in D)
    print(f'{lab}: {ok}/625')
