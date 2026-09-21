#  検定PL9（担体を大きくして）  二葉のあと両方の枝が伸び続ける配置はあるか
#   走らせ方  全ての（出発i・光j, i≠j）で grow(T=12)
#   OK   分岐のあと両方の枝がそれぞれ3コマ以上伸びる配置が一つ以上ある
#   NG   どの配置でも、分岐のあと片方または両方が3コマ未満で止まる
#   あわせて：分岐した星の向き a/b、発火の有無
import math, itertools, os
from collections import Counter, defaultdict
src=open('plant.py').read()
exec(src[:src.index('print("13/φ²')])
src=open('assign_zero.py').read(); exec(src[src.index('def star_tips'):src.index('def fuda')])
F, faces, _ = carrier(); ST=star_tips(faces)
tip={(round(c[0],6),round(c[1],6)): round(min(a)%72,3) for c,a in ST}
top=min(range(len(SC)), key=lambda i: SC[i][1])
o=lambda i: tip[(round(SC[i][0],6),round(SC[i][1],6))]
ori={i:'a' if o(i)==o(top) else 'b' for i in range(len(SC))}
runs=0; withsplit=0; fire=0; sp_ori=Counter(); good=[]; depth=Counter()
for i in range(len(SC)):
    for j in range(len(SC)):
        if i==j: continue
        runs+=1
        hist,fires,splits,addr,live=grow(i,SC[j],T=12,verbose=False)
        if fires: fire+=1
        if not splits: continue
        withsplit+=1
        # 実際に二本へ分かれた分岐だけ（行き先が異なる）
        for (t,s,br) in splits:
            kids=[h for h in hist if h[0]==t and h[1]==s and h[3]!=br]
            if len({h[2] for h in kids})<2: continue
            sp_ori[ori[s]]+=1
            nbs=[h[3] for h in kids]
            # 子孫の枝番号：nb から派生（br*2+1+w の木）
            def grown(nb):
                ts=set()
                for h in hist:
                    b=h[3]
                    while b>nb: b=(b-1)//2
                    if b==nb and h[0]>t: ts.add(h[0])
                return len(ts)
            g=min(grown(nb) for nb in nbs)
            depth[min(g,3)]+=1
            if g>=3: good.append((i,j,s,t,g))
print(f"行 {os.environ['B13_ROWS']}  配置 {runs}  発火あり {fire}  分岐あり {withsplit}"
      f"  二本に分かれた分岐の星 a{sp_ori['a']} b{sp_ori['b']}"
      f"  分岐後の短い側が伸びたコマ数 {dict(sorted(depth.items()))}（3は3以上）  → {'OK' if good else 'NG'}")
for g in good[:5]: print("   出発%d 光%d 分岐星%d(%s) %dコマ目 両枝%dコマ以上" % (g[0],g[1],g[2],ori[g[2]],g[3],g[4]))
