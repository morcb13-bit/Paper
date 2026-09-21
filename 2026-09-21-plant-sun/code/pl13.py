#  検定PL13  出発が光より下か上かで、分かれる星の向きが入れ替わるか（13行）
#   OK  出発が光より下 → a が多い、出発が光より上 → b が多い
#   NG  どちらでも同じ向きが多い
#   同じ段の出発は別に数える
from collections import Counter
src=open('pl11.py').read(); exec(src[:src.index('A=Counter()')])
C=Counter(); N=Counter()
for i in range(len(SC)):
    for j in range(len(SC)):
        if i==j: continue
        k='下' if SC[i][1]>SC[j][1]+1e-6 else ('上' if SC[i][1]<SC[j][1]-1e-6 else '同段')
        N[k]+=1
        for o_,ang in rec(i,j): C[(k,o_)]+=1
for k in ('下','上','同段'):
    print(f"出発が光より{k}  配置{N[k]:4d}  a{C[(k,'a')]:5d}  b{C[(k,'b')]:5d}")
a,b=C[('下','a')],C[('下','b')]; c,d=C[('上','a')],C[('上','b')]
print("判定", "OK" if a>b and d>c else "NG")
