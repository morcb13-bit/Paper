#  検定PL10  二葉はどちらの向きの星で起きるか
#   測る   全ての（星i・光=星j, j≠i）で、i のまわり5枚の最大が二つになるか
#          i の向き a（最上段の星と同じ向き）/ b ごとに件数
#   読み   件数を星の数で揃えて整数で比べる（cA*nB と cB*nA）
#          片方の向きでだけ起きる → 二葉は向きで決まる
#          両方で起きる        → 向きだけでは決まらない
#   再現   13行の総数は v257 の 98組
#   空試験 光が一様 → 5つ同点で二葉にならない
import math, itertools, os
from collections import Counter
exec(open('wind_core.py').read())
exec(open('pentagram_neuron.py').read().split('def run_tests')[0])
src=open('amoeba.py').read(); exec(src[src.index('def build'):src.index('def run_amoeba')])
src=open('assign_zero.py').read(); exec(src[src.index('def star_tips'):src.index('def fuda')])
SC, AROUND, _ = build()
F, faces, _SC = carrier()
ST=star_tips(faces)
tip={ (round(c[0],6),round(c[1],6)): round(min(a)%72,3) for c,a in ST}
top=min(range(len(SC)), key=lambda i: SC[i][1])
o=lambda i: tip[(round(SC[i][0],6),round(SC[i][1],6))]
ori={i: 'a' if o(i)==o(top) else 'b' for i in range(len(SC))}
def nmax(i,peak):
    v=[round(-((p[0]-peak[0])**2+(p[1]-peak[1])**2),6) for q,p in AROUND[i]]
    return v.count(max(v))
cnt=Counter(); tot=0
for i in range(len(SC)):
    for j in range(len(SC)):
        if i!=j and nmax(i,SC[j])==2: cnt[ori[i]]+=1; tot+=1
n=Counter(ori.values())
flat=all(len(AROUND[i])==5 for i in range(len(SC)))  # 一様光では値がすべて同じ→5同点
print(f"行 {os.environ['B13_ROWS']}  星 a{n['a']}/b{n['b']}  二葉の組 {tot}  a{cnt['a']} b{cnt['b']}"
      f"  星あたり比較 a:{cnt['a']*n['b']} b:{cnt['b']*n['a']}  空試験 {'OK' if flat else 'NG'}")
