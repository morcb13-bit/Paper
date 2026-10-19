# 3 又＋1 個：C9 の正四面体（星型正八面体の P の位置）と、辺を分け合う隣 3 個、誰とも分けない先 1 つ（H9）。整数だけ
from itertools import product, combinations
dot=lambda a,b:sum(x*y for x,y in zip(a,b)); sub=lambda a,b:tuple(x-y for x,y in zip(a,b))
SG=list(product((1,-1),repeat=3))
T={s:[tuple(6*s[k] if i==k else 0 for i in range(3)) for k in range(3)]+[tuple(6*x for x in s)] for s in SG}
C={s:tuple(3*x for x in s) for s in SG}
P=(1,1,-1); TP=set(T[P])
nb=[s for s in SG if s!=P and len(TP&set(T[s]))==2]
pt=[s for s in SG if s!=P and len(TP&set(T[s]))==1]
tip=[v for v in T[P] if all(v not in T[s] for s in SG if s!=P)][0]
print("C9",C[P],"辺を分け合う隣（3 又）",[C[s] for s in nb],"点だけ触れる",[C[s] for s in pt],"誰とも分けない先（H9）",tip)
u=[sub(C[s],C[P]) for s in nb]; h=sub(tip,C[P])
print("C9→隣の長さ²",[dot(x,x) for x in u],"互いの内積",[dot(a,b) for a,b in combinations(u,2)],"（0 は直角）")
print("C9→H の長さ²",dot(h,h),"C9→隣との内積",[dot(h,x) for x in u])
# H9 の軸に沿って見た 3 本：軸に垂直な成分の内積×2＋長さ²＝0 なら 120 度
w=[tuple(dot(h,h)*a-dot(x,h)*b for a,b in zip(x,h)) for x in u]
print("H9 の軸から見た 3 本：内積×2＋長さ²",[2*dot(a,b)+dot(a,a) for a,b in combinations(w,2)],"（0 なら 120 度）")
