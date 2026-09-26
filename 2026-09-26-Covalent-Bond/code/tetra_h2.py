# H₂ の正四面体模型：外接球が互いの中心を通る二つの正四面体（中心間＝R）を鏡映で組む。厳密計算（sympy）。
import sympy as sp, json, random
from itertools import combinations
T=[sp.Matrix(v) for v in [(1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1)]]
R2=3; E2=8
v1,v2=T[0],T[1]; g=v1.dot(v2)
al=sp.Rational(3,2)/(3+g); c=v1.cross(v2); be=sp.sqrt((R2-al**2*(6+2*g))/c.dot(c))
u=sp.simplify(al*(v1+v2)+be*c)
B=[sp.simplify(x-2*((x.dot(u)-sp.Rational(3,2))/3)*u) for x in T]
Z=lambda e:sp.simplify(e)==0
d2=lambda p,q:sp.simplify((p-q).dot(p-q))
# 頂点をまとめる：A0..A3、B の新しい頂点
pts=list(T); owner=['S' if i<2 else 'A' for i in range(4)]; mir=[None]*4
for i,b in enumerate(B):
    j=next((k for k,t in enumerate(T) if Z(d2(b,t))),None)
    if j is None: pts.append(b); owner.append('B'); mir[i]=len(pts)-1
    else: mir[i]=j
m=[None]*len(pts)
for i in range(4): m[i]=mir[i]; m[mir[i]]=i
N=len(pts)
E=[(i,j) for i,j in combinations(range(N),2) if Z(d2(pts[i],pts[j])-E2)]
print("頂点",N,"辺",len(E),"共有する頂点",[i for i in range(N) if owner[i]=='S'])
print("中心間²",sp.simplify(u.dot(u)),"  B の頂点と u の距離²",[d2(b,u) for b in B])
# 重ならないか：A はすべて v·u ≤ 3/2、B はすべて ≥ 3/2（等号は共有の2頂点だけ）
sa=[sp.sign(sp.simplify(x.dot(u)-sp.Rational(3,2))) for x in T]; sb=[sp.sign(sp.simplify(x.dot(u)-sp.Rational(3,2))) for x in B]
print("A 側の符号",sa,"  B 側の符号",sb)
# A と B をまたぐ余分な辺がないか
cross=[(i,j) for i,j in E if {owner[i],owner[j]}=={'A','B'}]
print("A と B を直接結ぶ辺",cross)
# 面（三角形）と角の欠損
F=[f for f in combinations(range(N),3) if all(tuple(sorted(p)) in set(E) for p in combinations(f,2)) and (all(owner[k]!='B' for k in f) or all(owner[k]!='A' for k in f))]
print("面",len(F),"  頂点−辺＋面 =",N-len(E)+len(F))
deg={k:sum(k in f for f in F) for k in range(N)}
defect={k:360-60*deg[k] for k in range(N)}
print("角の欠損",{k:(owner[k],defect[k]) for k in range(N)},"  総和",sum(defect.values()))
# 殻：相手の中心までの距離²
for i,p in enumerate(pts):
    oth=u if owner[i]!='B' else sp.zeros(3,1)
    dd=d2(p,oth); ang=sp.N(sp.acos(sp.sqrt(dd)/(2*sp.sqrt(3)))*180/sp.pi,6)
    print(f"  {owner[i]}{i}: 相手の中心までの距離² {dd}  見込む角 {ang}°  鏡像 {m[i]}")
# 鏡像の対応：二回で元へ・辺を辺へ・動かないのは共有頂点だけ
Es=set(E)
print("鏡像：二回で元へ",all(m[m[i]]==i for i in range(N)),"辺を辺へ",all(tuple(sorted((m[a],m[b]))) in Es for a,b in E),"動かない点",[i for i in range(N) if m[i]==i])
# 二つの電子の歩行（整数の擬似乱数）
adj={i:[] for i in range(N)}
for a,b in E: adj[a].append(b); adj[b].append(a)
random.seed(1); a=2; cnt={'A':0,'S':0,'B':0}; meet=0; meetS=0; bad=0
for t in range(200000):
    a=random.choice(adj[a]); b=m[a]
    cnt[owner[a]]+=1
    if a==b: meet+=1; meetS+= owner[a]=='S'
print("20万歩：電子 a の立ち寄り",cnt,"  重なり",meet,"うち共有頂点",meetS)
json.dump({"coords":[[float(sp.N(x)) for x in p] for p in pts],"edges":E,"owner":owner,"mirror":m,
           "OB":[float(sp.N(x)) for x in u],"start":2},open('tetra.json','w'))
