# 静止図：二層のシェルピンスキーから、ベンゼンとトリプチセンへ（段階を追う 6 枚）
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from itertools import product, combinations
for f in fm.findSystemFonts():
    if "NotoSansCJK" in f.replace(" ","") and "Regular" in f: fm.fontManager.addfont(f)
plt.rcParams["font.family"]=["Noto Sans CJK JP","Noto Sans CJK TC","sans-serif"]
SG=list(product((1,-1),repeat=3))
TET={s:[tuple(6*s[k] if i==k else 0 for i in range(3)) for k in range(3)]+[tuple(6*x for x in s)] for s in SG}
BIGA=[(6,6,6),(6,-6,-6),(-6,6,-6),(-6,-6,6)]; BIGB=[tuple(-x for x in v) for v in BIGA]
CEN=[tuple(3*x for x in s) for s in SG]
AX=[(1,1,1),(-1,-1,1),(-1,1,-1),(1,-1,-1)]
d2=lambda a,b:sum((x-y)**2 for x,y in zip(a,b))
HEX=[]
for ax in AX:
    rest=[c for c in CEN if abs(sum(x*y for x,y in zip(c,ax)))!=9]; r=[rest[0]]
    while len(r)<6: r.append([c for c in rest if c not in r and d2(c,r[-1])==36][0])
    HEX.append(r)
HC=["#7a3fd1","#c0392b","#1e8a5a","#c27c0e"]
P,Q,H9,H10=(3,3,-3),(-3,-3,3),(6,6,-6),(-6,-6,6)
PATHS=[[P,(3,3,3),(3,-3,3),Q],[P,(3,-3,-3),(-3,-3,-3),Q],[P,(-3,3,-3),(-3,3,3),Q]]; PH=[3,2,0]
def basis(view,up):
    v=np.array(view,float); v/=np.linalg.norm(v); u=np.array(up,float); u-=u.dot(v)*v; u/=np.linalg.norm(u); return np.cross(u,v),u,v
def fig(fn,title,note,view,up,rings,ring_alpha,paths=False,labels=(),spread=False):
    x,y,v=basis(view,up); pr=lambda p:(np.dot(p,x),np.dot(p,y),np.dot(p,v))
    f,ax=plt.subplots(figsize=(6,6.9),dpi=150); ax.set_aspect("equal"); ax.axis("off")
    it=[]
    L=lambda a,b,c,w,al,ls="-",dz=0:it.append(((pr(a)[2]+pr(b)[2])/2+dz,"l",(pr(a),pr(b),c,w,al,ls)))
    for B,c in ((BIGA,"#2f6f8f"),(BIGB,"#b0652a")):
        for a,b in combinations(B,2): L(a,b,c,1.1,0.45)
    for s in SG:
        c="#2f6f8f" if s[0]*s[1]*s[2]>0 else "#b0652a"
        for a,b in combinations(TET[s],2): L(a,b,c,0.7,0.22)
    for h in range(4):
        al=ring_alpha[h]
        if al<=0: continue
        off=np.array([np.cos(h*1.57+0.4),np.sin(h*1.57+0.4)])*0.28 if spread else np.zeros(2)
        for i in range(6):
            a,b=pr(HEX[h][i]),pr(HEX[h][(i+1)%6])
            a=(a[0]+off[0],a[1]+off[1],a[2]); b=(b[0]+off[0],b[1]+off[1],b[2])
            it.append(((a[2]+b[2])/2+0.2,"l",(a,b,HC[h],3.2 if spread else (4.2 if al>0.5 else 1.6),al,"-" if al>0.5 else "--")))
    if paths:
        for n,p in enumerate(PATHS):
            for i in range(3): L(p[i],p[i+1],HC[PH[n]],5,1,"-",0.3)
        L(P,H9,"#d08a0c",3,1,"-",0.3); L(Q,H10,"#d08a0c",3,1,"-",0.3)
        for hp in (H9,H10): it.append((pr(hp)[2]+0.5,"p",(pr(hp),"#d08a0c",80)))
    for c in CEN:
        big=paths and c in (P,Q)
        it.append((pr(c)[2]+0.5,"p",(pr(c),"#c0392b" if big else "#1d2422",120 if big else 60)))
    for z,k,d in sorted(it,key=lambda t:t[0]):
        if k=="l":
            a,b,c,w,al,ls=d; ax.plot([a[0],b[0]],[a[1],b[1]],color=c,lw=w,alpha=al,ls=ls,solid_capstyle="round")
        else:
            p,c,s=d; ax.scatter([p[0]],[p[1]],s=s,color=c,zorder=10,edgecolors="white",linewidths=1)
    for txt,p,dx,dy in labels:
        q=pr(p); ax.text(q[0]+dx,q[1]+dy,txt,fontsize=11,fontweight="bold",zorder=12)
    ax.set_xlim(-12,12); ax.set_ylim(-12,12)
    ax.set_title(title,fontsize=13.5,loc="left")
    f.text(0.04,0.03,note,fontsize=9.3,color="#3d4744")
    f.savefig(fn,bbox_inches="tight",facecolor="white"); plt.close(f)
side1=(1.15,-1,-0.15)
fig("cp_1_true.png","1. 本当の姿：二層のシェルピンスキーと 4 本の輪",
    "鏡像の大きな正四面体 2 つ（青・橙、辺² 288）の角に小さな正四面体 8 個（辺² 72）。\n中心は辺 6 の立方体の角。対角線 4 本ごとに輪が 1 本ずつ（4 色）。\nどの辺も 2 本の輪に入るので、輪は少しずらして描いた。",
    (1,0.35,0.6),(0,0,1),None,[0.9]*4,spread=True)
fig("cp_2_top.png","2. 対角線 (1,1,1) から写す：ベンゼンの六角形",
    "この向きから写すと、紫の輪だけが正六角形に見える（影の辺² 24）。\nふた 2 つは真ん中に重なる。教科書が描く平らな六角形はこの影。",
    (1,1,1),(0,0,1),None,[1,0,0,0])
fig("cp_3_side.png","3. 同じ輪を横から：への字 3 つ",
    "紫の輪を (1,1,1) が縦になる向きで横から見る。高さは +3・−3 と交互で、\nへの字 3 つの輪。向かい合う 2 つ（P・Q）は、への字の頂点と底。",
    side1,(1,1,1),None,[1,0,0,0],labels=(("P",P,0.6,0.4),("Q",Q,0.6,-1.2)))
fig("cp_4_next.png","4. 対角線 (−1,1,−1) から写す：別の輪が六角形に",
    "向きを替えると紫の輪は崩れ（点線）、緑の輪が正六角形として現れる。\n同じ点と同じ辺を、写す向きを選ぶだけで別の六角形が見える。",
    (-1,1,-1),(0,0,1),None,[0.45,0,1,0])
fig("cp_5_fork.png","5. P–Q から写す：3 又＋1 個",
    "P を C9、Q を C10 とする。2 つを両方含む輪 3 本（紫・緑・黄）から、\nC9→C10 の半分を 1 本ずつ取る。P–Q の向きから写すと 120 度の 3 又に、H が 1 つ。",
    (-1,-1,1),(0,0,1),None,[0,0,0,0],paths=True,labels=(("C9・C10",P,0.7,0.5),))
fig("cp_6_trip.png","6. 横から：トリプチセンの芯",
    "P–Q を縦にして横から見る。C9（への字の頂点）から C10（底）まで、\n辺 3 本の道が 3 本。どの結合も隣と辺を分け合い、8 個の正四面体を使い切る。",
    (1.15,-1,0.15),(-1,-1,1),None,[0,0,0,0],paths=True,labels=(("C9",P,0.7,-0.3),("C10",Q,0.7,0.3)))
print("ok")
