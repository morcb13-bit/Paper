# 静止図：3 又＋1 個（C9 の正四面体と、辺を分け合う隣 3 個、H9）
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from itertools import product, combinations
for f in fm.findSystemFonts():
    if "NotoSansCJK" in f.replace(" ","") and "Regular" in f: fm.fontManager.addfont(f)
plt.rcParams["font.family"]=["Noto Sans CJK JP","Noto Sans CJK TC","sans-serif"]
SG=list(product((1,-1),repeat=3))
T={s:[tuple(6*s[k] if i==k else 0 for i in range(3)) for k in range(3)]+[tuple(6*x for x in s)] for s in SG}
C={s:tuple(3*x for x in s) for s in SG}
P=(1,1,-1); nb=[(1,1,1),(1,-1,-1),(-1,1,-1)]; H9=(6,6,-6)
def basis(view,up):
    v=np.array(view,float); v/=np.linalg.norm(v); u=np.array(up,float); u-=u.dot(v)*v; u/=np.linalg.norm(u); x=np.cross(u,v); return x,u,v
def draw(view,up,title,fn,note):
    x,y,v=basis(view,up); pr=lambda p:(np.dot(p,x),np.dot(p,y),np.dot(p,v))
    fig,ax=plt.subplots(figsize=(6,6.6),dpi=150); ax.set_aspect("equal"); ax.axis("off")
    items=[]
    def tet(s,col,lw,alpha):
        V=T[s]
        for a,b in combinations(V,2):
            pa,pb=pr(a),pr(b); items.append(((pa[2]+pb[2])/2,"l",(pa,pb,col,lw,alpha,"-")))
    tet(P,"#c0392b",2.2,0.95)
    for s in nb: tet(s,"#2f6f8f",1.2,0.55)
    # 分け合う辺
    for s in nb:
        e=list(set(T[P])&set(T[s])); pa,pb=pr(e[0]),pr(e[1]); items.append(((pa[2]+pb[2])/2+0.1,"l",(pa,pb,"#7a3fd1",5,1,"-")))
    # 結合（中心どうし、C9–H）
    for s in nb:
        pa,pb=pr(C[P]),pr(C[s]); items.append(((pa[2]+pb[2])/2+0.2,"l",(pa,pb,"#1d2422",1.6,0.8,"--")))
    pa,pb=pr(C[P]),pr(H9); items.append(((pa[2]+pb[2])/2+0.2,"l",(pa,pb,"#1d2422",1.6,0.8,"--")))
    for s in nb: p=pr(C[s]); items.append((p[2]+0.5,"p",(p,"#1d2422",90,"a")))
    p=pr(C[P]); q=pr(H9); same=abs(p[0]-q[0])+abs(p[1]-q[1])<0.5
    items.append((p[2]+0.5,"p",(p,"#c0392b",130,"a" if same else "C9")))
    items.append((q[2]+0.5,"p",(q,"#d08a0c",110,"C9・H9（重なる）" if same else "H9")))
    for z,k,d in sorted(items,key=lambda t:t[0]):
        if k=="l":
            pa,pb,col,lw,al,ls=d; ax.plot([pa[0],pb[0]],[pa[1],pb[1]],color=col,lw=lw,alpha=al,ls=ls,solid_capstyle="round")
        else:
            p,col,sz,lab=d; ax.scatter([p[0]],[p[1]],s=sz,color=col,zorder=10,edgecolors="white",linewidths=1.2)
            if lab!="a": ax.text(p[0]+0.5,p[1]+0.5,lab,fontsize=12,fontweight="bold",zorder=11)
    ax.set_title(title,fontsize=14,loc="left")
    fig.text(0.04,0.04,note,fontsize=9.5,color="#5f6b67",wrap=True)
    fig.savefig(fn,bbox_inches="tight",facecolor="white"); plt.close(fig)
draw((1.6,0.4,0.9),(1,1,-1),"3 又＋1 個 ― 本当の姿（斜めから）","fork_1_true.png",
     "赤：C9 の正四面体（中心に C9）。青：辺を分け合う隣 3 個。紫：分け合う辺（結合）。\nC9→隣の長さ² 36、互いの内積 0（直角）。C9→H9 の長さ² 27。座標は整数。")
draw((1,1,-1),(0,0,1),"3 又＋1 個 ― H9 の軸から見た姿","fork_2_axis.png",
     "C9 と H9 を結ぶ軸に沿って見ると、3 本は互いに 120 度に開いて見える。\n立体の上では 3 本は互いに直角（内積 0）で、120 度は見る向きの影。")
print("ok")
