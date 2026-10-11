# 記事用の動かない説明図（座標は整数で作り、描画だけ浮動小数）
import itertools, math
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Circle
plt.rcParams['font.family']='Noto Sans CJK JP'
INK='#1b2333'; MUT='#5d6779'; LINE='#cfd5df'; A='#c2410c'; B='#0f766e'; D='#7c3aed'; F='#1d4ed8'; ACC='#3b5bdb'
ROW=[A,B,D,F]
OUT='/home/claude/mol/images/'

def view(p,yaw=0.6,pitch=0.35):
    x,y,z=p; cy,sy=math.cos(yaw),math.sin(yaw)
    X=x*cy-y*sy; Y0=x*sy+y*cy; cp,sp=math.cos(pitch),math.sin(pitch)
    return (X, z*cp-Y0*sp, z*sp+Y0*cp)
def axview(p,pitch=0.75):  # (1,1,1) を縦に
    x,y,z=p; a=(x-y)/math.sqrt(2); b=(x+y-2*z)/math.sqrt(6); c=(x+y+z)/math.sqrt(3)
    return view((a,b,c),0.35,pitch)
def tri2d(ax,pts,fc,ec,alpha=0.35,lw=1.6,ls='-'):
    ax.add_patch(Polygon(pts,closed=True,fc=fc,ec=ec,alpha=alpha,lw=0))
    ax.add_patch(Polygon(pts,closed=True,fill=False,ec=ec,lw=lw,ls=ls))
def setup(ax,title=None):
    ax.set_aspect('equal'); ax.axis('off')
    if title: ax.set_title(title,fontsize=15,color=INK,pad=10)

# 図1 水素：平らな三角形 → 正四面体（中心固定、m=1 の整数）
m=1
def hyd(k):
    base=[(4*m-k,-2*m-k,-2*m-k),(-2*m-k,4*m-k,-2*m-k),(-2*m-k,-2*m-k,4*m-k)]; ap=(3*k,3*k,3*k); return base,ap
fig,axs=plt.subplots(1,3,figsize=(15,5.4))
for ax,k,t in zip(axs,[0,0.5,1],['開く前：平らな正三角形（2次元）','途中','開ききる：正四面体（3次元）']):
    base,ap=hyd(k); P=[axview(p) for p in base]; Q=axview(ap)
    faces=[[P[0],P[1],P[2]],[P[0],P[1],Q],[P[1],P[2],Q],[P[2],P[0],Q]]
    for f in sorted(faces,key=lambda f:sum(p[2] for p in f)):
        tri2d(ax,[(p[0],p[1]) for p in f],A,A,0.12,1.4)
    for p in P: ax.plot(p[0],p[1],'o',color=A,ms=8)
    ax.plot(Q[0],Q[1],'o',color=INK,ms=9)
    o=axview((0,0,0)); ax.plot(o[0],o[1],'+',color=ACC,ms=14,mew=2)
    ax.set_xlim(-6,6); ax.set_ylim(-5,6.5); setup(ax,t)
fig.text(0.5,0.02,'底が k 下がる間に頂点が 3k 上がる。四点の和は常に (0,0,0)＝中心（青い十字）は動かない',ha='center',color=MUT,fontsize=12)
fig.savefig(OUT+'fig1_hydrogen.png',dpi=110,bbox_inches='tight',facecolor='white'); plt.close(fig)

# 図2 ヘリウム：六芒星（閉じた）と星型八面体（開いた、立方体の面に×）
fig,axs=plt.subplots(1,2,figsize=(12,6))
ax=axs[0]
for s,c,ls in[(1,A,'-'),(-1,B,'--')]:
    pts=[(s*math.cos(math.pi/2+2*math.pi*j/3),s*math.sin(math.pi/2+2*math.pi*j/3)) for j in range(3)]
    tri2d(ax,pts,c,c,0.2,2.2,ls)
ax.plot(0,0,'+',color=ACC,ms=14,mew=2); ax.set_xlim(-1.3,1.3); ax.set_ylim(-1.3,1.3)
setup(ax,'閉じたとき：一つの面に互い違いの二枚（六芒星）')
ax=axs[1]
V=[p for p in itertools.product((-1,1),repeat=3)]
for p,q in itertools.combinations(V,2):
    if sum(a!=b for a,b in zip(p,q))==1:
        P,Q=view(p),view(q); ax.plot([P[0],Q[0]],[P[1],Q[1]],color=LINE,lw=1.5,ls=':')
TA=[p for p in V if p[0]*p[1]*p[2]>0]; TB=[p for p in V if p[0]*p[1]*p[2]<0]
for T,c,ls in[(TA,A,'-'),(TB,B,'--')]:
    for p,q in itertools.combinations(T,2):
        P,Q=view(p),view(q); ax.plot([P[0],Q[0]],[P[1],Q[1]],color=c,lw=2.4,ls=ls)
    for p in T: P=view(p); ax.plot(P[0],P[1],'o',color=c,ms=8)
for p in[(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]:
    P=view(p); ax.plot(P[0],P[1],'o',color=ACC,ms=6)
ax.set_xlim(-1.9,1.9); ax.set_ylim(-1.9,1.9)
setup(ax,'開いたとき：星型八面体（立方体の各面に×）')
fig.text(0.5,0.02,'橙＝水素の正四面体、緑＝その点対称の鏡像。青い点は面の中心＝×の交点（芯の正八面体の頂点）',ha='center',color=MUT,fontsize=12)
fig.savefig(OUT+'fig2_helium.png',dpi=110,bbox_inches='tight',facecolor='white'); plt.close(fig)

# 図3 一つの面の網：一辺 n=4、段 1・3・5・7
fig,ax=plt.subplots(figsize=(8,7.4))
n=4; h=math.sqrt(3)/2
def P2(i,j): return (j-i/2, -i*h)   # i：段、j：段の中の位置
for i in range(n):
    for j in range(i+1):
        tri2d(ax,[P2(i,j),P2(i+1,j),P2(i+1,j+1)],ROW[i],ROW[i],0.55,1.2)
    for j in range(i):
        tri2d(ax,[P2(i,j),P2(i,j+1),P2(i+1,j+1)],ROW[i],ROW[i],0.22,1.2)
    ax.text(n/2+0.55-(n-1-i)*0.5+ (n-1-i)*0.5, -(i+0.5)*h, f'{2*i+1} 個', va='center',fontsize=13,color=ROW[i])
ax.set_xlim(-2.4,3.4); ax.set_ylim(-n*h-0.3,0.4)
setup(ax,'正四面体の一つの面を一辺4で割った網')
fig.text(0.5,0.03,'角から段ごとに 1・3・5・7 個。濃い＝上向き、淡い＝下向き。合わせて 16＝4²',ha='center',color=MUT,fontsize=12)
fig.savefig(OUT+'fig3_face_rows.png',dpi=110,bbox_inches='tight',facecolor='white'); plt.close(fig)

# 図4 球に閉じた網 n=1..4（中心の光源から投影）
TAi=[(1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1)]; FA=[(0,1,2),(0,1,3),(0,2,3),(1,2,3)]
def net(n):
    pts={}; tris=[]
    for f in FA:
        A_,B_,C_=(TAi[i] for i in f)
        def pt(a,b):
            c=n-a-b; return tuple(a*A_[k]+b*B_[k]+c*C_[k] for k in range(3))
        for a in range(n):
            for b in range(n-a):
                tris.append((pt(a,b),pt(a+1,b),pt(a,b+1)))
                if a+b<=n-2: tris.append((pt(a+1,b),pt(a,b+1),pt(a+1,b+1)))
    P=set(p for t in tris for p in t); E=set()
    for t in tris:
        for u,v in((t[0],t[1]),(t[1],t[2]),(t[2],t[0])): E.add(tuple(sorted((u,v))))
    return tris,P,E
fig,axs=plt.subplots(1,4,figsize=(18,5.2))
for ax,n in zip(axs,[1,2,3,4]):
    tris,P,E=net(n)
    def sph(p): r=math.sqrt(sum(x*x for x in p)); return view(tuple(x/r for x in p))
    for t in sorted(tris,key=lambda t:sum(sph(p)[2] for p in t)):
        Q=[sph(p) for p in t]; front=sum(q[2] for q in Q)>0
        tri2d(ax,[(q[0],q[1]) for q in Q],A,A,0.18 if front else 0.05,1.2 if front else 0.4)
    ax.set_xlim(-1.15,1.15); ax.set_ylim(-1.15,1.15)
    setup(ax,f'n＝{n}　点 {len(P)}＝{len(P)-2}＋2')
    ax.text(0,-1.32,f'面 {len(tris)}・辺 {len(E)}',ha='center',color=MUT,fontsize=12)
fig.text(0.5,-0.03,'正四面体の面を一辺 n で割り、中心の光源から球へ投影。点＝辺−面＋2＝2n²＋2、点−2 が 2・8・18・32',ha='center',color=MUT,fontsize=12)
fig.savefig(OUT+'fig4_sphere_nets.png',dpi=110,bbox_inches='tight',facecolor='white'); plt.close(fig)

# 図5 殻の大きさ：半径 ∝ n と ∝ n²（断面の同心円）
fig,axs=plt.subplots(1,2,figsize=(13,6.6))
for ax,pw,t in zip(axs,[1,2],['目盛りの長さをそろえる：半径 1：2：3：4','目盛りの長さも n 倍：半径 1：4：9：16']):
    R=max(k**pw for k in range(1,5))
    for k,c in zip(range(1,5),ROW):
        r=k**pw/R; ax.add_patch(Circle((0,0),r,fill=False,ec=c,lw=2.4))
        ax.text(r*0.72+0.02,r*0.72+0.02,'KLMN'[k-1],color=c,fontsize=14,fontweight='bold')
    ax.plot(0,0,'+',color=ACC,ms=12,mew=2); ax.set_xlim(-1.1,1.1); ax.set_ylim(-1.1,1.1); setup(ax,t)
fig.text(0.5,0.02,'どちらも殻 n の点の数は 2n²＋2。右は水素型の殻の半径（n² に比例）と同じ比',ha='center',color=MUT,fontsize=12)
fig.savefig(OUT+'fig5_scale.png',dpi=110,bbox_inches='tight',facecolor='white'); plt.close(fig)
print('ok')
