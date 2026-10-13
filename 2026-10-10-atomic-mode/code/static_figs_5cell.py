# 記事「水の角は4次元の影」の静止画4枚
import math, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
for f in fm.findSystemFonts():
    if "NotoSansCJK-Regular" in f: fm.fontManager.addfont(f)
plt.rcParams["font.family"]=["Noto Sans CJK JP","DejaVu Sans"]
INK="#1b1f27"; MUT="#5d6472"; BLUE="#3b55c9"; RED="#b2402f"; GOLD="#c99a1e"; GR="#9aa3b2"
OUT="/home/claude/mol/v285/images/five-cell/"
ang=lambda c: math.degrees(math.acos((2*c-5)/(2*c+10)))

# 図1 単体の列（中心から頂点への手）
fig,axs=plt.subplots(1,4,figsize=(12,3.4))
def arms2(ax,pts,title,sub):
    ax.set_aspect("equal"); ax.axis("off")
    for i in range(len(pts)):
        for j in range(i+1,len(pts)): ax.plot(*zip(pts[i],pts[j]),color=GR,lw=1,ls="--")
    for p in pts: ax.plot([0,p[0]],[0,p[1]],color=INK,lw=2.2); ax.plot(*p,"o",ms=11,mfc="#e7eef3",mec="#3a5866")
    ax.plot(0,0,"o",ms=12,color=BLUE)
    ax.set_title(title,fontsize=13,color=INK); ax.text(0.5,-0.08,sub,transform=ax.transAxes,ha="center",fontsize=11,color=MUT)
    ax.set_xlim(-1.4,1.4); ax.set_ylim(-1.3,1.4)
arms2(axs[0],[(-1,0),(1,0)],"n＝1 線分","(−1, 1)  180°  BeH₂")
arms2(axs[1],[(math.cos(a),math.sin(a)) for a in (math.pi/2,math.pi/2+2.0944,math.pi/2+4.1888)],"n＝2 正三角形","(−1, 2)  120°  BH₃")
# 正四面体の斜め投影
T=np.array([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]],float)/math.sqrt(3)
R=lambda y,p: np.array([[math.cos(y),0,math.sin(y)],[0,1,0],[-math.sin(y),0,math.cos(y)]])@np.array([[1,0,0],[0,math.cos(p),-math.sin(p)],[0,math.sin(p),math.cos(p)]])
P3=(T@R(0.5,0.35).T)[:,:2]*1.1
arms2(axs[2],[tuple(p) for p in P3],"n＝3 正四面体","(−1, 3)  109.47°  CH₄")
# 正5胞体：4次元座標を2次元へ（対称な図）
P5=[(math.cos(math.pi/2+k*2*math.pi/5),math.sin(math.pi/2+k*2*math.pi/5)) for k in range(5)]
arms2(axs[3],P5,"n＝4 正5胞体","(−1, 4)  104.48°  H₂O")
axs[3].text(0.5,-0.2,"（4次元の形を平面に写した図）",transform=axs[3].transAxes,ha="center",fontsize=9,color=MUT)
fig.tight_layout(rect=(0,0.04,1,0.95)); fig.savefig(OUT+"fig1_simplex_ladder.png",dpi=160); plt.close(fig)

# 図2 位相 c と影の角
fig,ax=plt.subplots(figsize=(8,5.2))
cs=np.linspace(0,1,400); ax.plot(cs,[ang(c) for c in cs],color=BLUE,lw=2.6,label="影の角（正5胞体を面Pで回す）")
pts=[("H₂O",1,104.4776,"c＝1"),("NH₃",0.82066,106.769,"c≈0.821"),("CH₄",0.625,109.4712,"c＝5/8"),("H₃O⁺",0.45932,111.95,"c≈0.459"),("平ら",0,120,"c＝0")]
for n,c,a,lab in pts:
    ax.plot(c,a,"o",ms=9,color=RED,zorder=3); ax.annotate(f"{n}  {a}°\n{lab}",(c,a),xytext=(8,6),textcoords="offset points",fontsize=10,color=INK)
ax.invert_xaxis(); ax.set_xlabel("c ＝ cos²θ（θ は4次元の面Pでの位相）"); ax.set_ylabel("H どうしの角（度）")
ax.set_ylim(103,122); ax.grid(color="#e3e6eb"); ax.legend(loc="upper left",frameon=False)
for s in ("top","right"): ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig(OUT+"fig2_phase_curve.png",dpi=160); plt.close(fig)

# 図3 三つの影（c＝1, 5/8, 0）
S30,S50,QR=math.sqrt(30),math.sqrt(50),math.sqrt(50/3)
H=[[-10/S30,0,QR*math.cos(2*math.pi*i/3),QR*math.sin(2*math.pi*i/3)] for i in range(3)]
L=[15/S30,25/S50,0,0]; F=[15/S30,-25/S50,0,0]
def shadow(v,th): c,s=math.cos(th),math.sin(th); return np.array([v[2],v[0]*c-v[1]*s,v[3]])
fig,axs=plt.subplots(1,3,figsize=(12,4.4))
for ax,(c,title) in zip(axs,[(1,"c＝1：本当の角 104.48°（水）\nL と F の影が重なる"),(0.625,"c＝5/8：正四面体 109.47°（CH₄）"),(0,"c＝0：平ら 120°")]):
    th=math.acos(math.sqrt(c)); Rm=R(0.55,0.3)
    pts3=[shadow(v,th)@Rm.T for v in H+[L,F]]
    ax.set_aspect("equal"); ax.axis("off")
    for i in range(5):
        for j in range(i+1,5): ax.plot([pts3[i][0],pts3[j][0]],[pts3[i][1],pts3[j][1]],color=GR,lw=0.9,ls="--")
    for k,p in enumerate(pts3):
        if k<3: ax.plot([0,p[0]],[0,p[1]],color=INK,lw=2.2); ax.plot(p[0],p[1],"o",ms=13,mfc="#e7eef3",mec="#3a5866")
        else:
            col=GOLD if k==3 else RED; ax.plot([0,p[0]],[0,p[1]],color=col,lw=1.6,ls=":"); ax.plot(p[0],p[1],"o",ms=9,color=col)
    ax.plot(0,0,"o",ms=12,color=BLUE); ax.set_title(title,fontsize=12,color=INK)
    ax.set_xlim(-5,5); ax.set_ylim(-5,5)
fig.text(0.5,0.02,"H（白）・埋まった席 L（金）・4次元目の支点 F（赤）。形は同じ正5胞体で、4次元目を落とす向きだけが違う",ha="center",fontsize=10,color=MUT)
fig.tight_layout(rect=(0,0.05,1,0.94)); fig.savefig(OUT+"fig3_three_shadows.png",dpi=160); plt.close(fig)

# 図4 イオンの向き
fig,ax=plt.subplots(figsize=(8,4.4))
rows=[("H₂O → H₂F⁺（L殻）",104.4776,112.30),("H₂S → H₂Cl⁺（M殻）",92.11,94.24),("NH₃ → H₃O⁺（L殻）",106.769,111.95)]
for i,(lab,a,b) in enumerate(rows):
    y=len(rows)-1-i
    ax.annotate("",xy=(b,y),xytext=(a,y),arrowprops=dict(arrowstyle="->",color=BLUE,lw=2.4))
    ax.plot(a,y,"o",color=INK,ms=8); ax.plot(b,y,"o",color=RED,ms=8)
    ax.text(a-0.4,y+0.18,f"{a}°",ha="center",fontsize=10,color=INK); ax.text(b+0.4,y+0.18,f"{b}°",ha="center",fontsize=10,color=RED)
ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows][::-1]); ax.set_xlim(89,115); ax.set_ylim(-0.6,2.6)
ax.set_xlabel("H どうしの角（度）　黒：中性　赤：陽イオン"); ax.grid(axis="x",color="#e3e6eb")
for s in ("top","right"): ax.spines[s].set_visible(False)
fig.tight_layout(); fig.savefig(OUT+"fig4_cation_opens.png",dpi=160); plt.close(fig)
print("ok")
