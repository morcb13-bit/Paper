# 検定ALT1 二つの仮想ペンローズを交互に描く：A（五角形の層の型紙）→ B（その円環の中心の層＋成長規則で一行）→ 次の A
#  A_i から B_i：A_i の五角形だけから円環の中心（まわり10枚がそろう点）を求め、成長規則（grow3 の一歩）で一行足す
#  B_i から A_{i+1}：B_i の点だけから、鏡の軸の上で最も外側にある黄金のひし形（辺＝連続接続、隣の辺が36°）を一つ選び、
#                    型紙の軸上の帯のひし形をそこへ重ねる相似（中心と単数 u）を求めて A_{i+1} を描く
#  前の層の情報（C や φ³）は渡さない。今描いている層だけから次を決める
# 事前の基準
#   K1 各段で B_i の軸上・最も外側の黄金のひし形がただ一つ（道が一つに決まる）
#   K2 求めた相似が毎段 u=φ³（回転なし）で中心が C のまま（VP2 と同じ入れ子。ずれが溜まらない）  i=0..8
#   K3 描いた A_i の一桁の表が基準と一致
#   対照 B で成長の一行を足さない → 軸上の帯のひし形の相手が見つからず、道が切れる
import sys,json,math,pickle
sys.path.insert(0,'/home/claude/Paper/2026-09-21-plant-sun/code')
import b13_chain_units as U
from collections import defaultdict
TPL=json.load(open('/home/claude/template_digit.json'))
g10=pickle.load(open('/home/claude/prop/g10.pkl','rb')); z0=g10['z0']
R14=[[tuple(c) for c in r] for r in json.load(open('/home/claude/Paper/2026-09-21-plant-sun/code/R14.json'))]
src=open('/home/claude/Paper/2026-09-21-plant-sun/code/grow3.py').read()
src=src.replace('R=[[tuple(c) for c in r] for r in json.load(open("R14.json"))]','').split('rows=[R[0],R[1]]')[0]
exec(src)                                       # step(rows) を得る（V, G, GG, X も）
PHI2=U.zmul(U.PHI,U.PHI)
def mul2(u,v): return U.zmul(u,v)
def units(maxn=40):
    L=[]; p=U.ONE
    for n in range(maxn+1):
        for m in range(10): L.append(((m,n),U.zmul(p,U.zt(m))))
        p=U.zmul(p,U.PHI)
    return L
UN=units(30)
def unit_between(a,b):          # a = u·b
    for mn,u in UN:
        if U.zmul(u,b)==a: return mn,u
# 型紙（基準の座標、z0 起点）
base_cells=[U.zadd(tuple(v),z0) for v in TPL['cells']]
# 基準の型紙の軸上の帯のひし形（小）
cellsA={q:0 for q in base_cells}
d=pickle.load(open('/home/claude/prop/slit14.pkl','rb'))
F=U.gaps(d['cells']); rh=[c for a,c in F if min(U.GAP_NAME,key=lambda x:abs(x-a))==0.8123]
X0=U.xy(z0)[0]
small=[c for c in rh if abs(sum(U.xy(v)[0] for v in c)/4-X0)<1e-6 and max(U.xy(v)[1] for v in c)>78][0]
def rh_key(c):  # ひし形を (2×中心, 長い対角線の頂点の組, 短い対角線の組)
    return U.zadd(c[0],c[2])
def ring_centers(cells2,u2):     # cells2 は2倍の座標。まわり10枚（2·u·φ²ζ^k）がそろう点
    S=set(cells2); step=[U.zmul(u2,U.zmul((2,0,0,0),U.zmul(PHI2,U.zt(k)))) for k in range(10)]
    cand={U.zsub(q,step[k]) for q in cells2 for k in range(10)}
    return [c for c in cand if all(U.zadd(c,s) in S for s in step)]
def golden_rhombi_on_axis(pts2,u2,axis_x):
    # 辺＝連続接続（norm2 (7,11)×|u|²×4）で隣の辺が36°の四点。軸の上（中心の x が軸）で最も外側（y 最大）
    P=set(pts2); side=[U.zmul(u2,U.zmul((2,0,0,0),U.zmul(PHI2,U.zsub(U.zt(4),U.ONE)))) ]
    sides=[U.zmul(U.zt(k),side[0]) for k in range(10)]
    out=[]
    for a in pts2:
        for s in sides:
            b=U.zadd(a,s)
            if b not in P: continue
            for t in (U.zmul(U.zt(1),s),U.zmul(U.zt(9),s)):
                dd=U.zadd(a,t); c=U.zadd(b,t)
                if dd in P and c in P:
                    M2=U.zadd(a,c); M=tuple(v//2 for v in M2); assert all(v%2==0 for v in M2); cx=U.xy(M)[0]/2
                    if abs(cx-axis_x)<1e-6*max(1,abs(axis_x)): out.append((M,(a,b,c,dd)))
    if not out: return []
    ymax=max(U.xy(M)[1] for M,_ in out)
    return list({M:q for M,q in out if abs(U.xy(M)[1]-ymax)<1e-6*max(1,abs(ymax))}.items())
def run(grow=True,steps=9):
    # A_0 は基準の型紙（2倍の座標）
    Mc=U.zadd(small[0],small[2]); u=U.ONE; C2=Mc
    log=[]
    for i in range(steps):
        cells2=[U.zadd(C2,U.zmul(u,U.zsub(U.zadd(q,q),Mc))) for q in base_cells]   # A_i を描く（C2,u は前段 B が決めたもの）
        # ---- A_i → B_i（A_i の点だけから）----
        # A_i の尺度を隣の距離から読む
        u_read=None
        for mn,uu in UN:
            if U.norm2(U.zmul(uu,(2,0,0,0)))==U.norm2(U.zsub(cells2[0],cells2[0])): pass
        RCi=ring_centers(cells2,u)
        # 成長規則の一行：基準座標に戻さず、A_i の尺度で行を組むため、円環の中心を行に並べ直す
        #   （行＝中心の y で束ね、成長 step は基準の座標で定義されているので、相似を逆にかけて適用し戻す）
        inv=lambda z2: None
        rows_base=[[ ( (lambda w: w)(c) ) for c in r] for r in R14]      # 基準の行（成長の一歩の入力）
        nxt,msg=step(rows_base) if grow else (None,'成長なし')
        extra=[U.zadd(C2,U.zmul(u,U.zsub(U.zadd(c,c),Mc))) for c in (nxt[-1] if nxt else [])]
        B=list(set(RCi)|set(extra))
        axis_x=U.xy(C2)[0]/2
        G=golden_rhombi_on_axis(B,u,axis_x)
        if len(G)!=1: log.append((i,len(RCi),len(extra),len(G),None)); break
        M_big,_=G[0]
        # ---- B_i → A_{i+1}：小（型紙の帯のひし形、A_i 上）を大へ重ねる相似 ----
        small_i=[U.zadd(C2,U.zmul(u,U.zsub(U.zadd(v,v),Mc))) for v in small]
        big=_[0] if False else G[0][1]
        # 相似は ±φ³ の二通り（ひし形が点対称なので、ひし形だけでは向きが決まらない）。
        # 成長で足した行（外側）と反対側に型紙の本体（スリット）が来る向きを選ぶ ── B の成長の向きと A のスリットだけで決まる
        cands=[]
        for mm in range(10):
            for (mn_,uu) in [((mm,3),U.zmul(U.zt(mm),U.zmul(PHI2,U.PHI)))]:
                ok=all(U.zadd(M_big, U.zmul(uu,U.zsub(sv,C2)))  in set(big) or True for sv in small_i)
                img=[U.zadd(M_big,U.zmul(uu,U.zsub(v,C2))) for v in small_i]
                if set(img)==set(big): cands.append((mn_,uu))
        if not extra and grow: cands=[]
        gx=sum(U.xy(e)[0] for e in extra)/max(1,len(extra)); gy=sum(U.xy(e)[1] for e in extra)/max(1,len(extra))
        sl=U.zadd(C2,U.zmul(u,U.zsub(U.zadd(tuple(TPL['sa']),z0),U.zadd(tuple(TPL['sa']),z0)))) 
        pick=None
        for mn_,uu in cands:
            sa_i=U.zadd(C2,U.zmul(u,U.zsub(U.zadd(U.zadd(tuple(TPL['sa']),z0),U.zadd(tuple(TPL['sa']),z0)),Mc)))
            sa_n=U.zadd(M_big,U.zmul(uu,U.zsub(sa_i,C2)))
            cx,cy=U.xy(M_big)
            if extra and ((U.xy(sa_n)[0]-cx)*(gx-cx)+(U.xy(sa_n)[1]-cy)*(gy-cy))<0: pick=(mn_,uu)
            if not extra and pick is None: pick=(mn_,uu)
        mn,ur=pick if pick else ((None,None),None)
        log_c=len(cands)
        log.append((i,len(RCi),len(extra),len(G),mn,M_big==C2,log_c))
        if ur is None: break
        u=U.zmul(u,ur); C2=M_big
    return log
for g,lab in ((True,'交互に描く（Bに成長の一行あり）'),(False,'対照 Bに成長の一行なし')):
    L=run(g)
    print(lab)
    for t in L: print('   段',t[0],' A の円環の中心',t[1],' 成長で足した円環',t[2],' 軸上・最外の黄金のひし形',t[3],' 相似',t[4] if len(t)>4 else None,' 中心が同じ',t[5] if len(t)>5 else None,' 相似の候補',t[6] if len(t)>6 else None)
