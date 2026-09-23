# 検定BR2 枝分かれ規則で半加算器を振り分ける配置を全探索する
# 送り先（整数）：明るい五角形 q から、BFS 歩数がいちばん近い別の五芒星の輪（同点は全部）
# 案A：参照スリット r を常に送る（位相オフセット 0/780/1560/2340）。強度 |参照+a+b|²。最大の枚数 1→伸びる 2→二葉 3以上→止まる
# 案B：参照なし。強度>0 の方向へ送る＝和の経路、強度=4 の方向へ送る＝桁上がりの経路（閾値読み）
# 事前の基準（走らせる前に固定）
#   OK   一つの五芒星 g が、送り先 h へ XOR、別の送り先 h' へ AND を出す（一星で半加算器）
#   望ましさの順 (1) OK の星がある (2) 片方を780ずらす対照でその星の OK が消える
#               (3) 参照が軸上の五芒星（配置全体が鏡で閉じる） (4) 部品が少ない（案B＜案A）
#   NG   全配置で OK の星が0
import pickle,itertools
from collections import Counter
exec(open('ha1.py').read().split('d=pickle.load')[0])
SL=list(stars)
x0=U.xy(place[0][0])[0]
dist_star={}
for g in SL: dist_star[g]=bfs(stars[g])
def nearest(q,g):
    best=None; out=[]
    for h in SL:
        if h==g: continue
        dd=min(dist_star[h][r] for r in [q]) if False else dist_star[h][q]
        if best is None or dd<best: best,out=dd,[h]
        elif dd==best: out.append(h)
    return out
NEAR={(q,g):nearest(q,g) for g in SL for q in stars[g]}
KEYS=((0,0),(0,1),(1,0),(1,1))
def amp4(d,off):
    ph=(STEP*d+off)%BASE; return AMP[ph//STEP]
def intens(srcs,q):
    re=im=0
    for dmap,off in srcs:
        v=amp4(dmap[q],off); re+=v[0]; im+=v[1]
    return re*re+im*im
def fn_of(s): return tuple(int(k in s) for k in KEYS)
XOR=(0,1,1,0); AND=(0,0,0,1)
def caseA(sa,sb,r,roff,offB=0):
    dA,dB,dR=dist_star[sa],dist_star[sb],dist_star[r]
    good=[]
    for g in SL:
        if g in (sa,sb,r): continue
        tab={}
        for a,b in KEYS:
            srcs=[(dR,roff)]+([(dA,0)] if a else [])+([(dB,offB)] if b else [])
            v=[intens(srcs,q) for q in stars[g]]; m=max(v)
            top=[q for q,x in zip(stars[g],v) if x==m]
            if len(top)>=3: continue
            for q in top:
                for h in NEAR[(q,g)]: tab.setdefault(h,set()).add((a,b))
        fs={h:fn_of(s) for h,s in tab.items()}
        hx=[h for h,f in fs.items() if f==XOR]; hc=[h for h,f in fs.items() if f==AND]
        if hx and hc: good.append((g,hx,hc))
    return good
def caseB(sa,sb,offB=0):
    dA,dB=dist_star[sa],dist_star[sb]; good=[]
    for g in SL:
        if g in (sa,sb): continue
        sumt={};cart={}
        for a,b in KEYS:
            srcs=([(dA,0)] if a else [])+([(dB,offB)] if b else [])
            for q in stars[g]:
                I=intens(srcs,q)
                for h in NEAR[(q,g)]:
                    if I>0: sumt.setdefault(h,set()).add((a,b))
                    if I==4: cart.setdefault(h,set()).add((a,b))
        hx=[h for h,s in sumt.items() if fn_of(s)==XOR]; hc=[h for h,s in cart.items() if fn_of(s)==AND]
        if hx and hc: good.append((g,hx,hc))
    return good
# 鏡の対と軸上の五芒星
pairs=[]
for g in SL:
    x,y=U.xy(g)
    if x<x0-1e-6:
        m=[h for h in SL if abs(U.xy(h)[0]-(2*x0-x))<1e-4 and abs(U.xy(h)[1]-y)<1e-4]
        if m: pairs.append((g,m[0]))
axis=[g for g in SL if abs(U.xy(g)[0]-x0)<1e-4]
print('鏡の対',len(pairs),'軸上の五芒星',len(axis))
resB=[(i,caseB(sa,sb)) for i,(sa,sb) in enumerate(pairs)]
resB_c=[len(caseB(sa,sb,780)) for sa,sb in pairs]
print('案B: OKの星を持つ対',sum(1 for _,g in resB if g),'/',len(pairs),' OKの星の総数',sum(len(g) for _,g in resB),' 780ずらした対照での総数',sum(resB_c))
rows=[]
for i,(sa,sb) in enumerate(pairs):
    for r in SL:
        if r in (sa,sb): continue
        for roff in (0,780,1560,2340):
            g=caseA(sa,sb,r,roff)
            if g:
                gc=caseA(sa,sb,r,roff,780)
                lost=[t for t in g if t[0] not in [u[0] for u in gc]]
                rows.append((len(lost)>0, r in axis, len(g), len(lost), i, r, roff, g, gc))
print('案A: 配置の総数',len(pairs)*34*4,' OKの星がある配置',len(rows),
      ' うち対照で消える',sum(r[0] for r in rows),' うち参照が軸上',sum(r[0] and r[1] for r in rows))
rows.sort(key=lambda t:(-t[0],-t[1],-t[3],-t[2]))
for t in rows[:8]:
    lost,ax,ng,nl,i,r,roff,g,gc=t
    print(f'   対{i:2d} 参照{tuple(round(v,2) for v in U.xy(r))} 軸上={ax} 参照位相{roff:4}  OKの星{ng} 対照で消える{nl}')
pickle.dump(dict(pairs=pairs,axis=axis,rows=rows,resB=resB),open('/home/claude/prop/br2.pkl','wb'))
