# 検定ALT3 三階建て：A_i（五角形）→ B_i（円環の中心＋成長の一行）→ C_i（五芒星の中心）→ A_{i+1}
#  ±φ³ の二つの候補のうち、C_i（A_i の五芒星の中心）と整数で噛み合う方を選ぶ。内積の判定は使わない
#  噛み合い＝C_i の点が、候補 A_{i+1} の五角形の中心と一致する数
# 事前の基準
#   L1 各段で、噛み合いが 0 でない候補がちょうど一つ
#   L2 選ばれた候補が ALT1（内積で選んだ +φ³）と同じ
#   L3 9段続けて中心と比がずれない
#   対照 C を使わない → 候補は二つのまま（決まらない）
import sys
exec(open('/home/claude/alt1.py').read().split("def run(grow=True")[0])
def star_centers(cells2,u):
    S=set(cells2); arm=[U.zmul(u,U.zmul((2,0,0,0),U.zmul(U.PHI,U.zt(k)))) for k in range(10)]
    cand={U.zsub(q,a) for q in cells2 for a in arm}
    out=[]
    for g in cand:
        for kk in (0,1):
            if all(U.zadd(g,arm[(kk+2*j)%10]) in S for j in range(5)): out.append(g); break
    return out
def run3(steps=9):
    Mc=U.zadd(small[0],small[2]); u=U.ONE; C2=Mc; log=[]
    for i in range(steps):
        cells2=[U.zadd(C2,U.zmul(u,U.zsub(U.zadd(q,q),Mc))) for q in base_cells]
        RCi=ring_centers(cells2,u); nxt,_=step([list(r) for r in R14])
        extra=[U.zadd(C2,U.zmul(u,U.zsub(U.zadd(c,c),Mc))) for c in nxt[-1]]
        B=list(set(RCi)|set(extra)); G=golden_rhombi_on_axis(B,u,U.xy(C2)[0]/2)
        if len(G)!=1: log.append((i,'B で道が一つに決まらない',len(G))); break
        M_big,big=G[0]
        small_i=[U.zadd(C2,U.zmul(u,U.zsub(U.zadd(v,v),Mc))) for v in small]
        cands=[]
        for mm in range(10):
            uu=U.zmul(U.zt(mm),U.zmul(PHI2,U.PHI))
            if set(U.zadd(M_big,U.zmul(uu,U.zsub(v,C2))) for v in small_i)==set(big): cands.append(((mm,3),uu))
        Ci=set(star_centers(cells2,u))
        mesh=[]
        for mn,uu in cands:
            nextA=set(U.zadd(M_big,U.zmul(uu,U.zsub(q,C2))) for q in cells2)
            mesh.append((mn,len(Ci&nextA)))
        nz=[(mn,c) for mn,c in mesh if c>0]
        log.append((i,len(Ci),mesh,M_big==C2))
        if len(nz)!=1: break
        mn=nz[0][0]; uu=[x for m_,x in cands if m_==mn][0]
        u=U.zmul(u,uu); C2=M_big
    return log
for t in run3():
    print('   段',t[0],' C の五芒星',t[1],' 候補ごとの噛み合い',t[2] if len(t)>2 else '',' 中心が同じ',t[3] if len(t)>3 else '')
print('--- 診断：段0で噛み合った点と、五芒星の層どうし・円環の層との噛み合い')
Mc=U.zadd(small[0],small[2]); u=U.ONE; C2=Mc
cells2=[U.zadd(q,q) for q in base_cells]
RCi=ring_centers(cells2,u); nxt,_=step([list(r) for r in R14]); extra=[U.zadd(c,c) for c in nxt[-1]]
B=list(set(RCi)|set(extra)); M_big,big=golden_rhombi_on_axis(B,u,U.xy(C2)[0]/2)[0]
small_i=[U.zadd(v,v) for v in small]; Ci=set(star_centers(cells2,u))
for mm in (0,5):
    uu=U.zmul(U.zt(mm),U.zmul(PHI2,U.PHI))
    nA=set(U.zadd(M_big,U.zmul(uu,U.zsub(q,C2))) for q in cells2)
    nC=set(U.zadd(M_big,U.zmul(uu,U.zsub(g,C2))) for g in Ci)
    nB=set(U.zadd(M_big,U.zmul(uu,U.zsub(g,C2))) for g in RCi)
    pts=sorted((round(U.xy(p)[0]/2,2),round(U.xy(p)[1]/2,2)) for p in Ci&nA)
    print(f'  候補 m={mm}: C_i∩A_(i+1)の五角形 {pts}')
    print(f'           B_i∩A_(i+1)の五角形 {len(set(B)&nA)}  C_i∩B_(i+1)の円環 {len(Ci&nB)}  B_i∩C_(i+1) {len(set(B)&nC)}  A_iの五角形∩C_(i+1) {len(set(cells2)&nC)}')
