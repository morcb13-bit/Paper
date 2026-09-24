# 検定ALT3b（ALT3 の診断のあと、走らせる前に基準を書き直したもの）
#  噛み合いを「次の桁の五芒星の中心 C_(i+1) が、今の桁の五角形 A_i に一致する数」に替える
#  基準  L1' 各段で噛み合いが 0 でない候補がちょうど一つ  L2 それが +φ³（ALT1 と同じ）  L3 9段で中心がずれない
exec(open('/home/claude/alt3.py').read().split("def run3(")[0])
Mc=U.zadd(small[0],small[2]); u=U.ONE; C2=Mc
for i in range(9):
    cells2=[U.zadd(C2,U.zmul(u,U.zsub(U.zadd(q,q),Mc))) for q in base_cells]
    RCi=ring_centers(cells2,u); nxt,_=step([list(r) for r in R14])
    extra=[U.zadd(C2,U.zmul(u,U.zsub(U.zadd(c,c),Mc))) for c in nxt[-1]]
    G=golden_rhombi_on_axis(list(set(RCi)|set(extra)),u,U.xy(C2)[0]/2)
    M_big,big=G[0]; small_i=[U.zadd(C2,U.zmul(u,U.zsub(U.zadd(v,v),Mc))) for v in small]
    Ci=star_centers(cells2,u); Ai=set(cells2); res=[]
    for mm in range(10):
        uu=U.zmul(U.zt(mm),U.zmul(PHI2,U.PHI))
        if set(U.zadd(M_big,U.zmul(uu,U.zsub(v,C2))) for v in small_i)!=set(big): continue
        nC=set(U.zadd(M_big,U.zmul(uu,U.zsub(g,C2))) for g in Ci)
        res.append(((mm,3),len(Ai&nC),uu))
    nz=[r for r in res if r[1]>0]
    print(f'段{i}  軸上の黄金のひし形 {len(G)}  候補と噛み合い {[(r[0],r[1]) for r in res]}  選ばれた {nz[0][0] if len(nz)==1 else "決まらない"}  中心が同じ {M_big==C2}')
    if len(nz)!=1: break
    u=U.zmul(u,nz[0][2]); C2=M_big
