#  final_cages.py ── 檻を行列の冪で確定させる（RT面30 と 笠10）
#      各 m で、多項式から出した候補を T^k=I で確かめ、真の約数まで落とす
import numpy as np, random, math
MARK="\n# "+"="*64+" 実行\n"
exec(open('cage_scan.py').read().split(MARK)[0])
random.seed(13)
nbr_rt,cap=rt_face_graph(); S=set(cap)
nbr_cap={f:[g for g in nbr_rt[f] if g in S] for f in cap}
def mat(rows,n,m):
    A=np.zeros((n,n))
    for i,r in enumerate(rows):
        for j,c in r: A[i,j]=c%m
    return A
def mm(A,B,m): return np.mod(A@B,m)
def mpow(A,e,m):
    R=np.eye(A.shape[0]); A=A.copy()
    while e:
        if e&1: R=mm(R,A,m)
        A=mm(A,A,m); e>>=1
    return R
def isid(A): return np.array_equal(A,np.eye(A.shape[0]))
MS=[5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61]
def facstr(k):
    f=fact(k); return " · ".join(f"{p}^{e}" if e>1 else str(p) for p,e in sorted(f.items()))
for name,nb in (("RT面30",nbr_rt),("笠10",nbr_cap)):
    arcs,rows=build(nb,2); n=len(arcs)
    print(f"\n【{name}】 t=2 有向辺{n}")
    print(f"  {'m':>3}  {'檻':>8}  {'素因数':<22} {'24':>4}  {'一致':>4}")
    for m in MS:
        cand=order_mod(minpoly_mod(rows,m),m)
        A=mat(rows,n,m)
        ok=isid(mpow(A,cand,m))
        k=cand
        if ok:
            for p in sorted(fact(cand)):
                while k%p==0 and isid(mpow(A,k//p,m)): k//=p
        print(f"  {m:>3}  {k:>8}  {facstr(k):<22} {'割る' if k%24==0 else '—':>4}  {'OK' if ok else 'NG':>4}")
