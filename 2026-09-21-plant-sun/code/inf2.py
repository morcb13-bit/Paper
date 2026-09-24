# 照合：置き換えで育てた赤10枚の輪の頂点に、10枚の担体の円環の中心が乗るか
# 写し：v → z0 + w·v、w = φ²(ζ⁴−1)ζ^m（辺の長さ＝円環の中心どうしの連続接続 CONT）
import b13_chain_units as U, inf1 as I
exec(open('gen10.py').read().split("pid={}")[0])
RC=set(place(c,k) for (w,r,c,k) in rings); z0=(2,-2,0,-3)
EDGE=U.norm2(U.CONT[0])
def d2(z): x,y=U.xy(U.zsub(z,z0)); return x*x+y*y
T=I.wheel
for n in range(1,10):
    T=I.step(T)
    R=1.618034**n
    for m in (0,1):
        w=U.zmul(U.zmul(U.PHI2,U.zsub(U.zt(4),U.ONE)),U.zt(m))
        f=lambda v:U.zadd(z0,U.zmul(w,v))
        V=set(); E=set()
        for _,A,B,C in T:
            for x,y in ((A,B),(A,C),(B,C)):
                a,b=f(x),f(y); V|={a,b}
                if U.norm2(U.zsub(a,b))==EDGE: E.add(frozenset((a,b)))
        lim=(0.8*R*4.9798)**2               # 輪の内側だけ数える（場所の選別にだけ浮動小数）
        inside=[c for c in RC if d2(c)<lim]
        hit=sum(1 for c in inside if c in V)
        S=set(inside); ce={frozenset((a,U.zadd(a,v))) for a in inside for v in U.CONT if U.zadd(a,v) in S}
        eh=sum(1 for e in ce if e in E)
        print(f"n={n} m={m}  円環の中心 {hit}/{len(inside)}  連続接続の辺 {eh}/{len(ce)}")
