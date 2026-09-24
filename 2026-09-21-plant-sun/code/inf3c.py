# 前提の確認：ひし形のペンローズで「距離＝辺の長さ」の対は必ず辺か
import b13_chain_units as U, inf1 as I
T=I.wheel
for _ in range(8): T=I.step(T)
V=set(); E=set(); UN=[U.zt(k) for k in range(10)]
for _,A,B,C in T:
    V|={A,B,C}
    for x,y in ((A,B),(A,C),(B,C)):
        if U.zsub(y,x) in UN: E.add(frozenset((x,y)))
pairs={frozenset((v,U.zadd(v,u))) for v in V for u in UN if U.zadd(v,u) in V}
print("距離1の対",len(pairs),"そのうち辺",len(pairs&E),"辺でない",len(pairs-E))
