# 一本道：閉じた三つ組の中から毎回一つだけ選んで進む。二歩で一組の剛体運動 M になり、列は M の繰り返し（ねじ運動）。
import sympy as sp
T0=[sp.Matrix(v) for v in [(1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1)]]
def rotP(P,a,b,th):
    ax=(b-a)/sp.sqrt((b-a).dot(b-a)); m=(a+b)/2; c,s=sp.cos(th),sp.sin(th)
    return [sp.simplify(m+(p-m)*c+ax.cross(p-m)*s+ax*(ax.dot(p-m))*(1-c)) for p in P]
def motion(P,Q):   # P→Q の剛体運動（頂点のラベル対応）を 4x4 で
    A=sp.Matrix.hstack(*[P[i]-P[0] for i in (1,2,3)]); B=sp.Matrix.hstack(*[Q[i]-Q[0] for i in (1,2,3)])
    R=sp.simplify(B*A.inv()); t=sp.simplify(Q[0]-R*P[0]); return R,t
for s1,s2,name in ((1,1,"同じ向きに回し続ける"),(1,-1,"向きを交互に変える")):
    T1=rotP(T0,T0[0],T0[1],s1*2*sp.pi/3)
    T2=rotP(T1,T1[2],T1[3],s2*2*sp.pi/3)
    R,t=motion(T0,T2)
    tr=sp.simplify(R.trace()); cth=sp.simplify((tr-1)/2)
    # 回転軸と軸方向の移動
    w,v=None,None
    ev=(R-sp.eye(3)).nullspace(); axv=sp.simplify(ev[0]/sp.sqrt(ev[0].dot(ev[0])))
    shift=sp.simplify(t.dot(axv))
    print(f"\n[{name}] 二歩の運動： det={sp.simplify(R.det())}  回転角 cos={cth}（{sp.N(sp.acos(cth)*180/sp.pi,8)}°）  軸方向の進み²={sp.simplify(shift**2)}")
    # 何回で回転が元に戻るか（回転角が 360° の有理数倍か）
    ang=sp.acos(cth)/(2*sp.pi); print("  回転角 / 360° =",sp.nsimplify(sp.N(ang,30),rational=False), " 数値",sp.N(ang,15))
    # 中心の列
    C=[sp.zeros(3,1)]; P=T0
    for k in range(6):
        P=[sp.simplify(R*p+t) for p in P]; C.append(sp.simplify(sum(P,sp.zeros(3,1))/4))
    print("  二歩ごとの中心の距離²:",[sp.simplify((C[k+1]-C[k]).dot(C[k+1]-C[k])) for k in range(3)])

# 一歩ごとの運動 S：T0（ラベル 0,1,2,3）を T1（次に回す辺が 0,1 に来るよう 2,3,0,1 と読む）へ
T1=rotP(T0,T0[0],T0[1],2*sp.pi/3)
R,t=motion(T0,[T1[2],T1[3],T1[0],T1[1]])
cth=sp.simplify((R.trace()-1)/2); ev=(R-sp.eye(3)).nullspace()
axv=sp.simplify(ev[0]/sp.sqrt(ev[0].dot(ev[0]))) if ev else None
print("\n[一歩] det",sp.simplify(R.det())," 回転角 cos=",cth,f"（{sp.N(sp.acos(cth)*180/sp.pi,8)}°）", " 軸方向の進み²=",sp.simplify(t.dot(axv)**2) if axv is not None else None)
P=T0; C=[sp.zeros(3,1)]
for k in range(4):
    P=[sp.simplify(R*p+t) for p in P]; C.append(sp.simplify(sum(P,sp.zeros(3,1))/4))
print("  一歩ごとの中心の距離²:",[sp.simplify((C[k+1]-C[k]).dot(C[k+1]-C[k])) for k in range(3)],
      "  二歩先との距離²:",[sp.simplify((C[k+2]-C[k]).dot(C[k+2]-C[k])) for k in range(2)])
