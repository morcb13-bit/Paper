import json
from judge_d0 import zm, za, zs, sgn
def z_int(T):
    S=sum(T); Q=sum(t*t for t in T); return S, S*S>9*Q
def z_phi(T):
    S=(0,0);Q=(0,0)
    for t in T: t=tuple(t); S=za(S,t); Q=za(Q,zm(t,t))
    return S, sgn(za(zm(S,S),zs(Q,-9)))>0
phi=(1+5**0.5)/2
for m in ['fan','off']:
    r=json.load(open('w1_%s.json'%m))
    SJ,bJ=z_phi(r['J']); SD,bD=z_int(r['D']); SY,bY=z_int(r['Ly']); SZ,bZ=z_int(r['Lz'])
    print(f"{m:4s} 帳簿={r['ok']}  噴流 ΣJ≈{(SJ[0]+SJ[1]*phi)/1e3:+.0f}千 雑音越え={bJ and sgn(SJ)>0}  抗力 ΣD={SD:+d} 雑音越え={bD and SD>0}  横 ΣLy={SY:+d} 越え={bY}  ΣLz={SZ:+d} 越え={bZ}")
