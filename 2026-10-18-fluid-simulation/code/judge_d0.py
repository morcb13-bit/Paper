import json
def zm(x,y): a,b=x;c,d=y; return (a*c+b*d,a*d+b*c+b*d)
def za(x,y): return (x[0]+y[0],x[1]+y[1])
def zs(x,k): return (x[0]*k,x[1]*k)
def sgn(x):
    p,q=x;A=2*p+q;B=q
    if A>=0 and B>=0: return 0 if A==0 and B==0 else 1
    if A<=0 and B<=0: return -1
    return (1 if A*A>5*B*B else -1) if A>0 else (1 if 5*B*B>A*A else -1)
def notbias(T):
    S=(0,0); Q=(0,0)
    for t in T: S=za(S,tuple(t)); Q=za(Q,zm(tuple(t),tuple(t)))
    return sgn(za(zs(Q,9),zs(zm(S,S),-1)))>0, S
import sys
for m in (sys.argv[1:] or ['normal','rev','nocycle']):
    try: r=json.load(open('d0_%s.json'%m))
    except Exception: print(m,'未完'); continue
    Pu,Pm,Pd=[tuple(x) for x in r['Pb']]; nu,nm,nd=r['nv']
    s=sgn(Pu)
    a=zs(Pd,nu*s); b=zs(Pu,nd*s)            # 頂点数で揃えた下流と上流（符号は流れの向きで揃える）
    keep=sgn(za(zs(a,5),zs(b,-4)))>0 and sgn(za(zs(a,5),zs(b,-6)))<0
    bu,Su=notbias(r['Tu']); bw,Sw=notbias(r['Tw'])
    print(f"{m:8s} 保存={'OK' if r['ok'] else 'NG'}  P∥の符号={'+' if s>0 else '-'}  流れが残る(0.8〜1.2倍)={'OK' if keep else 'NG'}  横u偏りなし={'OK' if bu else 'NG'} 横v偏りなし={'OK' if bw else 'NG'}")
    print('          P∥(上流,中流,下流)=',r['Pb'],' 頂点数',r['nv'],' ΣTu=',Su,' ΣTw=',Sw)
