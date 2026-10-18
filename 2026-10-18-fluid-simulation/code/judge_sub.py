import json
from judge_d0 import zm, za, zs, sgn
def load(m): return json.load(open('d0_%s.json'%m))
def test(T):
    S=(0,0);Q=(0,0)
    for t in T: S=za(S,t); Q=za(Q,zm(t,t))
    return sgn(za(zs(Q,9),zs(zm(S,S),-1)))>0, S
base=load('still1'); phi=(1+5**0.5)/2
for m in ['normal','rev','nocycle']:
    r=load(m); out=[]
    for key in ['Tu','Tw']:
        d=[(a[0]-b[0],a[1]-b[1]) for a,b in zip(r[key],base[key])]
        ok,S=test(d); out.append(f"{key}: {'OK' if ok else 'NG'} (ΣΔ≈{(S[0]+S[1]*phi)/1e6:+.1f}百万)")
    print(m,' / '.join(out))
