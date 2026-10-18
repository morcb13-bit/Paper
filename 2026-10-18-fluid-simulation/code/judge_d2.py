import json
from judge_d0 import zm, za, zs, sgn
phi=(1+5**0.5)/2
for m in ['mirror','twisted']:
    r=json.load(open('d2_%s.json'%m)); res=[]
    for a,nm in enumerate('xyz'):
        T=[tuple(w[a]) for w in r['W']]; S=(0,0); Q=(0,0)
        for t in T: S=za(S,t); Q=za(Q,zm(t,t))
        ok=sgn(za(zs(Q,9),zs(zm(S,S),-1)))>0
        res.append(f"P_{nm}: {'OK' if ok else 'NG'}（Σ≈{(S[0]+S[1]*phi)/1e6:+.2f}百万）")
    print(m,'帳簿',r['ok'],' / '.join(res))
