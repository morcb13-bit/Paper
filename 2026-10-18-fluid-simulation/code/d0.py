import numpy as np, sys, json, time
from tunnel3d import *
fname=sys.argv[1]; mode=sys.argv[2]
wl=(0,) if mode!='rev' else (6,)
T=T3(fname,base=4,wind=(0 if mode.startswith('still') else 8),windlab=wl,cycle=(mode!='nocycle'))
ok_all=True; Tu=[];Tw=[]; acc_u=np.zeros(2,np.int64); acc_w=np.zeros(2,np.int64)
bands=[np.where(T.band[:T.NI]==b)[0] for b in range(3)]
Pb=[np.zeros(2,np.int64) for _ in range(3)]
t0=time.time()
for t in range(1300):
    ok_all&=bool(T.run_step())
    if t>=300:
        n=T.n[:T.NI]; p,u,w=T.mom(n); acc_u+=u; acc_w+=w
        for b in range(3): Pb[b]+=T.mom(T.n[bands[b]])[0]
        if (t-300)%50==49: Tu.append(acc_u.tolist()); Tw.append(acc_w.tolist()); acc_u[:]=0; acc_w[:]=0
nv=[len(b) for b in bands]
json.dump(dict(mode=mode,ok=ok_all,Tu=Tu,Tw=Tw,Pb=[x.tolist() for x in Pb],nv=nv,sec=time.time()-t0),open('d0_%s.json'%mode,'w'))
print(mode,'conserve',ok_all,'Pband',[x.tolist() for x in Pb],'nv',nv,'sec',round(time.time()-t0))
