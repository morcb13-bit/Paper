import numpy as np, sys, json, time
from tunnel_mirror import MB
tw=(sys.argv[1]=='twisted')
T=MB('mbox43.npz',twisted=tw); print('fallback count (all phases)',T.fallback)
ok_all=True; W=[]; acc=np.zeros((3,2),np.int64); t0=time.time()
for t in range(1300):
    ok,P=T.run_step(); ok_all&=bool(ok)
    if t>=300:
        acc+=P
        if (t-300)%50==49: W.append(acc.tolist()); acc[:]=0
json.dump(dict(ok=ok_all,W=W),open('d2_%s.json'%sys.argv[1],'w'))
print(sys.argv[1],'帳簿',ok_all,'sec',round(time.time()-t0))
