# Rotationally symmetric Ricci flow on a 3-manifold (S^2 cross-sections, n=2), Angenent-Knopf form
#   psi_t = psi_ss - (1 - psi_s^2)/psi,  phi_t = 2 (psi_ss/psi) phi,  d/ds = (1/phi) d/dx
# x in [0,1]: bulb equator | neck | bulb equator, reflecting ends.
import numpy as np, json, sys
N=int(sys.argv[1]) if len(sys.argv)>1 else 1601
x=np.linspace(0,1,N); dx=x[1]-x[0]
psi=1-0.74*np.cos(np.pi*(x-0.5))**4; phi=np.full(N,4.0)
def ext(a): return np.concatenate(([a[1]],a,[a[-2]]))
def D(psi,phi):
    P=ext(psi); F=ext(phi)
    px=(P[2:]-P[:-2])/(2*dx); fx=(F[2:]-F[:-2])/(2*dx); pxx=(P[2:]-2*psi+P[:-2])/dx**2
    return px/phi, pxx/phi**2-px*fx/phi**3
t=0.0; log=[]; snaps=[]; marks=[0.20,0.15,0.10,0.07,0.05,0.035,0.025,0.018,0.013,0.01]
while True:
    ps,pss=D(psi,phi); h=psi.min()
    if marks and h<=marks[0]:
        ds=phi*dx; i=psi.argmin(); s=np.cumsum(ds)-np.cumsum(ds)[i]
        sel=np.abs(s)<=4*h
        log.append(dict(t=t,h=float(h),K_neck=float((1-ps[i]**2)/h**2),
            max_slope_4h=float(np.abs(ps[sel]).max()),
            len_2h_over_h=float(ds[psi<2*h].sum()/h)))
        snaps.append([t,psi.tolist(),phi.tolist()]); marks.pop(0)
        if not marks: break
    dt=min(0.2*(phi.min()*dx)**2,0.01*h*h)
    psi,phi=psi+dt*(pss-(1-ps**2)/psi), phi+dt*2*(pss/psi)*phi
    t+=dt
for r in log: print({k:round(v,4) for k,v in r.items()})
json.dump(snaps,open("neck2_snaps.json","w"))
