import json
from judge_d0 import zm, za, zs, sgn
def zabs(x): return x if sgn(x)>=0 else zs(x,-1)
R={}
for m in ['still_s','still1','still_l']:
    r=json.load(open('d0_%s.json'%m)); S=(0,0)
    for t in r['Tw']: S=za(S,tuple(t))
    R[m]=(zabs(S),sum(r['nv']))
    phi=(1+5**0.5)/2
    print(m,'N=',sum(r['nv']),'ΣT_v=',S,'(概算',round((S[0]+S[1]*phi)/1e6,1),'百万)  b概算',round(abs(S[0]+S[1]*phi)/sum(r['nv']),1))
def ratio_in(Bs,Ns,Bl,Nl,pred_num,pred_den):
    # b_s/b_l が pred*(0.6..1.4) に入るか:  b_s/b_l = Bs*Nl/(Bl*Ns)
    L=zs(Bs,Nl*pred_den*10); Rr=zs(Bl,Ns*pred_num)
    return sgn(za(L,zs(Rr,-6)))>0 and sgn(za(L,zs(Rr,-14)))<0
(Bs,Ns),(Bm,Nm),(Bl,Nl)=R['still_s'],R['still1'],R['still_l']
print('縁説 (比 2 と 1.5):', ratio_in(Bs,Ns,Bm,Nm,2,1), ratio_in(Bm,Nm,Bl,Nl,3,2))
print('中身説 (比 1 と 1):', ratio_in(Bs,Ns,Bm,Nm,1,1), ratio_in(Bm,Nm,Bl,Nl,1,1))
