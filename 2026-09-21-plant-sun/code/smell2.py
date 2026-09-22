# 刺激＝匂い（3次元を漂う分子）／感覚器＝隣と自分の大小／体＝core10 の規則
#  担体  翼0（105環）、隣＝五角形を共有する環。地面の層0に壁（段7、隙間は左の2環 or 右の2環）
#  分子  L層に重ねた同じ翼。各分子は毎コマ6口（面内4・上・下）から一つを等しく選ぶ。
#        面内の口の先：通れる環へ移る／層0の壁なら留まる／隣が4に足りない分は外へ
#        下：層0では留まる  上：最上層では外へ   餌（層0・71番）が毎コマ1個出す
#  層1枚の対照は面内4口だけ（D0 と同じ）
import numpy as np, collections
exec(open('maze.py').read().split("def mol_step")[0])
FOOD=71; SEED=(0,1,2); WALL_L=list(range(30,36)); WALL_R=list(range(28,34))
idx={i:n for n,i in enumerate(W)}; n0=len(W)
NB=np.full((n0,4),-1)
for i in W:
  for k,j in enumerate(sorted(ADJ[i])): NB[idx[i],k]=idx[j]
class Smell:
  def __init__(s,wall,L,seed):
    s.L=L; s.rng=np.random.default_rng(seed); s.C=np.zeros((L,n0),dtype=np.int64)
    s.wall=np.zeros(n0,bool); s.wall[[idx[w] for w in wall]]=True
    s.k=6 if L>1 else 4
  def step(s):
    L,C=s.L,s.C; C[0,idx[FOOD]]+=1
    N=np.zeros_like(C)
    for l in range(L):
      m=s.rng.multinomial(C[l],[1/s.k]*s.k)          # 各環の分子を口に分ける
      for k in range(4):
        tgt=NB[:,k]; ok=tgt>=0
        if l==0: blocked=ok&s.wall[np.where(ok,tgt,0)]; mv=ok&~blocked
        else: blocked=np.zeros(n0,bool); mv=ok
        np.add.at(N[l],tgt[mv],m[mv,k]); N[l][blocked]+=m[blocked,k]   # 壁は留まる、-1 は外へ
      if s.k==6:
        if l+1<L: N[l+1]+=m[:,4]                          # 上（最上層は外へ）
        if l>0: N[l-1]+=m[:,5]
        else: N[0]+=m[:,5]                                # 地面より下へは行けない
    N[0][s.wall]=0
    s.C=N
    return N[0]
def run(wall,L=3,T=20000,seed=1,drive=1,sense=True,body=True,M=None):
  open_=Wset-set(wall); sm=Smell(wall,L,seed)
  E={i:0 for i in W}; acc=dict.fromkeys(W,0); sl=dict.fromkeys(W,0)
  born_t={}
  for s_ in SEED: E[s_]=1; born_t[s_]=-1
  S=np.zeros(n0,dtype=np.int64); hist=[]; reach=-1; bpar=collections.defaultdict(list); last={}; births=deaths=0
  for t in range(T):
    g=sm.step(); S+=g
    if not body: continue
    c=lambda j:int(g[idx[j]]) if sense else 0
    for i in [i for i in W if E[i]]:
      fr=[j for j in ADJ[i] if j in open_ and not E[j]]
      ci=c(i); acc[i]+=drive+sum((c(j)>ci)-(c(j)<ci) for j in fr)
      while acc[i]>=3:
        acc[i]-=5
        if fr:
          m=max(c(j) for j in fr)
          for j in fr:
            if c(j)==m: sl[j]+=1; last[j]=i
      while acc[i]<=-3:
        acc[i]+=5; nb=[j for j in ADJ[i] if E[j]]
        if nb:
          m=min(c(j) for j in nb)
          for j in nb:
            if c(j)==m: sl[j]-=1
    for j in W:
      while sl[j]>=3:
        sl[j]-=5
        if not E[j] and j in open_:
          if M is not None and sum(E.values())>=M:          # 体の量は M まで：先が伸びたらしっぽ（いちばん古い環）が縮む
            tail=min((i for i in W if E[i]),key=lambda i:born_t[i]); E[tail]=0; deaths+=1
          E[j]=1;births+=1;born_t[j]=t;bpar[j].append((t,last.get(j)))
      while sl[j]<=-3:
        sl[j]+=5
        if E[j]: E[j]=0;deaths+=1
    hist.append(sum(E.values()))
    if reach<0 and E[FOOD]: reach=t
  path=None
  if reach>=0:
    path=[FOOD];cur=FOOD;tt=reach
    while cur not in SEED:
      rec=[(tb,p) for (tb,p) in bpar[cur] if tb<=tt]
      if not rec or rec[-1][1] is None: path=None;break
      tt,cur=rec[-1]; path.append(cur)
      if len(path)>200: path=None;break
  return dict(S={i:int(S[idx[i]]) for i in W},reach=reach,path=path,hist=hist,births=births,deaths=deaths,E=E)
