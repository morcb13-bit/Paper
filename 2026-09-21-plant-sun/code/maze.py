# 検定D0・A1  地図を被せた一翼（翼0、105環）で、餌の分子とアメーバの二つのオートマトン
# 事前登録（会話で決めたもの）
#  D0  分子が最初に届くコマの順 と 壁を避けた最短歩数の順 が一致（逆転する組が0）
#  A1  OK：餌に届く ∧ 先端がたどった歩数＝最短歩数 ∧ 伸びて縮む
#      NG：壁の前で止まる／全部埋まる／届かない
#  負の対照1 餌なし  2 分子の代わりに餌の向きだけ（太陽と同じ桁）  3 壁の隙間を反対側へ
import json,math,collections,sys
G=json.load(open('geo10.json'));R=G['R'];P=G['P']
W=[i for i,r in enumerate(R) if r[5]==0]; Wset=set(W)
own=collections.defaultdict(set)
for i in W:
  for k in R[i][4]: own[k].add(i)
ADJ={i:set() for i in W}
for s in own.values():
  for a in s:
    for b in s:
      if a!=b: ADJ[a].add(b)
TH=max(len(a) for a in ADJ.values())            # 配る閾値＝隣の最大数（4）。足りない分は外へ
def bfs(src,open_):
  d={src:0};q=collections.deque([src])
  while q:
    i=q.popleft()
    for j in ADJ[i]:
      if j in open_ and j not in d: d[j]=d[i]+1;q.append(j)
  return d
def mol_step(C,open_,src):
  C[src]+=1; out=0
  while True:
    F=[i for i in open_ if C[i]>=TH]
    if not F: break
    for i in F:
      C[i]-=TH
      for j in ADJ[i]:
        if j in open_: C[j]+=1
        else: out+=1
      out+=TH-len(ADJ[i])
  return out
def dig_angle(t):
  t=round(t,6); return 2 if t<36 else 1 if t<72 else 0 if t<108 else -1 if t<144 else -2
def ang(ax,ay,bx,by):
  d=abs(math.atan2(ay,ax)-math.atan2(by,bx))%(2*math.pi); return min(d,2*math.pi-d)*180/math.pi
def run(wall,food,T=3000,mode='mol',seed=(0,1,2)):
  open_=Wset-set(wall); C={i:0 for i in W}; first={}; outmol=0
  E={i:0 for i in W}; acc={i:0 for i in W}; sl={i:0 for i in W}
  for s in seed: E[s]=1
  births=deaths=0; hist=[]; reach=-1; bpar=collections.defaultdict(list); lastsend={}
  for t in range(T):
    if mode!='none': outmol+=mol_step(C,open_,food)
    for i in W:
      if C[i]>0 and i not in first: first[i]=t
    fx,fy=R[food][1],R[food][2]
    def val(j,i):   # 隣 j を 自分 i から見た桁
      if mode=='mol': return (C[j]>C[i])-(C[j]<C[i])
      if mode=='dir': return dig_angle(ang(R[j][1]-R[i][1],R[j][2]-R[i][2],fx-R[i][1],fy-R[i][2]))
      return 0
    alive=[i for i in W if E[i]]
    for i in alive:
      fr=[j for j in ADJ[i] if j in open_ and not E[j]]       # 表面＝まだ生きていない通れる隣
      s=1+sum(val(j,i) for j in fr)
      acc[i]+=s
      while acc[i]>=3:
        acc[i]-=5
        if fr:
          m=max(val(j,i) if mode!='mol' else C[j] for j in fr)
          for j in fr:
            if (val(j,i) if mode!='mol' else C[j])==m: sl[j]+=1; lastsend[j]=i
      while acc[i]<=-3:
        acc[i]+=5
        nb=[j for j in ADJ[i] if E[j]]
        if nb:
          m=min(C[j] for j in nb) if mode=='mol' else min(-val(j,i) for j in nb)
          for j in nb:
            if (C[j] if mode=='mol' else -val(j,i))==m: sl[j]-=1
    for j in W:
      while sl[j]>=3:
        sl[j]-=5
        if not E[j] and j in open_: E[j]=1;births+=1;bpar[j].append((t,lastsend.get(j)))
      while sl[j]<=-3:
        sl[j]+=5
        if E[j]: E[j]=0;deaths+=1
    n=sum(E.values()); hist.append(n)
    if reach<0 and E[food]: reach=t
  # 先端の路：餌から、誕生の送り手を時刻をさかのぼってたどる
  path=None
  if reach>=0:
    path=[food];cur=food;tt=reach
    while cur not in seed:
      rec=[p for (tb,p) in bpar[cur] if tb<=tt]
      if not rec or rec[-1] is None: path=None;break
      tb=max(tb for (tb,p) in bpar[cur] if tb<=tt); cur=rec[-1]; path.append(cur); tt=tb
      if len(path)>200: path=None;break
  return dict(first=first,reach=reach,path=path,hist=hist,births=births,deaths=deaths,E=E,open_=open_)
WALL_L=list(range(30,36))     # 段7：左端の2環（28,29）だけ隙間
WALL_R=list(range(28,34))     # 段7：右端の2環（34,35）だけ隙間
FOOD=97; SEED=(0,1,2)
def D0(wall):
  open_=Wset-set(wall); d=bfs(FOOD,open_); r=run(wall,FOOD,T=3000,mode='mol'); f=r['first']
  ks=[i for i in open_ if i in f and i in d]
  inv=sum(1 for a in ks for b in ks if d[a]<d[b] and f[a]>f[b])
  return len(open_),len(ks),inv,d
for name,wall in (('左の隙間',WALL_L),('右の隙間',WALL_R)):
  n,k,inv,d=D0(wall)
  ds=min(d[s] for s in SEED)
  print(f'D0 {name}  通れる環{n} 届いた環{k} 逆転する組{inv}  種→餌の最短歩数{ds}')
