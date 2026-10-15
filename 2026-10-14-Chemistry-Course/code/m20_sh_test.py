# m20: 検定 SH（SH_prereg.md のとおり）。整数と足し算・比較だけ。
import itertools
FULL={1:True,2:True,3:True,4:True,5:False,6:False,7:False}   # 記録（知識）
def nest(n):
    while n%2==0: n//=2
    while n%3==0: n//=3
    return n==1
def add(n):
    if n==1: return True
    return any(2*a+3*b==n for a in range(n) for b in range(n))
NS=range(1,8)
def score(pred): return sum(pred(n)==FULL[n] for n in NS)
N=score(nest); A=score(add)
print("NEST 予言:",{n:nest(n) for n in NS},"正答",N,"/7")
print("ADD  予言:",{n:add(n) for n in NS},"正答",A,"/7")
k=sum(nest(n) for n in NS)
tot=0;cnt=0;ge=0;dist={}
for S in itertools.combinations(NS,k):
    s=sum((n in S)==FULL[n] for n in NS); tot+=s; cnt+=1; ge+=s>=N; dist[s]=dist.get(s,0)+1
print(f"負の対照：{k} 個を選ぶ全 {cnt} 通り、総和 {tot}、分布 {dict(sorted(dist.items()))}、NEST 以上 {ge}/{cnt}")
c1=N>=A+2; c2=cnt*N>=tot+2*cnt
print(f"判定：{N} ≥ {A}+2 → {c1}　{cnt}·{N}={cnt*N} ≥ {tot}+{2*cnt}={tot+2*cnt} → {c2}　SH: {'OK' if c1 and c2 else 'NG'}")
