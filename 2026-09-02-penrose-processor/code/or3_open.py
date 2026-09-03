# 検定OR3-O  (1) W から次の W が決まらない型はどれだけあるか（直接の検定）
#             (2) M を多く抱える型は、語を多く持っているだけか（組み合わせか幾何か）
#             (3) 最大の型を開く
import json, math, sys
from collections import defaultdict, Counter
exec(open('or3_why.py').read().split('print("\\n■ (a)')[0].replace('NMAX=int(sys.argv[1]) if len(sys.argv)>1 else 16','NMAX=16'))
N=16
root=tuple((k,EI[allst[i]]) for k,i in enumerate(reps)); NR=len(reps)
Wsz={}; Ms=defaultdict(set); Aw=Counter(); kids=defaultdict(set); Wstarts={}
best=[0,None]
def dfs(al,d,w):
    if d==N:
        ws=tuple(i for i,_ in al); wh=hash(ws)
        Wsz[wh]=len(ws); Wstarts[wh]=ws; Aw[wh]+=1; Ms[wh].add(hash(al))
        tri=[]
        for x in (0,1,2):
            n2=tuple((i,t) for i,e in al if (t:=T[e][x])>=0)
            tri.append(hash(tuple(i for i,_ in n2)) if n2 else None)
        kids[wh].add(tuple(tri))
        return
    for x in (0,1,2):
        n2=tuple((i,t) for i,e in al if (t:=T[e][x])>=0)
        if n2: dfs(n2,d+1,w+(x-1,))
sys.setrecursionlimit(80); dfs(root,0,())
NW=len(Aw)
print("長さ%d・代表%d本：W %d 種／語 %d 本／M %d 種" % (N,NR,NW,sum(Aw.values()),sum(len(v) for v in Ms.values())))

print("\n■ (1) W から次の W が決まるか")
bad=sum(1 for wh in kids if len(kids[wh])>1)
badA=sum(Aw[wh] for wh in kids if len(kids[wh])>1)
print("  次の W の組が一通りに決まらない型 %d / %d（%.1f%%）" % (bad,NW,100*bad/NW))
print("  その型に属する語 %d / %d（%.1f%%）" % (badA,sum(Aw.values()),100*badA/sum(Aw.values())))
print("  → W(w)=W(w') でも W(wd)≠W(w'd) となる型が実在する" if bad else "  → 反例なし")

print("\n■ (2) M が多いのは語が多いだけか")
rows=[(len(Ms[wh]),Aw[wh],Wsz[wh],wh) for wh in Aw]
rows.sort(reverse=True)
print("  M の多い順に上位8つ")
print("   M の個数   語の本数  M/語   |W|/代表")
for m,a,s,wh in rows[:8]:
    print("   %8d %10d %6.3f %8.3f" % (m,a,m/a,s/NR))
b=defaultdict(list)
for m,a,s,wh in rows: b[min(int(5*s/NR),4)].append((m,a))
lab={0:"0〜20%",1:"20〜40%",2:"40〜60%",3:"60〜80%",4:"80〜100%"}
print("  型の広さごとの M/語（1に近いほど M が潰れていない）")
for k in sorted(b):
    v=b[k]; print("   %-9s W %5d  語/型 平均%8.1f  M/語 平均%.3f 中央%.3f"
        % (lab[k],len(v),sum(a for _,a in v)/len(v),
           sum(m/a for m,a in v)/len(v), sorted(m/a for m,a in v)[len(v)//2]))

print("\n■ (3) 最大の型を開く")
m,a,s,wh=rows[0]
ws=Wstarts[wh]
print("  M %d 種／語 %d 本／M/語 %.3f／|W| = %d 本（代表%d本中 %.1f%%）" % (m,a,m/a,len(ws),NR,100*len(ws)/NR))
sh=[R for R,sts in SS[:21] for _ in sts]
print("  生き残った出発の殻 %s" % dict(Counter(sh[reps[i]] for i in ws)))
print("  次の W の組の通り数 %d" % len(kids[wh]))
