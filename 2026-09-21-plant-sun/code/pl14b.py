src=open('pl14.py').read(); exec(src[:src.index('ymax=')])
def run2(start,sched,days):
    live={start}; body={start}; mx=1; out=[]
    for d in range(days):
        for deg in sched:
            if deg is None: continue
            nxt=set()
            for i in live:
                tg=tips_par(i,deg)
                if len(tg)>=3: nxt.add(i); continue
                dest={step_to(i,t) for t in tg}; dest.discard(None)
                nxt|= dest if dest else {i}
            live=nxt; body|=live; mx=max(mx,len(live))
        out.append((len(body),mx,sorted(live)))
    return out
ymax=max(y for _,y in SC); ground=[i for i in range(len(SC)) if abs(SC[i][1]-ymax)<1e-6]
top=min(range(len(SC)),key=lambda i:SC[i][1])
for name,sched in (("東→西",DAY),("もやし",[270]*11+[None,None])):
    print(name, "頂点の星", top)
    for s in ground:
        r=run2(s,sched,3)
        print(f"   出発{s:2d}  日ごと（体の星数, 同時の先端の最大, 夕方の先端）", [(b,m,l) for b,m,l in r])
