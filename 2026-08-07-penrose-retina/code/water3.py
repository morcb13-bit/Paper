#  検定TP3  閉じているかを水で読む
#      網の縁から、点いていないセルだけを通って水を入れる。届かないセルが残れば閉じている。
#      使うのは近接（沈みと同じ2.7）と伝播だけ。隣接表も新しい閾値も持たない。
#      測る量：届かなかったセルの数（＝内部）／内部のかたまりの数
#      OK なら：閉じた輪は内部あり、開いた輪・線・分かれ道は内部なし、二重の輪は内部2つ
#      NG なら：区別がつかない
exec(open('sink.py').read().split('print("■ 沈んだところへ移る')[0])
import math
R = 2.7

def paint(path):
    return set(min(CELLS, key=lambda q: (XYC[q][0]-x)**2+(XYC[q][1]-y)**2) for x, y in path)
def seg(a,b,n=120):
    return [(a[0]+(b[0]-a[0])*i/n, a[1]+(b[1]-a[1])*i/n) for i in range(n+1)]
def arc(c,Rr,t0,t1,n=200):
    return [(c[0]+Rr*math.cos(t0+(t1-t0)*i/n), c[1]+Rr*math.sin(t0+(t1-t0)*i/n)) for i in range(n+1)]
X, Y = O
FIG = {
    "一本の線":      seg((X-7,Y-4),(X+7,Y+4)),
    "離れた二本":    seg((X-8,Y-6),(X-1,Y-6))+seg((X+2,Y+5),(X+9,Y+5)),
    "閉じた輪(三角)": seg((X-6,Y-4),(X+6,Y-4))+seg((X+6,Y-4),(X,Y+7))+seg((X,Y+7),(X-6,Y-4)),
    "丸 R=6.5":     arc((X,Y),6.5,0,2*math.pi),
    "丸 R=5.0":     arc((X,Y),5.0,0,2*math.pi),
    "丸 R=4.0":     arc((X,Y),4.0,0,2*math.pi),
    "丸 R=3.0":     arc((X,Y),3.0,0,2*math.pi),
    "開いた輪(C字)": arc((X,Y),6.5,0.9,2*math.pi-0.9),
    "分かれ道(Y字)": seg((X,Y-7),(X,Y))+seg((X,Y),(X-6,Y+6))+seg((X,Y),(X+6,Y+6)),
    "二重の輪":      arc((X,Y),4.0,0,2*math.pi)+arc((X,Y),9.0,0,2*math.pi),
}
RNET = 14.0
NET = [q for q in CELLS if math.hypot(XYC[q][0]-X, XYC[q][1]-Y) <= RNET]
RIM = [q for q in NET if math.hypot(XYC[q][0]-X, XYC[q][1]-Y) >= RNET-2.0]

def flood(on):
    free = [q for q in NET if q not in on]
    wet = set(q for q in RIM if q not in on)
    st = list(wet)
    while st:
        u = st.pop()
        for b in free:
            if b not in wet and (XYC[u][0]-XYC[b][0])**2+(XYC[u][1]-XYC[b][1])**2 <= R*R:
                wet.add(b); st.append(b)
    dry = [q for q in free if q not in wet]
    seen = set(); comps = 0
    for a in dry:
        if a in seen: continue
        comps += 1; st=[a]; seen.add(a)
        while st:
            u = st.pop()
            for b in dry:
                if b not in seen and (XYC[u][0]-XYC[b][0])**2+(XYC[u][1]-XYC[b][1])**2 <= R*R:
                    seen.add(b); st.append(b)
    return len(dry), comps

print("■ 検定TP3　網の縁から水を入れる（網の半径 %.0f・網の中のセル %d枚）\n" % (RNET, len(NET)))
print("  図形             点いたセル   届かなかったセル   内部のかたまり")
for nm, path in FIG.items():
    on = paint(path) & set(NET)
    dry, c = flood(on)
    print("  %-14s     %3d           %3d              %2d" % (nm, len(on), dry, c))
