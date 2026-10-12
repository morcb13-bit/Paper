# フィボナッチで重ねる正四面体の殻：整数だけの検算
V=[(1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1)]
def comps(m):
    return [] if m<0 else [(a,b,c,m-a-b-c) for a in range(m+1) for b in range(m+1-a) for c in range(m+1-a-b)]
print("== 1. 辺 n の正四面体を辺 1 の正四面体と正八面体で埋める（6×体積：正四面体 1、正八面体 4） ==")
for n in range(1,6):
    up,down,octa=len(comps(n-1)),len(comps(n-3)),len(comps(n-2))
    print(f" 辺{n}: 上向き{up:3} 下向き{down:3} 八面体{octa:3}  6×体積 {up+down+4*octa:4} = n³ {n**3:4}  {'一致' if up+down+4*octa==n**3 else '不一致'}")
print("\n== 2. 正立と鏡像を交互に重ねる（足し算で作るフィボナッチ） ==")
F=[1,1]
while len(F)<8: F.append(F[-1]+F[-2])
for i in range(1,8):
    a,b=F[i-1],F[i]; print(f" {i}→{i+1}: 辺 {a}→{b}  向き {'A→B' if i%2 else 'B→A'}  b²−ab−a² = {b*b-a*b-a*a:2}  差 {b-a}")
print("\n== 3. 殻の席：一辺 n の網の点（a+b+c+d=n でどれか一つが 0） ==")
SH=["H","He","L","M","N","O","P","Q"]; cum=0
for i,sh in enumerate(SH):
    n=1 if i<2 else i//2+1
    pts=sum(1 for w in comps(n) if 0 in w)
    e=1 if i<2 else 2*n*n
    cum+=e
    print(f" {sh:2} 辺{F[i]:2} 層{n} 網の点 {pts:2} = 2n²+2 {2*n*n+2:2}  電子 {e:2}  累計 {cum:3}")
print("\n== 4. L と M の外接球（頂点 (s,s,s)） ==")
for nm,s,bohr in [("L",2,4),("M",3,9)]:
    print(f" {nm}: 辺 {s}  外接球の半径² 3s² = {3*s*s}  ボーアの n² = {bohr}")
print(" 12 : 27 と 4 : 9 の交差掛け算:", 12*9, "=", 27*4, "->", 12*9==27*4)
