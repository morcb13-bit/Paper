"""一歩の規則そのものを検査する。担体は使わない。次数だけで書ける。
   検定W1 二乗和が144倍になる次数の条件
   検定W2 出た側から入った側を戻せるか（後戻り）
   必ず落ちる設定：c を 24/d 以外にしたら W1 が落ちること"""
from fractions import Fraction as F
import random
rnd=random.Random(7)

def step(x, c):
    S=sum(x); return [c*S-12*xi for xi in x]

print("検定W1  out = c·S − 12·x の二乗和 / 入力の二乗和")
for d in range(2,13):
    for c in (F(24,d), F(24,d)+1, F(24,d)-1):
        ok=True
        for _ in range(50):
            x=[rnd.randint(-9,9) for _ in range(d)]
            if sum(v*v for v in x)==0: continue
            y=step(x,c)
            if sum(v*v for v in y)!=144*sum(v*v for v in x): ok=False; break
        tag="◀ 24/d" if c==F(24,d) else ""
        if ok or tag: print(f"  次数{d:>2} c={str(c):>6}: 144倍 {ok} {tag}")

print("\n  代数：Σ(cS−12x_i)² = (d·c²−24c)S² + 144Σx² なので、d·c² = 24c すなわち c = 24/d")
print("  c が整数であるための条件：d が 24 を割ること。担体の次数 2・3・4 はすべて割る。")

print("\n検定W2  出た側から入った側を戻す")
for d in (2,3,4):
    c=F(24,d); bad=0
    for _ in range(200):
        x=[rnd.randint(-50,50) for _ in range(d)]
        y=step(x,c)
        S=F(sum(y),12)                      # Σout = (d·c − 12)·S = 12·S
        back=[(c*S-yi)/12 for yi in y]
        if [F(v) for v in x]!=back: bad+=1
    print(f"  次数{d}: Σout = 12·S を使って復元。食い違い {bad} 件 / 200")

print("\n検定W3  出力の偶奇")
for d in (2,3,4):
    c=24//d; odd=0
    for _ in range(500):
        x=[rnd.randint(-99,99) for _ in range(d)]
        if any(v%2 for v in step(x,c)): odd+=1
    print(f"  次数{d}: 奇数が出た回数 {odd} / 500  （c={c}、cS も 12x も偶数か）")
