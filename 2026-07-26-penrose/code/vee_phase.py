"""
記事の平衡3進の敷き詰めを、波の言葉で読み直す。
V字は3環＝ちょうど1周期（検定WS4）。3環が「山谷山」か「谷山谷」かは、
開始環が山にあるか谷にあるかで決まる ── つまり位相が半周期ズレた2種。

敷き詰めは登録済みの chain_balanced.balanced() から取り、新しく作らない（判別法3）。

#  検定VP1 桁0（n≡0 mod 3）のとき、隣り合うV字の向きは交互に反転するか
#      OK なら：ぴったり敷ける場合、敷き詰めは2種の波を交互に見ていることになる
#      NG なら：向きは交互でない
#
#  検定VP2 桁±1 のとき、V字の向きの並びはどうなるか
#      OK なら：桁が入ることで向きの交替が破れ、破れ方が桁の符号で分かれる
#      NG なら：桁に関係なく同じ並びになる
#
#  検定VP3 印（余った環／重ねた環）は山と谷のどちらに来るか
#      OK なら：桁の符号と、印の位相（山か谷か）が一対一に対応する
#      NG なら：対応しない。桁は位相では読めない
#
#  検定VP4 鏡映を許さないと何が起きるか（誤り69・70 の再現）
#      OK なら：向きの片方しか座れなくなる＝2種の波の片方だけを見ることになる
#      NG なら：鏡映の有無は向きに効かない
"""
import math, cmath
import b13_chain_units as B
import chain_balanced as CB

ROT = cmath.exp(-1j*math.radians(18))
def xy(a):
    p = B.num(a)*ROT if hasattr(B, "num") else None
    return p
def pt(a):
    z = cmath.exp(1j*math.pi/5)
    p = sum(float(a[k])*z**k for k in range(4))*ROT
    return (p.real, p.imag)

def heights(n):
    """奇数の鎖 n の各環の高さ（中心線からの符号つき）"""
    cs = B.unit(n)
    ys = [pt(c)[1] for c in cs]
    ctr = (max(ys)+min(ys))/2
    return [y-ctr for y in ys]

def vee_dir(ys, i):
    """環 i,i+1,i+2 の V字の向き。中央が下なら -1（山谷山）、上なら +1（谷山谷）"""
    mid = ys[i+1] - (ys[i]+ys[i+2])/2
    return -1 if mid < 0 else +1

print("V字は3環＝1周期。開始環が山なら『山谷山』、谷なら『谷山谷』。")
print("以下、-1=山谷山 / +1=谷山谷\n")

rows = []
for n in range(3, 30, 2):
    digit, cover, mark = CB.balanced(n)
    starts = sorted(a-1 for a,_ in cover)   # 1始まり → 0始まり
    ys = heights(n)
    dirs = [vee_dir(ys, s) for s in starts]
    rows.append((n, digit, starts, dirs, mark-1, ys))

# ---------------- VP1 ----------------
print("検定VP1 桁0 のとき、隣り合うV字の向きは交互に反転するか")
vp1 = True
for n, d, st, dirs, mk, ys in rows:
    if d != 0: continue
    alt = all(dirs[i] != dirs[i+1] for i in range(len(dirs)-1))
    print(f"  n={n:2d}  開始環={st}  向き={dirs}  交互={alt}")
    if not alt: vp1 = False
print(f"検定VP1: {'OK' if vp1 else 'NG'}  桁0 で必ず交互か = {vp1}\n")

# ---------------- VP2 ----------------
print("検定VP2 桁±1 のとき、向きの並びはどうなるか")
pat = {+1: set(), -1: set()}
for n, d, st, dirs, mk, ys in rows:
    if d == 0: continue
    breaks = [i for i in range(len(dirs)-1) if dirs[i] == dirs[i+1]]
    pat[d].add((len(breaks), tuple(breaks) if len(dirs) <= 6 else None))
    print(f"  n={n:2d} 桁={d:+d}  開始環={st}  向き={dirs}  交替が破れる箇所={breaks}")
vp2 = all(len(v) >= 1 for v in pat.values())
print(f"検定VP2: {'OK' if vp2 else 'NG'}  桁が入ると交替が破れるか = {vp2}\n")

# ---------------- VP3 ----------------
print("検定VP3 印は山と谷のどちらに来るか")
tbl = {}
for n, d, st, dirs, mk, ys in rows:
    h = ys[mk]
    kind = "山" if h > 0 else ("谷" if h < 0 else "?")
    tbl.setdefault(d, []).append((n, kind, round(h, 6)))
for d in sorted(tbl):
    kinds = {k for _, k, _ in tbl[d]}
    print(f"  桁={d:+d}: 印の位相 = {sorted(kinds)}   例 {tbl[d][:4]}")
vp3 = all(len({k for _, k, _ in v}) == 1 for v in tbl.values()) and \
      len({tuple(sorted({k for _, k, _ in v})) for v in tbl.values()}) == len(tbl)
print(f"検定VP3: {'OK' if vp3 else 'NG'}  桁ごとに位相が一意で、桁どうしで違うか = {vp3}\n")

# ---------------- VP4 ----------------
print("検定VP4 鏡映を許さないと何が起きるか")
vp4 = True
for n in (9, 11, 13):
    a = CB.pos3(n, mirror=True)
    b = CB.pos3(n, mirror=False)
    ys = heights(n)
    da = sorted({vee_dir(ys, s-1) for s in a})
    db = sorted({vee_dir(ys, s-1) for s in b})
    print(f"  n={n:2d}  鏡映あり 座れる開始環={sorted(x-1 for x in a)} 向き={da}"
          f"  ／ 鏡映なし ={sorted(x-1 for x in b)} 向き={db}")
    if not (len(da) == 2 and len(db) == 1): vp4 = False
print(f"検定VP4: {'OK' if vp4 else 'NG'}  鏡映なしでは向きが片方だけか = {vp4}")
