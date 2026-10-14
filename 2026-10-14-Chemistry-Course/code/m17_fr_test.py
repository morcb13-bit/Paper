# m17: 検定 FR（事前登録 FR_prereg.md のとおり）。判定は整数の比較だけ。
import random

V = {}
for e in "Li Na K Rb Cs".split(): V[e] = 1
for e in "Be Mg Ca Sr Ba".split(): V[e] = 2
for e in "B Al Ga In Tl".split(): V[e] = 3
for e in "C Si Ge Sn Pb".split(): V[e] = 4
for e in "N P As Sb Bi".split(): V[e] = 5
for e in "O S Se Te Po".split(): V[e] = 6
for e in "F Cl Br I At".split(): V[e] = 7

ELEM = {}
for e in "Li Na K Rb Cs Be Mg Ca Sr Ba Al Ga In Tl Sn Pb Bi".split(): ELEM[e] = "金属"
for e in "C N P O S Se F Cl Br I".split(): ELEM[e] = "共有"
ELEM_UNDECIDED = "B Si Ge As Sb Te Po At".split()

COMP = [  # (名前, 元素1, 元素2, 記録)
    ("NaCl","Na","Cl","イオン"),("KBr","K","Br","イオン"),("LiF","Li","F","イオン"),("MgO","Mg","O","イオン"),
    ("CaF2","Ca","F","イオン"),("CaO","Ca","O","イオン"),("Na2O","Na","O","イオン"),("K2S","K","S","イオン"),
    ("MgCl2","Mg","Cl","イオン"),("BaCl2","Ba","Cl","イオン"),
    ("CO2","C","O","共有"),("CS2","C","S","共有"),("CCl4","C","Cl","共有"),("SiO2","Si","O","共有"),
    ("SiC","Si","C","共有"),("NCl3","N","Cl","共有"),("PCl3","P","Cl","共有"),("SO2","S","O","共有"),
    ("SF6","S","F","共有"),("OF2","O","F","共有"),("BN","B","N","共有"),("GaAs","Ga","As","共有"),
    ("NaK","Na","K","金属"),("CaMg2","Ca","Mg","金属"),("Mg17Al12","Mg","Al","金属"),
]
COMP_UNDECIDED = ["AlCl3","BeCl2","AlN","LiAl"]

def side(v): return "D" if v <= 3 else "S"
def pred_elem(v): return "金属" if side(v) == "D" else "共有"
def pred_comp(v1, v2):
    a, b = side(v1), side(v2)
    if a == b == "D": return "金属"
    if a == b == "S": return "共有"
    return "イオン"

def score(vmap, show=False):
    se = 0
    for e, rec in ELEM.items():
        p = pred_elem(vmap[e]); ok = p == rec; se += ok
        if show and not ok: print(f"   外れ 単体 {e}(v={vmap[e]}): 予言 {p} / 記録 {rec}")
    sc = 0
    for n, a, b, rec in COMP:
        p = pred_comp(vmap[a], vmap[b]); ok = p == rec; sc += ok
        if show and not ok: print(f"   外れ 化合物 {n}: 予言 {p} / 記録 {rec}")
    return se, sc

print("== FR")
fe, fc = score(V, show=True)
FR = fe + fc
print(f"  単体 {fe}/{len(ELEM)}  化合物 {fc}/{len(COMP)}  合計 {FR}/{len(ELEM)+len(COMP)}")

print("== 基準線（すべて共有）")
be = sum(r == "共有" for r in ELEM.values()); bc = sum(r[3] == "共有" for r in COMP)
BASE = be + bc
print(f"  単体 {be}  化合物 {bc}  合計 {BASE}")

print("== 負の対照（v の並べ替え、シード 13、10000 回）")
used = sorted(set(ELEM) | {x for c in COMP for x in c[1:3]})
rng = random.Random(13)
T = 10000
tot = tot_e = tot_c = 0; ge = 0; hist = {}; lo = 10**9; hi = -1
vals = [V[e] for e in used]
for _ in range(T):
    rng.shuffle(vals)
    vm = dict(zip(used, vals))
    a, b = score(vm); s = a + b
    tot += s; tot_e += a; tot_c += b
    hist[s] = hist.get(s, 0) + 1
    lo = min(lo, s); hi = max(hi, s)
    if s >= FR: ge += 1
print(f"  合計の総和 {tot}（平均 ≈ {tot//T}.{(tot%T)*100//T:02d}）  単体の総和 {tot_e}  化合物の総和 {tot_c}")
print(f"  最小 {lo}  最大 {hi}")
print("  分布:", " ".join(f"{k}:{hist[k]}" for k in sorted(hist)))
print(f"  併記：FR（{FR}）以上を取った回数 {ge}/{T}")

print("== 判定")
c1 = FR >= BASE + 3
c2 = T * FR >= tot + 3 * T
print(f"  基準線より 3 以上多い：{FR} ≥ {BASE}+3 → {c1}")
print(f"  負の対照の平均より 3 以上多い：{T}·{FR} = {T*FR} ≥ {tot} + {3*T} = {tot+3*T} → {c2}")
print("  FR:", "OK" if c1 and c2 else "NG")
print(f"  判定不能：単体 {len(ELEM_UNDECIDED)}（{' '.join(ELEM_UNDECIDED)}）、化合物 {len(COMP_UNDECIDED)}（{' '.join(COMP_UNDECIDED)}）")
