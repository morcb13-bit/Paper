# 結合の図を 4 枚（一図一ファイル）作る。座標は m15_share_vef.py と同じ整数。
import json, pathlib

A = [(3, 3, 3), (3, -3, -3), (-3, 3, -3), (-3, -3, 3)]
def add(a, b): return tuple(x + y for x, y in zip(a, b))
def sub(a, b): return tuple(x - y for x, y in zip(a, b))
def dot(a, b): return sum(x * y for x, y in zip(a, b))
def sc(k, a): return tuple(k * x for x in a)
def pref(p, m2): return sub(m2, p)
def plref(p, n, d):
    t = 2 * (dot(n, p) - d); nn = dot(n, n); assert t % nn == 0
    return sub(p, sc(t // nn, n))

C0 = dict(el="C", center=(0, 0, 0), verts=A)
def by_point(m2, el): return dict(el=el, center=m2, verts=[pref(p, m2) for p in A])
def by_face(el):
    n, d = (1, 1, -1), 3
    return dict(el=el, center=plref((0, 0, 0), n, d), verts=[plref(p, n, d) for p in A])

def finish(atoms, fill):
    # 共有した頂点 → bond、それ以外 → fill（H か lp）
    for i, a in enumerate(atoms):
        arms = []
        for v in a["verts"]:
            shared = any(v in b["verts"] for j, b in enumerate(atoms) if j != i)
            arms.append("bond" if shared else fill[a["el"]])
        a["arms"] = arms
    for a in atoms:
        a["center"] = list(a["center"]); a["verts"] = [list(v) for v in a["verts"]]
    return atoms

def tally(atoms):
    t = {"bond": 0, "H": 0, "lp": 0}
    for a in atoms:
        for s in a["arms"]: t[s] += 1
    return t

FIGS = []

def scene(atoms, fill):
    at = finish(atoms, fill)
    return dict(atoms=at, tally=tally(at), len2=dot(at[1]["center"], at[1]["center"]))

single = scene([dict(C0), by_point(sc(2, A[0]), "C")], {"C": "H"})
double = scene([dict(C0), by_point(add(A[0], A[1]), "C")], {"C": "H"})
triple = scene([dict(C0), by_face("C")], {"C": "H"})
co = scene([dict(C0), by_face("O")], {"C": "lp", "O": "lp"})
co2 = scene([dict(C0), by_point(add(A[0], A[1]), "O"), by_point(add(A[2], A[3]), "O")], {"C": "lp", "O": "lp"})

for name, s in [("single", single), ("double", double), ("triple", triple), ("co", co), ("co2", co2)]:
    print(name, s["tally"], "長さ²", s["len2"], [a["arms"] for a in s["atoms"]])

PAGES = [
  dict(file="bond1_single.html", title="単結合の腕",
       h1="単結合 ― 腕 1 本ずつで手をつなぐ",
       lead="炭素の正四面体には、中心から頂点へ向かう腕が 4 本ある。二つの炭素が頂点を一つ共有すると、それぞれの腕が 1 本ずつその頂点で出会う。これが単結合。残った腕 3 本ずつに水素がついて、エタン C₂H₆ の形になる。",
       scenes=[dict(cap="エタン C₂H₆", s=single)],
       count="腕の勘定：炭素 2 個で 8 本 ＝ 炭素どうし 2 本（1 組）＋ 水素 6 本"),
  dict(file="bond2_double.html", title="二重結合の腕",
       h1="二重結合 ― 腕 2 本ずつで辺をつなぐ",
       lead="二つの炭素が正四面体の辺を一本共有すると、それぞれの腕が 2 本ずつ、辺の両端で出会う。これが二重結合。残った腕 2 本ずつに水素がつく。水素 4 個は一枚の平面に並び、エチレン C₂H₄ の平らな形になる。",
       scenes=[dict(cap="エチレン C₂H₄", s=double)],
       count="腕の勘定：炭素 2 個で 8 本 ＝ 炭素どうし 4 本（2 組）＋ 水素 4 本"),
  dict(file="bond3_triple.html", title="三重結合の腕",
       h1="三重結合 ― 腕 3 本ずつで面をつなぐ",
       lead="二つの炭素が正四面体の面を一枚共有すると、それぞれの腕が 3 本ずつ、面の三つの角で出会う。これが三重結合。残った腕は 1 本ずつで、水素は炭素と同じ一直線の上につく。アセチレン C₂H₂ のまっすぐな形になる。",
       scenes=[dict(cap="アセチレン C₂H₂", s=triple)],
       count="腕の勘定：炭素 2 個で 8 本 ＝ 炭素どうし 6 本（3 組）＋ 水素 2 本"),
  dict(file="bond4_co.html", title="CO の腕",
       h1="一酸化炭素 CO ― 面を共有し、両側に埋まった腕が 1 本ずつ残る",
       lead="酸素の正四面体にも腕が 4 本ある。水素がいないので、相手とつながなかった腕は自分の電子 2 個で埋まる。CO では炭素と酸素が面を一枚共有し、腕が 3 本ずつ出会う。炭素にも酸素にも、埋まった腕が 1 本ずつ残り、どちらも分子の軸の上を外向きに向く。",
       scenes=[dict(cap="一酸化炭素 CO（面を共有）", s=co)],
       count="腕の勘定：炭素と酸素で 8 本 ＝ つないだ腕 6 本（3 組）＋ 埋まった腕 2 本。炭素の腕が 1 本、埋まったまま外を向いている"),
  dict(file="bond5_co2.html", title="CO₂ の腕",
       h1="二酸化炭素 CO₂ ― 向かい合う二辺で、炭素の腕を使い切る",
       lead="CO₂ では炭素が正四面体の向かい合う二辺を、左右の酸素とそれぞれ共有する。腕は 2 本ずつ出会い、炭素の腕 4 本はすべてつながる。酸素には埋まった腕が 2 本ずつ残る。左右の酸素は炭素をはさんで一直線に並び、二つの辺は互いに直交する。",
       scenes=[dict(cap="二酸化炭素 CO₂（辺を二本共有）", s=co2)],
       count="腕の勘定：炭素 1 個と酸素 2 個で 12 本 ＝ つないだ腕 8 本（4 組）＋ 埋まった腕 4 本（どれも酸素の側）。炭素に空きも埋まった腕も残らない"),
]

TPL = pathlib.Path(__file__).with_name("bond_fig_template.html").read_text()
for p in PAGES:
    data = dict(scenes=[dict(cap=s["cap"], atoms=s["s"]["atoms"], tally=s["s"]["tally"], len2=s["s"]["len2"]) for s in p["scenes"]])
    html = (TPL.replace("%%TITLE%%", p["title"]).replace("%%H1%%", p["h1"]).replace("%%LEAD%%", p["lead"])
               .replace("%%COUNT%%", p["count"]).replace("%%DATA%%", json.dumps(data, ensure_ascii=False))
               .replace("%%NCOL%%", str(len(p["scenes"]))))
    pathlib.Path(__file__).with_name(p["file"]).write_text(html)
    print("wrote", p["file"])
