#!/usr/bin/env python3
# 600胞体を Z[φ] の整数演算だけで組む。浮動小数は画面に落とすときだけ。
#
#  検定1 ノルム
#      OK なら：120点すべてが厳密に単位ノルム（2倍座標で 4）。座標が Z[φ] に乗っている
#      NG なら：頂点表が誤り。以降の数値は全部無効
#
#  検定2 一歩の隔たり
#      OK なら：最近接の内積が厳密に φ/2 の一種類だけ。辺が 720 本、各点 12 本
#      NG なら：辺の定義が浮動小数の閾値に依存している
#
#  検定3 掛け算で閉じる（四元数）
#      OK なら：120点は掛け算で閉じた群。回転が整数演算のまま書ける
#      NG なら：この120点は群ではなく、回転に浮動小数が要る
#
#  検定4 面と胞
#      OK なら：三角面 1200・四面体胞 600。名前どおりの立体
#      NG なら：隣接の取り方が誤り
#
#  検定5 w での層分け
#      OK なら：第4成分が取る値は5種（±)で、層の個数が 1,12,20,12,30,12,20,12,1
#      NG なら：層構造の読みが誤り
#
#  検定6 負の対照（偶置換を全置換にする）
#      NG が出るのが正しい。OK で終わったら検定が反証不能

from itertools import permutations
from fractions import Fraction
import json, math

# ---- Z[φ]：(a,b) は a+bφ、φ²=φ+1。係数は Fraction（2で割る操作が入るため）----
def add(x, y): return (x[0]+y[0], x[1]+y[1])
def sub(x, y): return (x[0]-y[0], x[1]-y[1])
def neg(x):    return (-x[0], -x[1])
def mul(x, y):
    a, b = x; c, d = y
    # (a+bφ)(c+dφ) = ac + (ad+bc)φ + bd φ² = (ac+bd) + (ad+bc+bd)φ
    return (a*c + b*d, a*d + b*c + b*d)
Z = (Fraction(0), Fraction(0))
ONE = (Fraction(1), Fraction(0))
PHI = (Fraction(0), Fraction(1))
HALF = Fraction(1, 2)
def sc(x, k): return (x[0]*k, x[1]*k)
def val(x): return float(x[0]) + float(x[1])*(1+5**0.5)/2

INV_PHI = (Fraction(-1), Fraction(1))          # 1/φ = φ-1
TWO = (Fraction(2), Fraction(0))

def even_perms(t):
    """偶置換のみ"""
    out = []
    base = list(range(4))
    for p in permutations(base):
        # 置換の符号
        s, seen = 1, [False]*4
        for i in range(4):
            if seen[i]: continue
            j, ln = i, 0
            while not seen[j]:
                seen[j] = True; j = p[j]; ln += 1
            if ln % 2 == 0: s = -s
        if s == 1:
            out.append(tuple(t[p[i]] for i in range(4)))
    return out

def build(all_perms=False):
    V = set()
    # (±1,0,0,0) 型 … 8
    for i in range(4):
        for s in (1, -1):
            v = [Z]*4; v[i] = sc(ONE, s); V.add(tuple(v))
    # (±1/2,±1/2,±1/2,±1/2) … 16
    for m in range(16):
        v = tuple(sc(ONE, HALF*(1 if (m >> k) & 1 else -1)) for k in range(4))
        V.add(v)
    # (±φ, ±1, ±1/φ, 0)/2 の偶置換 … 96
    base = (sc(PHI, HALF), sc(ONE, HALF), sc(INV_PHI, HALF), Z)
    perms = list(permutations(range(4))) if all_perms else None
    if all_perms:
        cand = [tuple(base[p[i]] for i in range(4)) for p in perms]
    else:
        cand = even_perms(base)
    for t in cand:
        for m in range(16):
            v = tuple(neg(t[k]) if (m >> k) & 1 else t[k] for k in range(4))
            V.add(v)
    return sorted(V, key=lambda v: tuple((float(c[0]), float(c[1])) for c in v))

def dot(u, v):
    r = Z
    for k in range(4): r = add(r, mul(u[k], v[k]))
    return r

def qmul(u, v):
    a1,b1,c1,d1 = u; a2,b2,c2,d2 = v
    return (
        sub(sub(sub(mul(a1,a2), mul(b1,b2)), mul(c1,c2)), mul(d1,d2)),
        sub(add(add(mul(a1,b2), mul(b1,a2)), mul(c1,d2)), mul(d1,c2)),
        add(add(sub(mul(a1,c2), mul(b1,d2)), mul(c1,a2)), mul(d1,b2)),
        add(sub(add(mul(a1,d2), mul(b1,c2)), mul(c1,b2)), mul(d1,a2)),
    )

def run(all_perms=False, quiet=False):
    V = build(all_perms)
    res = {}
    say = (lambda *a: None) if quiet else print

    res['n'] = len(V)
    say(f"頂点 {len(V)} 個")

    # 検定1
    ok1 = all(dot(v, v) == ONE for v in V)
    res[1] = ok1
    say(f"検定1 ノルム : {'OK' if ok1 else 'NG'}")

    # 検定2
    idx = {v: i for i, v in enumerate(V)}
    target = sc(PHI, HALF)          # φ/2
    edges, deg = [], [0]*len(V)
    for i in range(len(V)):
        for j in range(i+1, len(V)):
            if dot(V[i], V[j]) == target:
                edges.append((i, j)); deg[i] += 1; deg[j] += 1
    ok2 = (len(edges) == 720) and all(d == 12 for d in deg)
    res[2] = ok2; res['edges'] = len(edges); res['deg'] = sorted(set(deg))
    say(f"検定2 一歩    : {'OK' if ok2 else 'NG'}  辺 {len(edges)} / 次数 {sorted(set(deg))}")

    # 検定3
    S = set(V)
    ok3 = all(qmul(V[i], V[j]) in S for i in range(0, len(V), 7) for j in range(len(V)))
    res[3] = ok3
    say(f"検定3 掛け算  : {'OK' if ok3 else 'NG'}")

    # 検定4
    adj = [set() for _ in V]
    for i, j in edges: adj[i].add(j); adj[j].add(i)
    tri = 0; tet = 0
    for i, j in edges:
        if i > j: continue
        common = sorted(adj[i] & adj[j])
        tri += sum(1 for k in common if k > j)
        for a in range(len(common)):
            for b in range(a+1, len(common)):
                k, l = common[a], common[b]
                if k > j and l > j and l in adj[k]: tet += 1
    # {a<b<c<d} が数えられるのは辺(a,b)のときだけ（他の辺では k>j が破れる）ので重複なし
    ok4 = (tri == 1200 and tet == 600)
    res[4] = ok4; res['tri'] = tri; res['tet'] = tet
    say(f"検定4 面と胞  : {'OK' if ok4 else 'NG'}  面 {tri} / 胞 {tet}")

    # 検定5
    shells = {}
    for v in V: shells.setdefault(v[3], []).append(v)
    key = sorted(shells, key=lambda x: val(x))
    counts = [len(shells[k]) for k in key]
    ok5 = (counts == [1, 12, 20, 12, 30, 12, 20, 12, 1])
    res[5] = ok5; res['shell_counts'] = counts
    res['shell_w'] = [f"{k[0]}+{k[1]}φ" for k in key]
    say(f"検定5 w の層  : {'OK' if ok5 else 'NG'}  {counts}")

    return V, edges, res

if __name__ == "__main__":
    print("=== 本番（偶置換） ===")
    V, E, r = run()
    print()
    print("=== 検定6 負の対照（全置換にする） ===")
    _, _, r2 = run(all_perms=True)
    ok6 = not all(r2.get(k, False) for k in (1, 2, 3, 4, 5))
    print(f"検定6 負の対照: {'OK（NGが出た＝反証可能）' if ok6 else 'NG（検定が反証不能）'}")
    print()
    print("w の層の値 :", r['shell_w'])
    print("層の個数   :", r['shell_counts'])

    out = {
        "v": [[round(val(c), 12) for c in v] for v in V],
        "ve": [[[str(c[0]), str(c[1])] for c in v] for v in V],
        "e": E,
        "w": [f"{v[3][0]}+{v[3][1]}φ" for v in V],
    }
    with open("cell600.json", "w") as f: json.dump(out, f, separators=(",", ":"))
    print("→ cell600.json", len(out["v"]), "点 /", len(E), "辺")
