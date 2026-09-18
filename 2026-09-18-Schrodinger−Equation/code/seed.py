"""
seed.py --- 起点が目盛りを決めるか

対象を一行で書く（判別法9）：
  いま見ているのは「二階漸化式 s_{n+2} = lambda s_{n+1} - L^2 s_n の起点 (s_0, s_1)」。
  起点が二次形式 q の値を決め、それが以後の歩の単位になるかどうか。

  q(a, b) = b^2 - lambda a b + L^2 a^2
  q0 = q(s_0, s_1)

事前登録（判別法6）

  検定M1  起点が単位になるか
      OK なら：q(s_n, s_{n+1}) = L^(2n) * q0 が、種を変えても固有値を変えても厳密成立。
               q0 は起点だけで決まり、以後の歩はすべてその倍数として測られる
      NG なら：q0 が起点以外にも依存する

  検定M2  種 (1,1) と種 (1,0)
      OK なら：種 (1,1) では q0 = 10 - lambda となり、
               lambda = 6, 2, 1, -2, -3 に対して 4, 8, 9, 12, 13 の五つに分かれる。
               負の対照として種 (1,0) では五つとも 9 = L^2 に潰れる
      NG なら：分かれない、あるいは値が違う

  検定M3  足せないこと
      異なる q0 を持つ二つの成分（別の固有値に属する）を足したものは、
      一つの lambda に対する q の等式を満たさないはず
      OK なら：和は L^2 倍の形を保たない。単位の違うものは足せない
      NG なら：足しても形が保たれる。足せないという読みが誤り

  検定M4  同じ固有値どうしは足せること
      OK なら：同じ lambda に属する二つの数列の和は、やはり L^2 倍の形を保つ
               （線形性。ただし q0 は和の起点から決まり直す）
      NG なら：同じ固有値でも足せない

  検定M5  L がノードの中で一定であること（速さの側）
      OK なら：q の倍率は起点にも成分にも依らず L^2 の一つだけ
      NG なら：倍率が起点によって変わる
"""

import random

random.seed(13)
results = []


def report(name, ok, note=""):
    results.append((name, ok, note))
    print(f"{'OK ' if ok else 'NG '} {name}  {note}")


L = 3
EIGS = [6, 2, 1, -2, -3]


def q(a, b, lam):
    return b * b - lam * a * b + L * L * a * a


def run(s0, s1, lam, steps):
    seq = [s0, s1]
    for _ in range(steps):
        seq.append(lam * seq[-1] - L * L * seq[-2])
    return seq


# ---------------- 検定M1 ----------------
m1_ok = True
checked = 0
for lam in EIGS:
    seeds = [(1, 1), (1, 0), (0, 1), (2, -3), (5, 7), (-4, 1)]
    for s0, s1 in seeds:
        seq = run(s0, s1, lam, 12)
        q0 = q(seq[0], seq[1], lam)
        for n in range(len(seq) - 1):
            checked += 1
            if q(seq[n], seq[n + 1], lam) != (L * L) ** n * q0:
                m1_ok = False
report("検定M1 起点が単位になる", m1_ok,
       f"q(s_n, s_{{n+1}}) = L^(2n)·q0 を {checked} 件で厳密確認（5固有値 × 6種 × 12歩）")


# ---------------- 検定M2 ----------------
q11 = {lam: q(1, 1, lam) for lam in EIGS}
q10 = {lam: q(1, 0, lam) for lam in EIGS}
print()
print("  種ごとの q0")
print("    lambda    種 (1,1)    種 (1,0)")
for lam in EIGS:
    print(f"      {lam:>3}         {q11[lam]:>3}         {q10[lam]:>3}")
m2_ok = (all(q11[lam] == 10 - lam for lam in EIGS)
         and sorted(q11.values()) == [4, 8, 9, 12, 13]
         and set(q10.values()) == {9})
report("検定M2 種 (1,1) と種 (1,0)", m2_ok,
       f"(1,1) → {sorted(q11.values())} の五つ / (1,0) → 五つとも {list(set(q10.values()))[0]} に潰れる")


# ---------------- 検定M3 ----------------
# 別の固有値に属する二つの数列を足す。和に対して単一の lambda で q の等式が立つか
m3_rows = []
m3_ok = True
for a_lam, b_lam in [(6, 2), (2, 1), (1, -3), (-2, -3), (6, -3)]:
    sa = run(1, 1, a_lam, 10)
    sb = run(1, 1, b_lam, 10)
    tot = [x + y for x, y in zip(sa, sb)]
    holds = []
    for lam in (a_lam, b_lam):
        q0 = q(tot[0], tot[1], lam)
        ok = all(q(tot[n], tot[n + 1], lam) == (L * L) ** n * q0
                 for n in range(len(tot) - 1))
        holds.append(ok)
    broke = not any(holds)
    m3_rows.append((a_lam, b_lam, holds, broke))
    if not broke:
        m3_ok = False
print()
print("  異なる固有値の成分を足したとき、単一の lambda で q が保たれるか")
for a_lam, b_lam, holds, broke in m3_rows:
    print(f"    λ={a_lam:>2} と λ={b_lam:>2} の和   λ={a_lam} で {holds[0]} / λ={b_lam} で {holds[1]}   "
          f"{'足せない' if broke else '足せてしまう'}")
report("検定M3 足せないこと", m3_ok,
       "単位の違う成分の和は、どちらの lambda でも L^2 倍の形を保たない")


# ---------------- 検定M4 ----------------
m4_ok = True
for lam in EIGS:
    for _ in range(20):
        s = run(random.randint(-6, 6), random.randint(-6, 6), lam, 10)
        t = run(random.randint(-6, 6), random.randint(-6, 6), lam, 10)
        tot = [x + y for x, y in zip(s, t)]
        q0 = q(tot[0], tot[1], lam)
        if not all(q(tot[n], tot[n + 1], lam) == (L * L) ** n * q0
                   for n in range(len(tot) - 1)):
            m4_ok = False
report("検定M4 同じ固有値どうしは足せる", m4_ok,
       "和も L^2 倍の形を保つ。q0 は和の起点から決まり直す")


# ---------------- 検定M5 ----------------
ratios = set()
for lam in EIGS:
    for s0, s1 in [(1, 1), (1, 0), (3, -2), (7, 4)]:
        seq = run(s0, s1, lam, 8)
        q0 = q(seq[0], seq[1], lam)
        if q0 == 0:
            continue
        for n in range(len(seq) - 2):
            num = q(seq[n + 1], seq[n + 2], lam)
            den = q(seq[n], seq[n + 1], lam)
            if den != 0 and num % den == 0:
                ratios.add(num // den)
m5_ok = (ratios == {L * L})
report("検定M5 倍率は L^2 の一つだけ", m5_ok,
       f"起点・固有値を変えても倍率は {sorted(ratios)}")


# ---------------- 参考 ----------------
print()
print("  種 (1,1) から 8歩（λ ごと）。q は L^2 = 9 倍ずつ")
for lam in EIGS:
    seq = run(1, 1, lam, 8)
    print(f"    λ={lam:>3}  q0={q(seq[0], seq[1], lam):>3}  {seq[:7]}")

print()
n_ok = sum(1 for _, ok, _ in results if ok)
print(f"{n_ok}/{len(results)} OK")
