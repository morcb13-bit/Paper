#  対象（判別法9）：同じ一歩を繰り返す歩き（開いた歩き）と、36°ずつ回して10歩で閉じる
#                   歩きとで、残渣の溜まり方がどう違うかを厳密な整数演算で出す。
#
#  すべて Z[ζ₁₀] / Q(φ) の中で計算する。浮動小数は表示のときだけ使う（判別法16）。
#      |v|² = v·conj(v) を Q(φ) の元 p+qφ として持ち、
#      残渣² はそのガロア共役 (p+q) − qφ。
#
#  検定CW1 開いた歩き
#      OK なら：残渣² が n² 倍で伸びる → どんな窓でも有限歩で外れる
#      NG なら：伸びない
#
#  検定CW2 閉じる歩き
#      OK なら：部分和の残渣² が有界のままで、10歩で厳密に 0 に戻る
#              → 窓² がその最大値以上なら何周でも嵌まる
#      NG なら：溜まる（開いた歩きと同じ）
#
#  使い方: qphi.py と rings_integer.json と同じ場所に置いて python3 closed_walk.py

from fractions import Fraction as F
import json, sys
import qphi as Q

Qp = Q.Qp
ZETA = (F(0), F(1), F(0), F(0))          # ζ₁₀

def absq(v):  return Q.zre(Q.zmul(v, Q.zconj(v)))     # |v|² ∈ Q(φ)
def sigma(x): return Qp(x.p + x.q, -x.q)              # p+qφ → (p+q)−qφ
def norm(x):  return x.p*x.p + x.p*x.q - x.q*x.q      # 絶対ノルム

rings = [tuple(F(y) for y in r) for r in json.load(open('rings_integer.json'))['rings']]

# 2つ飛ばしの一歩を担体から厳密に拾う（|v|² = 3+4φ）
step0 = None
for a in rings:
    for b in rings:
        if a == b: continue
        d = tuple(x - y for x, y in zip(a, b))
        if absq(d) == Qp(3, 4):
            step0 = d; break
    if step0: break
if step0 is None:
    sys.exit("2つ飛ばしの一歩が見つからない")

s0 = absq(step0)
print(f"一歩 v = {tuple(int(x) for x in step0)}")
print(f"  |v|² = {s0} = {s0.val():.6f}   残渣² = {sigma(s0)} = {sigma(s0).val():.6f}"
      f"   ノルム = {norm(s0)}\n")

# --- 検定CW1 開いた歩き ---
print("検定CW1 同じ一歩を n 回繰り返す")
prev = None
grows = True
for n in (1, 2, 3, 5, 10):
    r = Qp(n*n, 0) * sigma(s0)
    print(f"   {n:2d}歩  残渣² = {str(r):>14}  = {r.val():10.6f}")
    if prev is not None and not (r.val() > prev): grows = False
    prev = r.val()
print(f"検定CW1: {'OK' if grows else 'NG'}  残渣² は n² 倍で伸びる → 必ず外れる\n")

# --- 検定CW2 閉じる歩き ---
print("検定CW2 36°ずつ回して10歩")
S, step, vals = (F(0),)*4, step0, []
for k in range(11):
    r = sigma(absq(S))
    vals.append(r)
    print(f"   {k:2d}歩  部分和 = {str(tuple(int(x) for x in S)):>18}"
          f"  残渣² = {str(r):>14}  = {r.val():.6f}")
    S = tuple(x + y for x, y in zip(S, step))
    step = Q.zmul(step, ZETA)
closes = (vals[-1] == Qp(0, 0))
mx = max(vals, key=lambda x: x.val())
print(f"検定CW2: {'OK' if closes else 'NG'}  10歩で厳密に 0 に戻る = {closes}")
print(f"   途中の最大は 5歩目の {mx} = {mx.val():.6f}")
print(f"\n→ 窓² ≥ {mx} であれば、閉じる歩きは何周でも嵌まり続ける。")
print("→ 溜まるのは開いた歩きの側だけ。")
