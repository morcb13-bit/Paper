"""
b13_check.py --- 平衡13進の桁列の検査

事前登録（判別法6。測る前に書く）

  検定B1 往復
      OK なら：整数と桁列の対応が可逆。表現が一意
      NG なら：変換か正規化が壊れている

  検定B2 桁の範囲
      OK なら：加算・減算・乗算のあとも全桁が -6..+6 に収まる
      NG なら：繰り上がりが足りない

  検定B3 繰り上がりの完全性（負の対照つき）
      OK なら：13 の繰り返し減算による繰り上がりと、商による繰り上がりが一致する。
               さらに 2026-02 版の ±1 固定繰り上がりは |合計| が大きい所で崩れる
      NG なら：どちらかが誤り。旧版が崩れなければ、この検査は NG を返せない検査

  検定B4 桁ずらし = 13 の冪倍
      OK なら：13^k 倍が桁の移動だけで書ける（乗算を要しない）
      NG なら：桁ずらしが 13 の冪倍になっていない

  検定B5 加減乗が整数演算と一致
      OK なら：桁の繰り返し加算と桁ずらしだけで乗算が作れている
      NG なら：どこかで乗算命令に相当する操作が要る

  検定B6 Z[phi]
      OK なら：phi^2 = phi + 1、sigma と N の乗法性、単数のノルムが ±1、
               符号判定が u^2 と 5v^2 の比較だけで正しい
      NG なら：整数対の算術が閉じていない
"""

import random
from b13num import (B13, Zphi, BASE, HALF, ZERO, ONE,
                    _carry_by_subtraction, _carry_fast, _normalize)

random.seed(13)
results = []


def report(name, ok, note=""):
    results.append((name, ok, note))
    print(f"{'OK ' if ok else 'NG '} {name}  {note}")


# ---------------- 検定B1 ----------------
vals = list(range(-2000, 2001))
vals += [random.randint(-13**12, 13**12) for _ in range(500)]
vals += [13**k for k in range(15)] + [-13**k for k in range(15)]
bad = [v for v in vals if B13.from_int(v).to_int() != v]
report("検定B1 往復", not bad, f"{len(vals)}件 / 不一致 {len(bad)}")

uniq_ok = True
for v in random.sample(vals, 200):
    x = B13.from_int(v)
    y = B13(list(x.c) + [0, 0])       # 末尾の0を足しても同じ桁列になるか
    if not x.eq(y):
        uniq_ok = False
report("検定B1b 表現の一意性", uniq_ok, "末尾0の付加で桁列が変わらない")


# ---------------- 検定B2 ----------------
rng_ok = True
for _ in range(500):
    x = B13.from_int(random.randint(-13**8, 13**8))
    y = B13.from_int(random.randint(-13**8, 13**8))
    for z in (x.add(y), x.sub(y), x.mul(y), x.shift(3), x.times_small(6)):
        if not z.digits_in_range():
            rng_ok = False
report("検定B2 桁の範囲", rng_ok, "加減乗・桁ずらし・桁倍のあと全桁 |c|<=6")


# ---------------- 検定B3 ----------------
same = all(_carry_by_subtraction(t) == _carry_fast(t) for t in range(-5000, 5001))
report("検定B3a 繰り上がり二法の一致", same, "t = -5000..5000 で完全一致")


def _normalize_old(digits):
    """2026-02 版の再現。繰り上がりを ±1 に固定している"""
    coeffs = list(digits)
    carry = 0
    for i in range(len(coeffs)):
        total = coeffs[i] + carry
        if total > HALF:
            coeffs[i] = total - BASE
            carry = 1
        elif total < -HALF:
            coeffs[i] = total + BASE
            carry = -1
        else:
            coeffs[i] = total
            carry = 0
    if carry != 0:
        coeffs.append(carry)
    while len(coeffs) > 1 and coeffs[-1] == 0:
        coeffs.pop()
    return coeffs


def digits_to_int(ds):
    v = 0
    for d in reversed(ds):
        v = v * BASE + d
    return v


old_out_of_range = 0
old_value_wrong = 0
for _ in range(400):
    raw = [random.randint(-40, 40) for _ in range(5)]
    truth = digits_to_int(raw)
    got = _normalize_old(raw)
    if digits_to_int(got) != truth:
        old_value_wrong += 1
    if any(d < -HALF or d > HALF for d in got):
        old_out_of_range += 1
    assert digits_to_int(_normalize(raw)) == truth
    assert all(-HALF <= d <= HALF for d in _normalize(raw))
report("検定B3b 負の対照（旧版・値）", old_value_wrong == 0,
       f"値は保たれる（400件中 崩れ {old_value_wrong}）── 旧版の欠陥は値ではない")
report("検定B3b2 負の対照（旧版・桁の範囲）", old_out_of_range > 0,
       f"400件中 {old_out_of_range} 件で桁が -6..+6 を外れる")

# 同じ値に二つの桁列が付く＝正準でない
x_raw = [20, 0]
y_raw = [7, 1]
canon_broken = (digits_to_int(x_raw) == digits_to_int(y_raw)
                and _normalize_old(x_raw) != _normalize_old(y_raw)
                and _normalize(x_raw) == _normalize(y_raw))
report("検定B3c 正準性", canon_broken,
       f"旧版 {_normalize_old(x_raw)} と {_normalize_old(y_raw)} / 新版 {_normalize(x_raw)}")


# ---------------- 検定B4 ----------------
shift_ok = True
for _ in range(200):
    v = random.randint(-13**6, 13**6)
    x = B13.from_int(v)
    for k in range(1, 5):
        if x.shift(k).to_int() != v * BASE**k:
            shift_ok = False
report("検定B4 桁ずらし = 13^k 倍", shift_ok, "k = 1..4")


# ---------------- 検定B5 ----------------
ops_ok = True
for _ in range(400):
    u = random.randint(-13**7, 13**7)
    v = random.randint(-13**7, 13**7)
    x, y = B13.from_int(u), B13.from_int(v)
    if x.add(y).to_int() != u + v:
        ops_ok = False
    if x.sub(y).to_int() != u - v:
        ops_ok = False
    if x.mul(y).to_int() != u * v:
        ops_ok = False
report("検定B5 加減乗の一致", ops_ok, "乗算は桁倍（繰り返し加算）と桁ずらしの和")


# ---------------- 検定B6 ----------------
phi = Zphi.from_ints(0, 1)
one = Zphi.from_ints(1, 0)
ok_phi2 = phi.mul(phi).eq(phi.add(one))
report("検定B6a phi^2 = phi + 1", ok_phi2)

mult_ok = True
norm_ok = True
for _ in range(300):
    x = Zphi.from_ints(random.randint(-200, 200), random.randint(-200, 200))
    y = Zphi.from_ints(random.randint(-200, 200), random.randint(-200, 200))
    if not x.mul(y).conj().eq(x.conj().mul(y.conj())):
        mult_ok = False
    if x.mul(y).norm().to_int() != x.norm().to_int() * y.norm().to_int():
        norm_ok = False
report("検定B6b sigma の乗法性", mult_ok)
report("検定B6c N の乗法性", norm_ok)

u = one
unit_ok = True
for k in range(1, 25):
    u = u.mul(phi)
    if abs(u.norm().to_int()) != 1:
        unit_ok = False
report("検定B6d 単数 phi^k のノルム", unit_ok, "k = 1..24 で |N| = 1")

sign_ok = True
PHI = (1 + 5 ** 0.5) / 2
for _ in range(2000):
    a = random.randint(-300, 300)
    b = random.randint(-300, 300)
    z = Zphi.from_ints(a, b)
    approx = a + b * PHI
    want = 0 if abs(approx) < 1e-9 else (1 if approx > 0 else -1)
    if z.sign() != want:
        sign_ok = False
report("検定B6e 符号判定", sign_ok, "u^2 と 5v^2 の比較のみ。浮動小数は照合側だけ")

print()
n_ok = sum(1 for _, ok, _ in results if ok)
print(f"{n_ok}/{len(results)} OK")
