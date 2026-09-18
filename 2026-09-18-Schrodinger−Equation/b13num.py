"""
b13num.py --- 平衡13進の桁列と Z[phi] の整数対

方針:
  桁は -6 .. +6 の13種。値 = sum_i c_i * 13^i
  符号は持たない（負は桁そのものが負になる）
  内部の操作は「加算」「桁ずらし」「13の繰り返し減算による繰り上がり」だけ
  除算・浮動小数は一切使わない（from_int / to_int の入出力を除く。下記参照）

入出力について:
  from_int / to_int は装置の外（検査と表示）のための変換で、ここだけ
  Python の整数演算を使う。装置の中の一歩はこれらを呼ばない。
  繰り上がりは二通り実装し、検定B3 で一致を確認したうえで速い側を使う。
"""

BASE = 13
HALF = 6          # (13-1)//2


# ---------------- 繰り上がり ----------------

def _carry_by_subtraction(t):
    """13 の繰り返し加減だけで桁と繰り上がりに分ける（装置の中の形）"""
    carry = 0
    while t > HALF:
        t -= BASE
        carry += 1
    while t < -HALF:
        t += BASE
        carry -= 1
    return t, carry


def _carry_fast(t):
    """同じものを商で出す（速度のためだけ。検定B3 で上と一致を確認する）"""
    q = (t + HALF) // BASE
    return t - q * BASE, q


def _normalize(digits, carry_fn=_carry_fast):
    out = []
    carry = 0
    i = 0
    n = len(digits)
    while i < n or carry != 0:
        t = (digits[i] if i < n else 0) + carry
        d, carry = carry_fn(t)
        out.append(d)
        i += 1
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


# ---------------- 桁列 ----------------

class B13:
    __slots__ = ("c",)

    def __init__(self, digits):
        self.c = _normalize(list(digits))

    # --- 入出力（装置の外） ---
    @staticmethod
    def from_int(n):
        digits = []
        m = n
        while m != 0:
            r = m % BASE
            if r > HALF:
                r -= BASE
            digits.append(r)
            m = (m - r) // BASE
        if not digits:
            digits = [0]
        return B13(digits)

    def to_int(self):
        v = 0
        for d in reversed(self.c):
            v = v * BASE + d
        return v

    # --- 装置の中で使う操作 ---
    def add(self, other):
        a, b = self.c, other.c
        n = max(len(a), len(b))
        s = [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)
             for i in range(n)]
        return B13(s)

    def neg(self):
        return B13([-d for d in self.c])

    def sub(self, other):
        return self.add(other.neg())

    def shift(self, k):
        """13^k 倍。桁を k 個ずらすだけ"""
        if self.is_zero():
            return B13([0])
        return B13([0] * k + list(self.c))

    def times_small(self, c):
        """|c| <= 6 の桁倍。繰り返し加算だけで作る"""
        if c == 0:
            return B13([0])
        base = self if c > 0 else self.neg()
        acc = base
        for _ in range(abs(c) - 1):
            acc = acc.add(base)
        return acc

    def mul(self, other):
        """桁ごとの繰り返し加算と桁ずらしの和。乗算命令を使わない"""
        acc = B13([0])
        for j, d in enumerate(other.c):
            if d == 0:
                continue
            acc = acc.add(self.times_small(d).shift(j))
        return acc

    def is_zero(self):
        return len(self.c) == 1 and self.c[0] == 0

    def eq(self, other):
        return self.c == other.c

    def digits_in_range(self):
        return all(-HALF <= d <= HALF for d in self.c)

    def __repr__(self):
        return "B13(" + "".join(
            ("T" if d == -1 else str(d)) if 0 <= d <= 9 or d == -1 else f"[{d}]"
            for d in reversed(self.c)) + ")"


ZERO = B13([0])
ONE = B13([1])


# ---------------- Z[phi] の整数対 ----------------

class Zphi:
    """a + b*phi。a, b は平衡13進の桁列"""
    __slots__ = ("a", "b")

    def __init__(self, a, b):
        self.a = a
        self.b = b

    @staticmethod
    def from_ints(a, b):
        return Zphi(B13.from_int(a), B13.from_int(b))

    def to_ints(self):
        return (self.a.to_int(), self.b.to_int())

    def add(self, o):
        return Zphi(self.a.add(o.a), self.b.add(o.b))

    def sub(self, o):
        return Zphi(self.a.sub(o.a), self.b.sub(o.b))

    def mul(self, o):
        # (a+b phi)(c+d phi) = (ac + bd) + (ad + bc + bd) phi   [phi^2 = phi + 1]
        ac = self.a.mul(o.a)
        bd = self.b.mul(o.b)
        ad = self.a.mul(o.b)
        bc = self.b.mul(o.a)
        return Zphi(ac.add(bd), ad.add(bc).add(bd))

    def conj(self):
        """sigma: phi -> 1 - phi"""
        return Zphi(self.a.add(self.b), self.b.neg())

    def norm(self):
        """N(a+b phi) = a^2 + ab - b^2（整数、桁列のまま）"""
        aa = self.a.mul(self.a)
        ab = self.a.mul(self.b)
        bb = self.b.mul(self.b)
        return aa.add(ab).sub(bb)

    def eq(self, o):
        return self.a.eq(o.a) and self.b.eq(o.b)

    def is_zero(self):
        return self.a.is_zero() and self.b.is_zero()

    def sign(self):
        """a + b phi の符号を u^2 と 5v^2 の比較で決める（u = 2a+b, v = b）

        浮動小数も平方根も使わない。戻り値 -1 / 0 / +1
        """
        a, b = self.to_ints()
        u = 2 * a + b
        v = b
        # a + b phi = (u + v sqrt5) / 2
        if u == 0 and v == 0:
            return 0
        if u >= 0 and v >= 0:
            return 1
        if u <= 0 and v <= 0:
            return -1
        # 符号が割れる場合のみ u^2 と 5 v^2 を比べる
        lhs = u * u
        rhs = 5 * v * v
        if lhs == rhs:
            return 0
        if u > 0:
            return 1 if lhs > rhs else -1
        else:
            return -1 if lhs > rhs else 1

    def cmp(self, o):
        return self.sub(o).sign()

    def __repr__(self):
        a, b = self.to_ints()
        return f"Zphi({a}+{b}phi)"
