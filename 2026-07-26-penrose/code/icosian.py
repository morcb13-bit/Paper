#  icosian.py  単位イコシアン120個（＝600胞体の頂点＝二項20面体群）
#
#      係数は Q(√5) を Fraction で持つ。浮動小数はどこにも入らない（val() は表示用）。
#      φ² = φ+1 だけで積が閉じる。
#
#      構成（Conway–Smith『On Quaternions and Octonions』ch.8 と同じ並び）
#        ±1, ±i, ±j, ±k                        8個
#        (±1±i±j±k)/2                          16個
#        (0, ±1, ±1/φ, ±φ)/2 の偶置換          96個
#                                             ─────
#                                              120個
#
#      検定（期待値は引継書 v223 の記録）
#        1 120個ちょうど・ノルムが1でないもの0（厳密等号）
#        2 14400通りの積がすべて集合の中（四元数の積で閉じた群）
#        3 位数の分布 {1:1, 2:1, 3:20, 4:30, 5:24, 6:20, 10:24}
#        4 各頂点の最近接は内積 φ/2 でちょうど12個・辺720本
#        5 ab≠ba の対が 13320/14400（担体に非可換が入っている）

from fractions import Fraction as F
from itertools import permutations


class Q5:
    """a + b·φ（φ² = φ + 1）。係数は Fraction。"""
    __slots__ = ("a", "b")

    def __init__(self, a=0, b=0):
        self.a = a if isinstance(a, F) else F(a)
        self.b = b if isinstance(b, F) else F(b)

    def __add__(self, o):
        o = _q5(o)
        return Q5(self.a + o.a, self.b + o.b)

    __radd__ = __add__

    def __sub__(self, o):
        o = _q5(o)
        return Q5(self.a - o.a, self.b - o.b)

    def __rsub__(self, o):
        return _q5(o) - self

    def __neg__(self):
        return Q5(-self.a, -self.b)

    def __mul__(self, o):
        o = _q5(o)
        # (a+bφ)(c+dφ) = ac + (ad+bc)φ + bd·φ² = (ac+bd) + (ad+bc+bd)φ
        return Q5(self.a * o.a + self.b * o.b,
                  self.a * o.b + self.b * o.a + self.b * o.b)

    __rmul__ = __mul__

    def __eq__(self, o):
        if not isinstance(o, Q5):
            try:
                o = _q5(o)
            except TypeError:
                return NotImplemented
        return self.a == o.a and self.b == o.b

    def __hash__(self):
        return hash((self.a, self.b))

    def val(self):
        """表示・作図用の実数値。判定には使わない。"""
        return float(self.a) + float(self.b) * 1.618033988749894848204586834

    def __repr__(self):
        return "Q5(%s,%s)" % (self.a, self.b)


def _q5(x):
    if isinstance(x, Q5):
        return x
    if isinstance(x, (int, F)):
        return Q5(x, 0)
    raise TypeError("Q5 に混ぜられない型: %r" % type(x))


ZERO = Q5(0, 0)
ONE = Q5(1, 0)
PHI = Q5(0, 1)
INV_PHI = Q5(-1, 1)          # 1/φ = φ − 1


class H:
    """四元数。成分は Q5 の4つ。"""
    __slots__ = ("c",)

    def __init__(self, c):
        self.c = tuple(_q5(x) for x in c)

    def __mul__(self, o):
        a0, a1, a2, a3 = self.c
        b0, b1, b2, b3 = o.c
        return H([a0 * b0 - a1 * b1 - a2 * b2 - a3 * b3,
                  a0 * b1 + a1 * b0 + a2 * b3 - a3 * b2,
                  a0 * b2 - a1 * b3 + a2 * b0 + a3 * b1,
                  a0 * b3 + a1 * b2 - a2 * b1 + a3 * b0])

    def conj(self):
        return H([self.c[0], -self.c[1], -self.c[2], -self.c[3]])

    def norm(self):
        """q q̄ の実部。単位イコシアンなら ONE。"""
        return (self.c[0] * self.c[0] + self.c[1] * self.c[1]
                + self.c[2] * self.c[2] + self.c[3] * self.c[3])

    def __eq__(self, o):
        return isinstance(o, H) and self.c == o.c

    def __hash__(self):
        return hash(self.c)

    def __repr__(self):
        return "H%s" % (self.c,)


def build():
    """単位イコシアン120個。順序は 8 → 16 → 96。"""
    out = []
    # ±1, ±i, ±j, ±k
    for k in range(4):
        for s in (1, -1):
            c = [ZERO, ZERO, ZERO, ZERO]
            c[k] = Q5(s, 0)
            out.append(H(c))
    # (±1±i±j±k)/2
    half = Q5(F(1, 2), 0)
    for s0 in (1, -1):
        for s1 in (1, -1):
            for s2 in (1, -1):
                for s3 in (1, -1):
                    out.append(H([half * s0, half * s1, half * s2, half * s3]))
    # (0, ±1, ±1/φ, ±φ)/2 の偶置換
    base = [ZERO, Q5(F(1, 2), 0), Q5(F(-1, 2), F(1, 2)), Q5(0, F(1, 2))]
    even = [p for p in permutations(range(4)) if _parity(p) == 0]
    seen = set()
    for p in even:
        for s1 in (1, -1):
            for s2 in (1, -1):
                for s3 in (1, -1):
                    sg = [1, s1, s2, s3]
                    c = [ZERO] * 4
                    for i in range(4):
                        c[p[i]] = base[i] * sg[i]
                    q = H(c)
                    if q not in seen:
                        seen.add(q)
                        out.append(q)
    return out


def _parity(p):
    n = 0
    p = list(p)
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]:
                n += 1
    return n % 2


def _order(q, ident):
    r = q
    for k in range(1, 41):
        if r == ident:
            return k
        r = r * q
    return None


def main():
    from collections import Counter
    V = build()
    ident = H([ONE, ZERO, ZERO, ZERO])
    S = set(V)

    print("== 検定1 120個・ノルムが厳密に1 ==")
    bad = sum(1 for q in V if q.norm() != ONE)
    print("   個数 %d ／ ノルムが1でないもの %d → %s"
          % (len(V), bad, "OK" if len(V) == 120 and bad == 0 else "NG"))

    print("== 検定2 積で閉じているか ==")
    out = sum(1 for a in V for b in V if (a * b) not in S)
    print("   %d 通りのうち集合の外に出たもの %d → %s"
          % (len(V) ** 2, out, "OK" if out == 0 else "NG"))

    print("== 検定3 位数の分布 ==")
    dist = Counter(_order(q, ident) for q in V)
    want = {1: 1, 2: 1, 3: 20, 4: 30, 5: 24, 6: 20, 10: 24}
    print("   %s → %s" % (dict(sorted(dist.items())),
                          "OK" if dict(dist) == want else "NG"))

    print("== 検定4 最近接12個・辺720本 ==")
    half_phi = Q5(0, F(1, 2))
    deg = []
    for a in V:
        n = 0
        for b in V:
            ip = (a.c[0] * b.c[0] + a.c[1] * b.c[1]
                  + a.c[2] * b.c[2] + a.c[3] * b.c[3])
            if ip == half_phi:
                n += 1
        deg.append(n)
    edges = sum(deg) // 2
    print("   最近接（内積 φ/2）の個数 %s ／ 辺 %d 本 → %s"
          % (set(deg), edges,
             "OK" if set(deg) == {12} and edges == 720 else "NG"))

    print("== 検定5 非可換 ==")
    nc = sum(1 for a in V for b in V if (a * b) != (b * a))
    print("   ab≠ba の対 %d / %d → %s"
          % (nc, len(V) ** 2, "OK" if nc == 13320 else "NG"))


if __name__ == "__main__":
    main()
