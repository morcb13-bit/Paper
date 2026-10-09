# 立体の検定（整数だけ）
#  対象：正四面体の中心に置いた原子と、4頂点に置いた4種の置換基 A,B,C,D。
#        頂点は立方体の一つおきの角 (1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1)。
#  鏡：平面 x=y に映す（x と y を入れ替える）。右型の B と C の席が入れ替わったものが左型。
#  検定S1 正四面体を自分に重ねる回転（整数の行列）はちょうど12通り
#  検定S2 左型を12通りのどれで回しても、右型と一致する置換基は最大2個（4個にならない）
#      NG なら：右と左は回せば重なる（作り分ける意味がない）
#  検定S3 負の対照：右型を12通りで回すと、4個一致する回し方が1通りだけある
#  検定S4 型紙（B,C,D の3つの受け口）に、右型は3つともはまり、左型はどの回し方でも2つまで
#  検定S5 負の対照：受け口を2つ（B,C）にした型紙では、左型も2つともはまる（作り分けられない）
from itertools import permutations, product
P = [(1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1)]
def apply(M,v): return tuple(M[i][0]*v[0]+M[i][1]*v[1]+M[i][2]*v[2] for i in range(3))
def det(M):
    return (M[0][0]*(M[1][1]*M[2][2]+-M[1][2]*M[2][1]) + -M[0][1]*(M[1][0]*M[2][2]+-M[1][2]*M[2][0])
            + M[0][2]*(M[1][0]*M[2][1]+-M[1][1]*M[2][0]))
rots = []
for perm in permutations(range(3)):
    for s in product((1,-1), repeat=3):
        M = [[0]*3 for _ in range(3)]
        for i in range(3): M[i][perm[i]] = s[i]
        if det(M) == 1 and sorted(apply(M,v) for v in P) == sorted(P): rots.append(M)
print("S1 回転の数", len(rots), "OK" if len(rots) == 12 else "NG")
R = {P[0]:"A", P[1]:"B", P[2]:"C", P[3]:"D"}
mirror = lambda v: (v[1], v[0], v[2])
L = {mirror(v): k for v, k in R.items()}
print("   左型", [L[v] for v in P], "（右型", [R[v] for v in P], "）")
def turned(mol, M): return {apply(M, v): k for v, k in mol.items()}
mL = [sum(1 for v in P if turned(L, M)[v] == R[v]) for M in rots]
mR = [sum(1 for v in P if turned(R, M)[v] == R[v]) for M in rots]
print("S2 左型の一致の最大", max(mL), "分布", sorted(mL), "OK" if max(mL) == 2 else "NG")
print("S3 右型で4個一致する回し方", mR.count(4), "OK" if mR.count(4) == 1 else "NG")
def socket(mol, need):
    best = 0
    for M in rots:
        t = turned(mol, M); best = max(best, sum(1 for v in need if t[v] == R[v]))
    return best
print("S4 受け口3つ: 右型", socket(R, P[1:]), "左型", socket(L, P[1:]),
      "OK" if (socket(R, P[1:]), socket(L, P[1:])) == (3, 2) else "NG")
print("S5 受け口2つ: 右型", socket(R, P[1:3]), "左型", socket(L, P[1:3]),
      "OK(作り分けられない)" if socket(L, P[1:3]) == 2 else "NG")
