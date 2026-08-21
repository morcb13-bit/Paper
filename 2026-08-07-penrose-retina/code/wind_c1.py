#  検定C1  巻き取りの残渣の散り方 ── 網膜を通したあと歪みと雑音を分けられるか
#
#      対象：三角形を ε だけ変形（または ε だけ雑音を載せ）てから網膜に着地させ、
#            観測された輪郭から理想の三角形を巻き取った残り（＝変位場）が
#            巻き数のどこに集まるかを見る。
#            変位場 = （歪み または 雑音）＋ 着地。着地は両方に等しく載る
#      測る量：r = Σ_{2≤|k|≤4}|c_k|² ÷ Σ_{2≤|k|≤N/2}|c_k|²
#              （k=0 と |k|=1 は型紙の位置・大きさ・向きが吸う。判定に使わない）
#      比べ方：歪みは軸の向き10通り（ζ10^m）、雑音は種40本。
#              歪みの帯が雑音の帯より上に完全に離れていれば分けられる
#
#      OK なら：ある ε から上で、歪み10本の最小値が雑音40本の最大値を超える。
#               その ε が「網膜を通して見分けられる歪みの下限」
#      NG なら：どの ε でも二つの帯が重なる。着地が歪みを潰しており、
#               巻き取りでは分けられない（v227 §2.2 の帰結どおり閾値は置けない）
#
#      必ず落ちる設定：ε=0 では歪みも雑音も「着地だけ」になり、
#                      二つの帯は完全に重なるはず（重ならなければ物差しが壊れている）
#
#      厳密性の内訳は wind_blank.py と同じ。
#      着地の判定は Q(φ) の符号だけ。巻き取りの評価だけ50桁の十進

exec(open('wind_core.py').read())

import statistics

SIDE_R = 11.0          # 辺長19.11。空試験で着地が雑音の帯の中に入っていた大きさ
CELLGAP = D("1.618033988749894848204586834365638117720309179805762862135")

V = []
for k in range(3):
    th = math.radians(90 + 120 * k)
    o, _ = find_offset(SIDE_R * math.cos(th), SIDE_R * math.sin(th))
    V.append(U.zadd(CTR, o))
SIDE = math.hypot(U.xy(V[1])[0] - U.xy(V[0])[0], U.xy(V[1])[1] - U.xy(V[0])[1])
M = max(3, int(round(SIDE / 1.618034)))
N = 3 * M
BASE = []
for k in range(3):
    a, b = V[k], V[(k + 1) % 3]
    for i in range(M):
        t = F(i, M)
        BASE.append(tuple(F(a[j]) + t * (b[j] - a[j]) for j in range(4)))
CEN = tuple(sum(p[j] for p in BASE) / N for j in range(4))
W = [(dcos(2 * PI * D(m) / D(N)), dsin(2 * PI * D(m) / D(N))) for m in range(N)]
FLAT = 6.0 / (2 * (N // 2 - 2) + 1)

def scale_to(ds, target_rms):
    """変位場の二乗平均を target_rms に合わせる有理数倍"""
    cur = ms(ds, N).sqrt()
    if cur == 0:
        return [tuple(F(0) for _ in range(4)) for _ in ds]
    k = F(int(target_rms / cur * 10 ** 12), 10 ** 12)
    return [tuple(k * c for c in d) for d in ds]

def observe(delta):
    """真の輪郭 = BASE + delta を網膜に着地させ、理想の輪郭からの変位場を返す"""
    out = []
    for p, d in zip(BASE, delta):
        q = tuple(F(x) for x in [c + e for c, e in zip(p, d)])
        out.append(U.zsub(tuple(F(x) for x in land(q)), p))
    return out

def r_of(ds):
    return float(ratio(wind(ds, N, W), N)) / FLAT

# 歪み：アフィン変形 z → z + ε·ζ^m·conj(z−重心)。m が剪断の軸を回す
def distort(m, eps):
    raw = [U.zrot(U.zconj(U.zsub(p, CEN)), m) for p in BASE]
    return scale_to(raw, eps)

def noise(seed, eps):
    random.seed(seed)
    raw = [tuple(F(random.randint(-1000, 1000), 1000) for _ in range(4)) for _ in range(N)]
    return scale_to(raw, eps)

print("辺長 %.2f・標本 %d 点（間隔は五角形の間隔にそろえてある）" % (SIDE, N))
print("取り分は一様値 %.4f で割って示す\n" % FLAT)
print("  ε/1.618   歪み10本（軸の向き）      雑音40本               判定")

def band(vals):
    vals = sorted(vals)
    return vals[0], statistics.median(vals), vals[-1]

for ratio_eps in (D(0), D("0.1"), D("0.25"), D("0.5"), D(1), D(2), D(4)):
    eps = ratio_eps * CELLGAP
    rD = [r_of(observe(distort(m, eps))) for m in range(10)]
    rB = [r_of(observe(noise(2000 + s, eps))) for s in range(40)]
    dl, dm, dh = band(rD)
    bl, bm, bh = band(rB)
    same = abs(dl - bl) < 1e-12 and abs(dh - bh) < 1e-12
    sep = "同じ" if same else ("分かれる" if dl > bh else "重なる")
    print("  %6.2f   %5.2f〜%5.2f (中央%5.2f)  %5.2f〜%5.2f (中央%5.2f)   %s"
          % (float(ratio_eps), dl, dh, dm, bl, bh, bm, sep))


#  検定C2  巻き取った残りをアフィン成分へ射影する
#      C1 は「低い巻き数の取り分」という帯だけを見た。歪みの形は分かっているのだから、
#      帯ではなく成分そのものへ当てるほうが強いはず（型紙合わせと同じ手）。
#      残りから相似の自由度（平行移動・大きさ・向き＝1 と z−重心）を先に抜き、
#      残った方向 conj(z−重心) への取り分 q を見る。
#      q = |<d,m>|² / (<m,m>·<d,d>)。雑音なら期待値 1/(N−3) = 1/33、歪みなら 1 に近い
#
#      OK なら：歪みと雑音が分かれる ε の下限が C1 より小さくなる
#      NG なら：小さくならない。帯で見ても成分で見ても同じで、C1 の限界が装置の限界
#      必ず落ちる設定：ε=0 では歪みも雑音も着地だけになり、二つは同じ値になる

def cdot(u, v):
    re = D(0); im = D(0)
    for (a, b), (c, d) in zip(u, v):
        re += a * c + b * d
        im += b * c - a * d
    return (re, im)

def cabs2(u):
    return cdot(u, u)[0]

def csub_proj(u, e):
    """u から e 方向を抜く（e は正規化前でよい）"""
    n = cabs2(e)
    if n == 0: return u
    re, im = cdot(u, e)
    return [(a - (re * c - im * d) / n, b - (re * d + im * c) / n) for (a, b), (c, d) in zip(u, e)]

P = [to_xy(p) for p in BASE]
CP = to_xy(CEN)
E0 = [(D(1), D(0))] * N
E1 = [(x - CP[0], y - CP[1]) for x, y in P]
E2 = [(x - CP[0], -(y - CP[1])) for x, y in P]
E1 = csub_proj(E1, E0)
MODE = csub_proj(csub_proj(E2, E0), E1)
UNIF = 1.0 / (N - 3)

def q_of(ds):
    d = [to_xy(x) for x in ds]
    d = csub_proj(csub_proj(d, E0), E1)
    n = cabs2(d)
    if n == 0: return 0.0
    re, im = cdot(d, MODE)
    return float((re * re + im * im) / (cabs2(MODE) * n))

print("\n検定C2 アフィン成分への取り分 q（雑音の期待値 1/%d = %.4f）\n" % (N - 3, UNIF))
print("  ε/1.618   歪み10本                雑音40本                判定")
for ratio_eps in (D(0), D("0.1"), D("0.25"), D("0.5"), D(1), D(2)):
    eps = ratio_eps * CELLGAP
    qD = [q_of(observe(distort(m, eps))) for m in range(10)]
    qB = [q_of(observe(noise(2000 + s, eps))) for s in range(40)]
    dl, dm, dh = band(qD)
    bl, bm, bh = band(qB)
    sep = "分かれる" if dl > bh else ("同じ" if abs(dl - bl) < 1e-12 and abs(dh - bh) < 1e-12 else "重なる")
    print("  %6.2f   %.4f〜%.4f (中央%.4f)  %.4f〜%.4f (中央%.4f)  %s"
          % (float(ratio_eps), dl, dh, dm, bl, bh, bm, sep))
