# 検定AL（交互桁）
#
#   監督の指定：桁は平衡3進のまま。+1 が現れるたびに右・左と向きを入れ替える。
#              戻るときはその逆。
#
#   盤面 V0  いまの盤面（交互なし）＝対照
#   盤面 V1  +1 だけ交互
#   盤面 V2  +1 と −1 の両方が交互
#
#   測るもの
#     (1) 歩けるか（受理率）と着地の個数        ── 盤面が壊れていないこと
#     (2) a=(1,0,1)、b=(1,−1,0)、ab、ba の着地
#     (3) 足し算になっているか：変位 disp(ab) = disp(a) + disp(b) か
#     (4) 可換か：disp(ab) = disp(ba) か
#     (5) 逆順で戻るか：語を逆順にして符号を反転して歩き、出発の辺に戻るか
#
#   OK(3)  disp(ab) = disp(a)+disp(b) が出発の大半で成り立つ
#   NG(3)  成り立たない（連結は平行移動の合成ではない）
#   必ず落ちる設定  空語では動かない（disp = 0）
#   負の対照  桁の付け替えをした盤面で同じ規則性が出ないこと
import pickle, math, random
from collections import Counter, defaultdict
import b13_chain_units as U

B = pickle.load(open("cm3_base.pkl", "rb"))
O, ST, CEN = B["O"], B["ST"], B["CEN"]
n = len(CEN)
XY = {}
for (p, q) in O:
    if p not in XY: XY[p] = U.xy(p)
    if q not in XY: XY[q] = U.xy(q)
CXY = [(U.xy(c)[0] / 10, U.xy(c)[1] / 10) for c in CEN]

def walk(w, e0, mode):
    """交互桁で歩く。mode 0=交互なし 1=+1のみ交互 2=両方交互。
       返り値は最後の有向辺、歩けなければ None"""
    prev, cur = e0
    used = Counter()
    for x in w:
        if x != 0 and (mode == 2 or (mode == 1 and x == 1)):
            eff = x if used[x] % 2 == 0 else -x
            used[x] += 1
        else:
            eff = x
        o = O[(prev, cur)]
        if eff not in o: return None
        prev, cur = cur, o[eff]
    return (prev, cur)

a = (1, 0, 1)
b = (1, -1, 0)
ab = a + b
ba = b + a

print("検定AL（交互桁）  五芒星 %d 個" % n)
print("a = %s   b = %s" % (str(a), str(b)))
print()

def run(w, mode):
    """全ての五芒星の全ての出発辺で歩く。変位のリストを返す"""
    out = []
    for i in range(n):
        for e0 in ST[i]:
            r = walk(w, e0, mode)
            if r is None: continue
            dx = XY[r[1]][0] - CXY[i][0] * 10
            dy = XY[r[1]][1] - CXY[i][1] * 10
            out.append(((i, e0), (round(dx, 6), round(dy, 6))))
    return dict(out)

for mode in (0, 1, 2):
    name = {0: "V0 交互なし", 1: "V1 +1だけ交互", 2: "V2 両方交互"}[mode]
    Da, Db, Dab, Dba = (run(w, mode) for w in (a, b, ab, ba))
    tot = sum(len(ST[i]) for i in range(n))
    print("%s   歩けた出発  a %d / b %d / ab %d / ba %d  （全 %d）"
          % (name, len(Da), len(Db), len(Dab), len(Dba), tot))
    # 必ず落ちる設定
    D0 = run((), mode)
    print("   必ず落ちる設定 空語 disp=0 :",
          "OK" if all(abs(v[0]) < 1e-6 and abs(v[1]) < 1e-6 for v in D0.values()) else "NG（出発は辺なので中心からずれる）")
    # (3) 足し算
    keys = [k for k in Dab if k in Da and k in Db]
    same = sum(1 for k in keys
               if abs(Dab[k][0] - (Da[k][0] + Db[k][0])) < 1e-6
               and abs(Dab[k][1] - (Da[k][1] + Db[k][1])) < 1e-6)
    print("   disp(ab) = disp(a)+disp(b) :  %d / %d" % (same, len(keys)))
    # (4) 可換
    k2 = [k for k in Dab if k in Dba]
    comm = sum(1 for k in k2 if abs(Dab[k][0] - Dba[k][0]) < 1e-6 and abs(Dab[k][1] - Dba[k][1]) < 1e-6)
    print("   disp(ab) = disp(ba) :         %d / %d" % (comm, len(k2)))
    # 着地の散らばり
    print("   ab の着地  相異なる変位 %d 種 / %d 本"
          % (len(set(Dab.values())), len(Dab)))
    print()

# --- (5) 逆順で戻るか ---
print("逆順で戻るか（語を逆順にし符号を反転して歩き、出発の辺に戻るか）")
for mode in (0, 1, 2):
    name = {0: "V0", 1: "V1", 2: "V2"}[mode]
    for w, wn in ((a, "a"), (ab, "ab")):
        rev = tuple(-x for x in reversed(w))
        back = 0; tried = 0
        for i in range(n):
            for e0 in ST[i]:
                r = walk(w, e0, mode)
                if r is None: continue
                tried += 1
                # 着いた辺を反転してから逆順の語を歩く
                r2 = walk(rev, (r[1], r[0]), mode)
                if r2 is not None and (r2[1], r2[0]) == e0: back += 1
        print("  %s  %-3s  戻った %d / %d" % (name, wn, back, tried))
print()

# --- 負の対照：桁の付け替え ---
print("負の対照（桁の付け替え：各頂点の出口の並びを無作為に振り直す）")
rng = random.Random(20260904)
O2 = {}
for e, o in O.items():
    ks = list(o.keys()); vs = list(o.values()); rng.shuffle(vs)
    O2[e] = dict(zip(ks, vs))
OSAVE = O
O = O2
for mode in (0, 1, 2):
    Da, Db, Dab = (run(w, mode) for w in (a, b, ab))
    keys = [k for k in Dab if k in Da and k in Db]
    same = sum(1 for k in keys
               if abs(Dab[k][0] - (Da[k][0] + Db[k][0])) < 1e-6
               and abs(Dab[k][1] - (Da[k][1] + Db[k][1])) < 1e-6)
    print("  mode %d   disp(ab)=disp(a)+disp(b)  %d / %d" % (mode, same, len(keys)))
O = OSAVE
