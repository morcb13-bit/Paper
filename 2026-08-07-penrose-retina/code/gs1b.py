#  検定GS1b  数回の視線の多数決で、閉じているかが決まるか
#      OK なら：閉じた形は多数決で閉じ、開いた形は多数決で開く（K を増やすと分離が上がる）
#      NG なら：K を増やしても混ざる
#      必ず落ちる設定：対象なし・一本の線が閉じと読まれないこと
exec(open('gs1.py').read().split('print("  図形')[0])
import random
V = {}
for nm, path in FIG2.items():
    V[nm] = [read(path, g)[1] for g in GAZ]
rnd = random.Random(13)
print("  多数決（視線を K 通りランダムに選び、内部≥1 が過半なら「閉じている」・試行2000回）")
print("  %-14s %7s %7s %7s %7s" % ("図形", "K=1", "K=3", "K=5", "K=7"))
for nm in FIG2:
    row = []
    for K in (1, 3, 5, 7):
        h = 0
        for _ in range(2000):
            s = rnd.sample(V[nm], K)
            if sum(1 for c in s if c >= 1) * 2 > K:
                h += 1
        row.append(h / 2000)
    print("  %-14s %7.3f %7.3f %7.3f %7.3f" % (nm, *row))
print("\n  内部の数の側（最頻値と、正しい値が出る割合）")
for nm, truth in (("閉じた輪(三角)", 1), ("二重の輪", 2)):
    c = {}
    for x in V[nm]:
        c[x] = c.get(x, 0) + 1
    mode = max(c, key=lambda k: c[k])
    print("  %-14s 最頻 %d（%d通り）／ %d が出る割合 %.3f"
          % (nm, mode, c[mode], truth, c.get(truth, 0) / len(GAZ)))
