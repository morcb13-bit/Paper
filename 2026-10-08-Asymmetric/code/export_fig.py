import io, contextlib, json
with contextlib.redirect_stdout(io.StringIO()):
    exec(open("chiral_penrose.py").read().split('deg = {}')[0])
addr = {tuple(int(t) for t in k.split(",")): v for k, v in json.load(open("rings_integer.json"))["cells"].items()}
xy = []
for c in cells:
    z = T.to_xy(c)  # 観察者側の表示
    xy.append([round(z.real, 4), round(z.imag, 4), addr[c] % 2])
out = {"xy": xy, "nb": nb, "conj": conj, "pairs": pairs}
json.dump(out, open("fig_data.json", "w"), separators=(",", ":"))
# 照合用：いくつかの着地での最終
ref = {}
for e, s, k in ((5, 1, 0), (5, -1, 0), (100, 1, 0), (250, 1, 0), (5, 1, 11), (300, 1, 11)):
    x = seed(e, s); _, _, h = run(x, k=k); ref[f"{e},{s},{k}"] = h[-1]
x = seed(); _, _, h = run(x); ref["sym"] = h[-1]
print(json.dumps(ref, ensure_ascii=False))
