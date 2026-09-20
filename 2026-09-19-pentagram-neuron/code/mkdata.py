import math, itertools, json
exec(open('wind_core.py').read())
exec(open('pentagram_neuron.py').read().split('def run_tests')[0])
src = open('amoeba.py').read()
exec(src[src.index('def build'):src.index('def run_amoeba')])
SC, AROUND, _ = build()
F, faces, SCx = carrier()
XY = {q: tuple(float(t) for t in U.xy(q)) for q in F}
cells = [[round(x, 4), round(y, 4)] for x, y in XY.values()]
stars = [[round(x, 4), round(y, 4)] for x, y in SC]
around = [[[round(p[0], 4), round(p[1], 4)] for q, p in A] for A in AROUND]
links = [[i, j] for i, j in itertools.combinations(range(30), 2)
         if math.dist(SC[i], SC[j]) <= 10.5784 + 1e-3]
axis = [i for i in range(30) if abs(SC[i][0]) < 1e-9]
json.dump({"cells": cells, "stars": stars, "around": around,
           "links": links, "axis": axis},
          open('carrier.json', 'w'), separators=(',', ':'))
print(len(cells), len(stars), len(links), axis)
