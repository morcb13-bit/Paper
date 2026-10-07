# 扇10枚の担体の隙間（細ひし形・五芒星・舟・正十角形）を、床の五角形に対応づける
import sys, json, math, pickle
sys.path.insert(0, '/home/claude/Paper/2026-09-21-plant-sun/code')
from collections import Counter
import b13_chain_units as U
CODE = '/home/claude/Paper/2026-09-21-plant-sun/code/'
R = [[tuple(c) for c in r] for r in json.load(open(CODE + "R14.json"))]
z0 = (2, -2, 0, -3); t = (-11, 4, -4, 11)
def place(c, k):
    p = U.zadd(U.zrot(U.zsub(c, z0), k), z0)
    return U.zadd(p, U.zrot(t, k - 1)) if k % 2 else p
cells = {}
for k in range(10):
    for row in R:
        for c in row:
            for q, a in U.ring_cells(place(c, k)): cells[q] = a
cx, cy = U.xy(z0)
def build(rho=20.0):
    near = {q: a for q, a in cells.items() if math.hypot(U.xy(q)[0] - cx, U.xy(q)[1] - cy) < rho + 6}
    faces = U.gaps(near)
    return near, faces
if __name__ == "__main__":
    near, faces = build()
    pickle.dump((cells, z0), open('/home/claude/icolearn/cells10.pkl', 'wb'))
    print(len(cells), len(near), Counter(round(a, 4) for a, _ in faces))
