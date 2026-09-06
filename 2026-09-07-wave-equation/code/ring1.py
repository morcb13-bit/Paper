from b13_two_tilings import zt, zadd, zrot, PHI2, ONE, zsub, norm2, phi_cmp

def build_ring1():
    cells = {}                      # 位置Zeta -> 番地
    for k in range(10):
        cells[zadd((0,0,0,0), zrot(PHI2, k))] = k
    verts, edges = set(), set()
    for p, a in cells.items():
        o = a % 2
        cs = [zadd(p, zt(2*j+o)) for j in range(5)]
        for c in cs: verts.add(c)
        for j in range(5):
            e = tuple(sorted([cs[j], cs[(j+1) % 5]]))
            edges.add(e)
    return sorted(verts), sorted(edges)

if __name__ == "__main__":
    V, E = build_ring1()
    deg = {v: 0 for v in V}
    for a, b in E:
        deg[a] += 1; deg[b] += 1
    from collections import Counter
    print("五角形10 / 頂点", len(V), "/ 無向辺", len(E), "/ 有向辺", 2*len(E))
    print("次数分布", dict(sorted(Counter(deg.values()).items())))
