# 検定S2：友達リンク（SNAP ego-Facebook）上で「一刻の遅れ」の紹介規則を走らせる
# 規則F：step s で知り合った人は step s+2 から、毎刻 1 人ずつ、まだ届いていない友達を紹介する
# 対照D：遅れなし（step s+1 から毎刻 1 人）
# 対照B：遅れなし・全員に一度に（1刻に1歩の匂い＝幅優先）
# 負の対照：次数を保ったまま辺を張り替えたグラフ（二重辺交換）で同じことをする
# 判定はすべて整数の数え上げ
import gzip, random, sys, json
import networkx as nx

def load():
    G = nx.Graph()
    with gzip.open('fb.txt.gz', 'rt') as f:
        for line in f:
            a, b = line.split()
            G.add_edge(int(a), int(b))
    return G

def spread(adj, start, delay, order_seed, cap=5000):
    """delay=2: 規則F, delay=1: 対照D。返り値：各刻の届いた人数の列"""
    rng = random.Random(order_seed)
    nbr = {v: rng.sample(adj[v], len(adj[v])) for v in adj}   # 紹介の順番（固定の乱順）
    ptr = {v: 0 for v in adj}
    reached = {start: 1}          # 人 -> 届いた刻
    order = [start]
    counts = [1]                  # counts[0] = 刻1
    n = len(adj)
    t = 1
    while len(reached) < n and t < cap:
        t += 1
        new = []
        for v in order:
            if t - reached[v] < delay:
                continue
            L = nbr[v]; p = ptr[v]
            while p < len(L) and L[p] in reached:
                p += 1
            ptr[v] = p
            if p < len(L):
                u = L[p]; ptr[v] = p + 1
                reached[u] = t
                new.append(u)
        order.extend(new)
        counts.append(len(reached))
    return counts

def fib_run(counts):
    """counts が 1,1,2,3,5,... と一致している刻の長さ"""
    a, b = 1, 1
    k = 0
    for c in counts:
        if c != a:
            break
        k += 1
        a, b = b, a + b
    return k

def first_at_least(counts, num, den, n):
    # den*counts >= num*n となる最初の刻（整数比較）
    for i, c in enumerate(counts):
        if den * c >= num * n:
            return i + 1
    return None

def bfs_ecc(adj, s):
    dist = {s: 0}; fr = [s]; d = 0
    while fr:
        d += 1; nx_ = []
        for v in fr:
            for u in adj[v]:
                if u not in dist:
                    dist[u] = d; nx_.append(u)
        fr = nx_
    return d - 1

def summarize(xs):
    xs = sorted(xs)
    m = len(xs)
    return dict(min=xs[0], q1=xs[m // 4], med=xs[m // 2], q3=xs[(3 * m) // 4], max=xs[-1])

def run(G, label, starts):
    adj = {v: list(G[v]) for v in G}
    n = len(adj)
    out = {}
    for name, delay in (('F', 2), ('D', 1)):
        fr, half, full = [], [], []
        for s in starts:
            c = spread(adj, s, delay, order_seed=s)
            fr.append(fib_run(c) if delay == 2 else 0)
            half.append(first_at_least(c, 1, 2, n))
            full.append(len(c))
        out[name] = dict(fib_run=summarize(fr) if delay == 2 else None,
                         half=summarize(half), full=summarize(full))
    out['B_ecc'] = summarize([bfs_ecc(adj, s) for s in starts])
    print(label, json.dumps(out, ensure_ascii=False))
    return out

if __name__ == '__main__':
    G = load()
    print('人数', G.number_of_nodes(), '友達関係', G.number_of_edges(),
          '連結成分', nx.number_connected_components(G))
    rng = random.Random(13)
    starts = rng.sample(sorted(G.nodes()), 200)
    run(G, '実データ', starts)
    R = G.copy()
    nx.double_edge_swap(R, nswap=10 * R.number_of_edges(), max_tries=10**8, seed=13)
    comp = max(nx.connected_components(R), key=len)
    print('張り替え後の最大成分', len(comp))
    R = R.subgraph(comp).copy()
    starts_r = [s for s in starts if s in R]
    run(R, '張り替え（次数を保つ）', starts_r)
