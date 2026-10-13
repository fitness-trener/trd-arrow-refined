from math import comb
def connected_triangle(bfile='b399030.txt'):  # needs the b-file of A399030 (rows 1..10)
    t = [int(l.split()[1]) for l in open(bfile) if l.strip() and l[0] != '#']
    V, K, i, T = 10, 100, 0, {}
    for n in range(1, V + 1):
        for k in range(1, n * n + 1): T[n, k] = t[i]; i += 1
    P = [[0] * (K + 1) for _ in range(V + 1)]; P[0][0] = 1; D = {}
    for v in range(1, V + 1):
        lev = {(e, v): T[v, e] - P[v][e] for e in range(1, v * v + 1) if T[v, e] != P[v][e]}
        for (e, _), m in lev.items():
            Q = [r[:] for r in P]
            for n in range(V + 1):
                for k in range(K + 1):
                    j = 1
                    while P[n][k] and n + j * v <= V and k + j * e <= K:
                        Q[n + j * v][k + j * e] += P[n][k] * comb(m + j - 1, j); j += 1
            P = Q
        D.update(lev)
    return [[D.get((k, n), 0) for k in range(1, n * n + 1)] for n in range(1, V + 1)]
