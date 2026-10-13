"""Package the connected triangle for OEIS, and verify it five ways before any field is written.

Verification, not an experiment: the values are the inverse bivariate Euler transform of A399030,
already recovered by 2026-09-11-bfile-inversion.py. Nothing here is predicted.

usage: python 2026-09-11-connected-triangle.py
"""
import importlib.util
import io
import json
import os
from itertools import permutations
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, '..', 'sources', '2026-09-11-b399030.txt')
A296105 = [1, 2, 5, 25, 157, 1325, 14358, 199763, 3549001, 80673244, 2352747542]
V = 10

D = {tuple(map(int, k.split(','))): c
     for k, c in json.load(open(os.path.join(HERE, '2026-09-11-bfile-inversion-D.json'))).items()}
rows = {n: [D.get((k, n), 0) for k in range(1, n * n + 1)] for n in range(1, V + 1)}
flat = [x for n in range(1, V + 1) for x in rows[n]]
pub = [int(l.split()[1]) for l in io.open(SRC, encoding='utf-8') if l.strip() and l[0] != '#']
A = {}
i = 0
for n in range(1, V + 1):
    for k in range(1, n * n + 1):
        A[n, k] = pub[i]; i += 1

ok = True
def check(name, cond, extra=''):
    global ok
    ok &= bool(cond)
    print('  %-58s %s %s' % (name, 'PASS' if cond else 'FAIL', extra))

print('checks')
check('385 entries, row n has n^2 of them', len(flat) == 385)
check('all entries non-negative', min(flat) >= 0)
check('band: T(n,k) = 0 whenever k < n-1',
      all(rows[n][k - 1] == 0 for n in rows for k in range(1, n * n + 1) if k < n - 1))
check('row sums = A296105(n) for 2 <= n <= 10',
      all(sum(rows[n]) == A296105[n] for n in range(2, V + 1)),
      '(row 1 = %d; A296105(1) = 2 counts the empty relation)' % sum(rows[1]))
check('k > (n-1)^2+1 forces connected: T(n,k) = A399030(n,k)',
      all(rows[n][k - 1] == A[n, k] for n in rows for k in range((n - 1) ** 2 + 2, n * n + 1)))

# round trip: bivariate Euler transform of the triangle must give back A399030 exactly
P = [[0] * (V * V + 1) for _ in range(V + 1)]
P[0][0] = 1
for (e, v), m in sorted(D.items()):
    Q = [r[:] for r in P]
    for n in range(V + 1):
        for k in range(V * V + 1):
            j = 1
            while P[n][k] and n + j * v <= V and k + j * e <= V * V:
                Q[n + j * v][k + j * e] += P[n][k] * comb(m + j - 1, j); j += 1
    P = Q
check('Euler transform of the triangle = A399030, all 385 terms',
      all(P[n][k] == A[n, k] for (n, k) in A))

# independent brute force of rows 1..4, every k
def brute_row(n):
    U = [(a, b) for a in range(n) for b in range(n)]
    perms = list(permutations(range(n)))
    seen = {}
    for mask in range(1, 1 << len(U)):
        R = [U[i] for i in range(len(U)) if mask >> i & 1]
        S = set(R)
        if {x for p in R for x in p} != set(range(n)): continue
        if any((a, d) not in S for a, b in R for c, d in R if b == c): continue
        adj = {x: set() for x in range(n)}
        for a, b in R: adj[a].add(b); adj[b].add(a)
        st, sn = [0], {0}
        while st:
            for y in adj[st.pop()] - sn: sn.add(y); st.append(y)
        if len(sn) < n: continue
        c = min(tuple(sorted((p[a], p[b]) for a, b in R)) for p in perms)
        seen.setdefault(len(R), set()).add(c)
    return [len(seen.get(k, ())) for k in range(1, n * n + 1)]
bf = {n: brute_row(n) for n in range(1, 5)}
check('brute force rows 1..4, every k, = triangle', all(bf[n] == rows[n] for n in bf),
      '(%d cells)' % sum(n * n for n in bf))

# the OEIS PROG must print exactly these rows
spec = importlib.util.spec_from_file_location('prog', os.path.join(HERE, '2026-09-11-oeis-prog-triangle.py'))
prog = importlib.util.module_from_spec(spec); spec.loader.exec_module(prog)
check('OEIS PROG output = triangle', prog.connected_triangle(SRC) == [rows[n] for n in range(1, V + 1)])

print('\nALL CHECKS: %s' % ('PASS' if ok else 'FAIL'))
assert ok

with open(os.path.join(HERE, '2026-09-11-connected-triangle-bfile.txt'), 'w', newline='\n') as f:
    for idx, x in enumerate(flat, 1):
        f.write('%d %d\n' % (idx, x))

data, s = [], ''
for x in flat:
    nxt = (s + ', ' if s else '') + str(x)
    if len(nxt) > 260: break
    s = nxt; data.append(x)
print('\nDATA (%d terms, %d chars):\n%s' % (len(data), len(s), s))
print('\nEXAMPLE rows:')
for n in range(1, 6):
    print('  n=%d: %s;' % (n, ', '.join(map(str, rows[n]))))
print('\nrow sums n=1..10: %s' % [sum(rows[n]) for n in range(1, V + 1)])
print('column k=n-1 (spanning trees, first entries): %s' % [rows[n][n - 2] for n in range(2, V + 1)])
print('b-file written: %d lines' % len(flat))
