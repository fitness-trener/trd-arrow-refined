"""Invert the two-variable Euler product on the PUBLISHED A399030 b-file.

Predictions are pre-registered in wiki/concepts/plan-bfile-inversion.md and are not restated here.

Input: raw/sources/2026-09-11-b399030.txt, fetched from oeis.org/A399030/b399030.txt - rows n = 1..10,
every k = 1..n^2, 385 terms. No enumeration of any kind is done.

Peeling: every component has v >= 1, so the coefficient of x^e y^v in the product of the factors with
smaller v is exactly the contribution of multi-component objects; D(e,v) is what is left of T(v,e).
Two components of the same v cannot coexist in an object of v points, so all D(e,v) at one level v
are read off before any of them is multiplied in.

usage: python 2026-09-11-bfile-inversion.py
"""
import io
import json
import os
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
BFILE = os.path.join(HERE, '..', 'sources', '2026-09-11-b399030.txt')
CACHE = os.path.join(HERE, '2026-09-09-connected-classes.json')
NMAX = 10
KMAX = NMAX * NMAX
DMAX = NMAX - 2                    # C_d needs v <= d+2

vals = [int(line.split()[1]) for line in io.open(BFILE, encoding='utf-8')
        if line.strip() and not line.startswith('#')]
T, idx = {}, 0
for n in range(1, NMAX + 1):
    for k in range(1, n * n + 1):
        T[(n, k)] = vals[idx]
        idx += 1
assert idx == len(vals) == 385, (idx, len(vals))

P = [[0] * (KMAX + 1) for _ in range(NMAX + 1)]
P[0][0] = 1
D = {}
for v in range(1, NMAX + 1):
    level = {}
    for e in range(1, v * v + 1):
        d = T[(v, e)] - P[v][e]
        assert d >= 0, 'P2 FAIL: D(%d,%d) = %d is negative' % (e, v, d)
        if d:
            level[(e, v)] = d
    for (e, vv), cnt in sorted(level.items()):
        new = [row[:] for row in P]
        for n in range(NMAX + 1):
            for k in range(KMAX + 1):
                if not P[n][k]:
                    continue
                j = 1
                while n + j * vv <= NMAX and k + j * e <= KMAX:
                    new[n + j * vv][k + j * e] += P[n][k] * comb(cnt + j - 1, j)
                    j += 1
        P = new
    D.update(level)

print('inverted %d nonzero D(e,v) with v <= %d, %d connected classes in total'
      % (len(D), NMAX, sum(D.values())))

# P2 - the band v <= e+1
out_of_band = {ev: c for ev, c in D.items() if ev[1] > ev[0] + 1}
print('P2 band v <= e+1: %s' % ('PASS' if not out_of_band else 'FAIL %r' % out_of_band))

# P1 - against the brute-force enumeration, every e <= 6
brute = {tuple(map(int, k.split(','))): c for k, c in json.load(open(CACHE)).items()}
p1_bad = []
for e in range(1, 7):
    for v in range(1, NMAX + 1):
        if D.get((e, v), 0) != brute.get((e, v), 0):
            p1_bad.append((e, v, D.get((e, v), 0), brute.get((e, v), 0)))
print('P1 against brute force, e <= 6: %s'
      % ('PASS (%d classes agree)' % sum(brute.values()) if not p1_bad else 'FAIL %r' % p1_bad))

# C_d
C = {}
for d in range(0, DMAX + 1):
    tot = 0
    for e in range(max(1, -(-(d + 1) // 2)), d + 2):
        v = 2 * e - d
        assert 1 <= v <= NMAX, (d, e, v)
        tot += D.get((e, v), 0)
    C[d] = tot

# E_N and D_N
MAXJ = 2 * DMAX
E = [[0] * (MAXJ + 1) for _ in range(DMAX + 1)]
E[0][0] = 1
for (e, v), cnt in sorted(D.items()):
    d = 2 * e - v
    if d < 1 or d > DMAX:
        continue
    new = [row[:] for row in E]
    for dd in range(DMAX + 1):
        for j in range(MAXJ + 1):
            if not E[dd][j]:
                continue
            m = 1
            while dd + m * d <= DMAX and j + m * e <= MAXJ:
                new[dd + m * d][j + m * e] += E[dd][j] * comb(cnt + m - 1, m)
                m += 1
    E = new
DN = [sum(E[N]) for N in range(DMAX + 1)]

print('\nC_d, d = 0..%d : %s' % (DMAX, ', '.join(str(C[d]) for d in range(DMAX + 1))))
print('D_N, N = 0..%d : %s' % (DMAX, ', '.join(map(str, DN))))

print('\nP3 overlap with the recorded run:')
print('  C_0..C_5 %s' % ('PASS' if [C[d] for d in range(6)] == [1, 3, 5, 11, 22, 54] else 'FAIL'))
print('  D_0..D_5 %s' % ('PASS' if DN[:6] == [1, 3, 11, 36, 115, 357] else 'FAIL'))

print('\nP4 top two coefficients of E_N (Theorems 4.6, 4.7):')
p4_ok = True
for N in range(1, DMAX + 1):
    top = max(j for j in range(MAXJ + 1) if E[N][j])
    lead, second = E[N][2 * N], E[N][2 * N - 1]
    ok = top == 2 * N and lead == N + 1 and second == 4 * N - 3
    p4_ok &= ok
    print('  N=%d  deg=%-2d lead=%-3d (N+1=%-2d)  next=%-3d (4N-3=%-2d)  %s'
          % (N, top, lead, N + 1, second, 4 * N - 3, 'OK' if ok else 'FAIL'))
print('P4: %s' % ('PASS' if p4_ok else 'FAIL'))

print('\nE_N(x) for the new terms:')
for N in range(6, DMAX + 1):
    print('  E_%d = %s' % (N, ' + '.join('%d x^%d' % (E[N][j], j)
                                         for j in range(MAXJ + 1) if E[N][j])))

json.dump({'%d,%d' % ev: c for ev, c in sorted(D.items())},
          open(os.path.join(HERE, '2026-09-11-bfile-inversion-D.json'), 'w'))
