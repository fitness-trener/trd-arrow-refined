"""Low-arrow corner and defect diagonals of A399030, WITH the validity mask v1 lacked.

v1 (2026-09-09-corner-table.py) printed columns k = 7, 8 without marking that they are only correct
where the point count forces every component to have at most MAX_E = 6 arrows. It does not: v1's
T(4,7) = 4 against a published 35. The output of v1 is kept as the record of that mistake.

VALIDITY. An object with n points and k arrows has r components with n <= sum(e_i + 1) = k + r, so
r >= n - k, and the largest single component has at most e_max = k - (r_min - 1) arrows. The cell is
determined by the enumerated components iff e_max <= MAX_E.

usage: python 2026-09-09-corner-table-v2.py
"""
import csv
import json
import os
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
MAX_E = 6
MAX_N, MAX_K = 11, 8

D = {tuple(map(int, k.split(','))): v
     for k, v in json.load(open(os.path.join(HERE, '2026-09-09-connected-classes.json'))).items()}
assert max(e for e, v in D) == MAX_E


def valid(n, k):
    if k < 1 or n < 1:
        return False
    r_min = max(1, n - k)
    return k - (r_min - 1) <= MAX_E


def assemble(max_n, max_k):
    T = [[0] * (max_k + 1) for _ in range(max_n + 1)]
    T[0][0] = 1
    for (e, v), cnt in sorted(D.items()):
        new = [row[:] for row in T]
        for n in range(max_n + 1):
            for k in range(max_k + 1):
                if not T[n][k]:
                    continue
                j = 1
                while n + j * v <= max_n and k + j * e <= max_k:
                    new[n + j * v][k + j * e] += T[n][k] * comb(cnt + j - 1, j)
                    j += 1
        T = new
    return T


T = assemble(MAX_N, MAX_K)

ship = {}
with open(os.path.join(HERE, '2026-08-04-trd-bivariate-n10.csv')) as f:
    for r in csv.DictReader(f):
        ship.setdefault(int(r['objects']), {})[int(r['arrows'])] = int(r['arrow_type_semigroupoids'])

print('ORACLE: every cell the mask calls valid, n <= 10, against the shipped table')
bad = checked = 0
for n in range(1, 11):
    for k in range(1, MAX_K + 1):
        if not valid(n, k):
            continue
        checked += 1
        if ship.get(n, {}).get(k, 0) != T[n][k]:
            bad += 1
            print('  MISMATCH n=%d k=%d: got %d, shipped %d'
                  % (n, k, T[n][k], ship.get(n, {}).get(k, 0)))
print('  %d checked, %d mismatches -> %s' % (checked, bad, 'PASS' if bad == 0 else 'FAIL'))
assert bad == 0

print('\nCOUNTER-CHECK: cells the mask calls invalid must generally DISagree with the shipped table')
dis = same = 0
for n in range(1, 11):
    for k in range(1, MAX_K + 1):
        if valid(n, k):
            continue
        w = ship.get(n, {}).get(k, 0)
        if w != T[n][k]:
            dis += 1
        else:
            same += 1
print('  %d invalid cells: %d differ from shipped (as they must), %d coincide' % (dis + same, dis, same))

print('\nlow-arrow corner, "." = not determined by components with e <= %d' % MAX_E)
print('   n |' + ''.join('%9d' % k for k in range(1, MAX_K + 1)))
for n in range(1, MAX_N + 1):
    cells = ''.join(('%9d' % T[n][k]) if valid(n, k) else '%9s' % '.'
                    for k in range(1, MAX_K + 1))
    print('  %2d |' % n + cells + ('   <- beyond the published b-file' if n == 11 else ''))

print('\ndefect diagonals T(2k-N, k), valid cells only; D_N is the stable value')
DN = [1, 3, 11, 36, 115, 357]
for N in range(0, 6):
    row = [(k, T[2 * k - N][k]) for k in range(1, MAX_K + 1)
           if 1 <= 2 * k - N <= MAX_N and valid(2 * k - N, k)]
    print('  N=%d  D_N=%-5d %s' % (N, DN[N],
          ',  '.join('T(%d,%d)=%d%s' % (2 * k - N, k, v, '*' if k >= 2 * N else '')
                     for k, v in row)))
    mono = all(row[i][1] <= row[i + 1][1] for i in range(len(row) - 1))
    stable = {v for k, v in row if k >= 2 * N}
    print('       monotone: %s | k>=2N values: %s | equals D_N: %s'
          % (mono, stable or 'none in range',
             'yes' if stable == {DN[N]} else ('n/a' if not stable else 'NO')))
    below = [(k, v) for k, v in row if k < 2 * N]
    if below:
        k, v = below[-1]
        print('       largest k below the threshold: k=%d (2N-%d), T=%d %s D_N=%d'
              % (k, 2 * N - k, v, '<' if v < DN[N] else '>=', DN[N]))
