"""E_N(x): defect-N multisets of positive-defect connected objects, refined by total arrow count.

T(2k-N, k) = sum of the coefficients of x^j in E_N(x) over j <= k, because an object of defect N with
k arrows is a multiset of positive-defect components (total arrows j) padded with k-j single arrows,
and the padding exists iff j <= k. D_N = E_N(1), and deg E_N = 2N, which is the theorem.

Every partial sum is checked against the published A399030 values wherever the cell is determined by
components with at most MAX_E arrows.

usage: python 2026-09-09-defect-polynomials.py
"""
import csv
import json
import os
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
MAX_E = 6
D = {tuple(map(int, k.split(','))): v
     for k, v in json.load(open(os.path.join(HERE, '2026-09-09-connected-classes.json'))).items()}

MAXD = MAX_E - 1                       # C_d complete only for d <= MAX_E - 1
MAXJ = 2 * MAXD

# P[d][j] = number of multisets of positive-defect connected objects, total defect d, total arrows j
P = [[0] * (MAXJ + 1) for _ in range(MAXD + 1)]
P[0][0] = 1
for (e, v), cnt in sorted(D.items()):
    d = 2 * e - v
    if d < 1 or d > MAXD:
        continue
    new = [row[:] for row in P]
    for dd in range(MAXD + 1):
        for j in range(MAXJ + 1):
            if not P[dd][j]:
                continue
            m = 1
            while dd + m * d <= MAXD and j + m * e <= MAXJ:
                new[dd + m * d][j + m * e] += P[dd][j] * comb(cnt + m - 1, m)
                m += 1
    P = new

ship = {}
with open(os.path.join(HERE, '2026-08-04-trd-bivariate-n10.csv')) as f:
    for r in csv.DictReader(f):
        ship.setdefault(int(r['objects']), {})[int(r['arrows'])] = int(r['arrow_type_semigroupoids'])


def valid(n, k):
    return k >= 1 and n >= 1 and k - (max(1, n - k) - 1) <= MAX_E


print('E_N(x), coefficient of x^j = defect-N multisets with j arrows in total')
for N in range(MAXD + 1):
    terms = [(j, P[N][j]) for j in range(MAXJ + 1) if P[N][j]]
    print('  E_%d = %s' % (N, ' + '.join('%d x^%d' % (c, j) for j, c in terms)))
    print('        D_%d = E_%d(1) = %d, degree %d %s'
          % (N, N, sum(c for _, c in terms), terms[-1][0],
             '(= 2N)' if terms[-1][0] == 2 * N else '(NOT 2N)'))

print('\npartial sums T(2k-N, k) = sum_{j<=k} [x^j] E_N, against published A399030')
bad = ok = 0
for N in range(MAXD + 1):
    line = []
    for k in range(1, 12):
        n = 2 * k - N
        if n < 1:
            continue
        pred = sum(P[N][j] for j in range(min(k, MAXJ) + 1))
        pub = ship.get(n, {}).get(k)
        if pub is None or not valid(n, k):
            line.append('T(%d,%d)=%d' % (n, k, pred))
            continue
        mark = 'OK' if pub == pred else 'MISMATCH(pub=%d)' % pub
        if pub == pred:
            ok += 1
        else:
            bad += 1
        line.append('T(%d,%d)=%d[%s]' % (n, k, pred, mark))
    print('  N=%d  %s' % (N, '  '.join(line)))
print('\n%d published cells matched, %d mismatched -> %s'
      % (ok, bad, 'PASS' if bad == 0 else 'FAIL'))
assert bad == 0

print('\nWhat each published cell pins')
for N in range(1, MAXD + 1):
    tail = sum(P[N][j] for j in range(MAXJ + 1)) - 0
    for k in range(1, 12):
        n = 2 * k - N
        if n < 1 or not valid(n, k) or ship.get(n, {}).get(k) is None:
            continue
        miss = sum(P[N][j] for j in range(k + 1, MAXJ + 1))
        if miss:
            print('  D_%d = T(%d,%d) + %d = %d + %d = %d   [published cell + tail]'
                  % (N, n, k, miss, ship[n][k], miss, tail))
            break
