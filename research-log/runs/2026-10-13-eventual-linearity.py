"""Eventual linearity of the top coefficients of E_N  (written 2026-10-13, after the data existed).

STATUS OF THIS CHECK, stated before it runs: the pattern eps_N(2N-2) = 3, 15, 33, 51, 69, 87, 105
was SEEN in the E_N tables first; the proof below (cost grading) was derived afterwards and the
slope A_2 = 18 was computed by hand before this script ran.  So the comparison below verifies the
proof against data the author had already looked at.  It is NOT a pre-registered prediction and
must not be written up as one.

Claim.  Component type (e,v) has defect d = 2e-v and COST c = 2d - e = 3e - 2v.  Positive-defect
components have c >= d-1 >= 0, and c = 0 only for the two stars (e=2, v=3).  Total cost of a
multiset counted by eps_N(j) is sum(2d_i - e_i) = 2N - j.  Hence

    eps_N(2N - i) = sum over multisets M of NON-STAR positive-defect components with total cost i
                    of  (N - def(M) + 1)         [stars: multisets of N - def(M) items from 2 types]
                  = A_i * (N + 1) - B_i          for N >= 2i - 1,

with A_i = #M, B_i = sum of def(M).  Cost i needs types with c <= i, and c = 3e-2v with v <= e+1
gives e <= c+2, v <= c+3, so the b-file (v <= 10) fixes A_i, B_i exactly for i <= 7.

Input: raw/runs/2026-09-11-bfile-inversion-D.json  (keys "e,v" -> D(e,v), v <= 10).
"""
import json
import os
from collections import defaultdict
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
D = {tuple(map(int, k.split(','))): int(v)
     for k, v in json.load(open(os.path.join(HERE, '2026-09-11-bfile-inversion-D.json'))).items()}

IMAX = 7                                   # cost range that is exact with v <= 10
NMAX = 8                                   # defect range that is exact with v <= 10

# ---------- E_N(x) straight from D, no theory: product over (e, d) with d = 2e - v >= 1
# poly[(N, j)] = number of multisets with total defect N, total arrows j
types = []
for (e, v), cnt in sorted(D.items()):
    d = 2 * e - v
    if 1 <= d <= NMAX:
        types.append((e, d, cnt))

# coefficient of z^m in (1 - z)^(-cnt) is C(cnt + m - 1, m); cnt can be 10^8, so no loop over cnt.


def mult_type(poly, e, d, cnt):
    new = defaultdict(int)
    for (N, j), val in poly.items():
        m = 0
        while N + m * d <= NMAX:
            new[(N + m * d, j + m * e)] += val * comb(cnt + m - 1, m)
            m += 1
    return new


poly = defaultdict(int)
poly[(0, 0)] = 1
for e, d, cnt in types:
    poly = mult_type(poly, e, d, cnt)

E = {N: {j: poly[(N, j)] for j in range(0, 2 * N + 1) if poly.get((N, j))} for N in range(0, NMAX + 1)}

PAPER = {
    1: {1: 1, 2: 2},
    2: {2: 3, 3: 5, 4: 3},
    3: {3: 8, 4: 15, 5: 9, 6: 4},
    4: {3: 1, 4: 22, 5: 41, 6: 33, 7: 13, 8: 5},
    5: {4: 8, 5: 59, 6: 118, 7: 98, 8: 51, 9: 17, 10: 6},
    6: {4: 1, 5: 33, 6: 177, 7: 325, 8: 301, 9: 165, 10: 69, 11: 21, 12: 7},
    7: {5: 6, 6: 126, 7: 520, 8: 931, 9: 870, 10: 530, 11: 232, 12: 87, 13: 25, 14: 8},
    8: {6: 42, 7: 444, 8: 1577, 9: 2648, 10: 2553, 11: 1604, 12: 774, 13: 299, 14: 105, 15: 29, 16: 9},
}
print('E_N from D(e,v) against the tables typed into the paper')
for N in range(1, NMAX + 1):
    ok = E[N] == PAPER[N]
    print('  N=%d  sum=%-6d  %s' % (N, sum(E[N].values()), 'OK' if ok else 'MISMATCH'))
    assert ok

# ---------- cost grading: A_i, B_i
nonstar = []
for (e, v), cnt in D.items():
    d = 2 * e - v
    c = 3 * e - 2 * v
    if d >= 1 and c >= 1 and c <= IMAX:
        nonstar.append((c, d, cnt))
    if d >= 1 and c < 0:
        raise SystemExit('negative cost with positive defect: (%d,%d)' % (e, v))
assert all(not (2 * e - v >= 1 and 3 * e - 2 * v == 0) or (e, v) == (2, 3) for e, v in D)

# gf[(i, dsum)] = number of multisets of non-star types with total cost i and total defect dsum
gf = defaultdict(int)
gf[(0, 0)] = 1
DMAX = 2 * IMAX + 2
for c, d, cnt in nonstar:
    new = defaultdict(int)
    for (i, ds), val in gf.items():
        m = 0
        while i + m * c <= IMAX and ds + m * d <= DMAX:
            new[(i + m * c, ds + m * d)] += val * comb(cnt + m - 1, m)
            m += 1
    gf = new

A, B = [], []
for i in range(IMAX + 1):
    A.append(sum(v for (ii, ds), v in gf.items() if ii == i))
    B.append(sum(v * ds for (ii, ds), v in gf.items() if ii == i))
print('\ncost grading, i = 0..%d' % IMAX)
print('  A_i =', A)
print('  B_i =', B)
print('  max defect of a counted multiset at cost i:',
      [max(ds for (ii, ds) in gf if ii == i) for i in range(IMAX + 1)], ' (<= 2i, hence N >= 2i-1)')

# ---------- compare:  eps_N(2N - i) = A_i (N+1) - B_i  for N >= 2i - 1
print('\neps_N(2N-i) against A_i (N+1) - B_i   (checked where N >= 2i-1 and N <= %d)' % NMAX)
checked = failed = 0
for N in range(1, NMAX + 1):
    for i in range(0, IMAX + 1):
        if N < 2 * i - 1 or 2 * N - i < 0:
            continue
        got = E[N].get(2 * N - i, 0)
        want = A[i] * (N + 1) - B[i]
        checked += 1
        if got != want:
            failed += 1
        print('  N=%d i=%d  eps=%-5d  formula=%-5d  %s' % (N, i, got, want, 'ok' if got == want else 'FAIL'))
print('checked %d, failed %d' % (checked, failed))
assert failed == 0

# where the formula is expected to fail (N < 2i-1): record, do not hide
print('\nbelow the threshold (expected to differ somewhere):')
for N in range(1, NMAX + 1):
    for i in range(0, IMAX + 1):
        if N >= 2 * i - 1 or 2 * N - i < 0:
            continue
        got = E[N].get(2 * N - i, 0)
        want = A[i] * (N + 1) - B[i]
        print('  N=%d i=%d  eps=%-5d  formula=%-5d  %s' % (N, i, got, want, 'same' if got == want else 'differs'))
