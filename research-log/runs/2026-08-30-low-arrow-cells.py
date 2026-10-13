"""T(11,6), T(11,7), T(11,8) directly - the three cells the n=11 run exists to produce.

Pre-registered predictions live in wiki/concepts/plan-unlabeled-n11.md §3 and are NOT restated
here, so this script cannot quietly become the place they get edited.

WHY THIS WORKS WITHOUT THE 8-DAY RUN. A TRD with no isolated point is a multiset of connected TRDs,
each carrying at least one arrow. A connected relation on v points needs at least v-1 arrows in its
underlying graph, so v <= e+1: every component "wastes" e+1-v points of capacity.

For n = 11 with k arrows in r components: 11 = sum v_i <= sum (e_i + 1) = k + r, and r <= k since
each component has e_i >= 1. So

    k = 6 -> r in {5,6},  max single component e = 6 - (r-1) <= 2
    k = 7 -> r in {4..7}, max e <= 4
    k = 8 -> r in {3..8}, max e <= 6

**Every component needed has e <= 6, hence v <= 7.** That is a brute-force-sized problem
(C(49,6) = 14M subsets) instead of a 46.7-million-poset walk.

METHOD
  1. D(e,v) = connected transitive relations with e arrows on v points, no isolated point, up to
     isomorphism. Brute force over e-subsets of the v^2 ordered pairs, then canonical form.
  2. Multiset assembly. With D(e,v) distinguishable classes of type (e,v), the generating function
     is prod_{e,v} (1 - x^e y^v)^{-D(e,v)}; T(n,k) = [x^k y^n]. Done as a DP, choosing j copies from
     D classes with repetition in C(D+j-1, j) ways.

ORACLE, and it is broad. Every T(n,k) with k <= 6 is *fully determined* by D(e<=6, v<=7), because a
decomposition of k arrows cannot contain a component with more than k. So the whole k <= 6 corner of
the shipped table 2026-08-04-trd-bivariate-n10.csv must be reproduced exactly, for every n <= 10,
before any n = 11 value is believed.

usage: python 2026-08-30-low-arrow-cells.py
"""
import csv
import os
import sys
import time
from itertools import combinations, permutations
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
SHIPPED = os.path.join(HERE, '2026-08-04-trd-bivariate-n10.csv')

MAX_E = 6          # justified above: no component with more than 6 arrows can appear
MAX_V = MAX_E + 1  # connected => v <= e+1


def is_transitive(pairs):
    s = set(pairs)
    for (a, b) in pairs:
        for (c, d) in pairs:
            if b == c and (a, d) not in s:
                return False
    return True


def touched(pairs):
    t = set()
    for (a, b) in pairs:
        t.add(a)
        t.add(b)
    return t


def weakly_connected(pairs, pts):
    adj = {p: set() for p in pts}
    for (a, b) in pairs:
        adj[a].add(b)
        adj[b].add(a)
    start = next(iter(pts))
    seen, stack = {start}, [start]
    while stack:
        x = stack.pop()
        for y in adj[x]:
            if y not in seen:
                seen.add(y)
                stack.append(y)
    return seen == pts


def canon(v, pairs):
    best = None
    for p in permutations(range(v)):
        e = tuple(sorted((p[a], p[b]) for (a, b) in pairs))
        if best is None or e < best:
            best = e
    return best


def connected_classes():
    """D[(e,v)] = number of connected TRD isomorphism classes with e arrows on v points."""
    D = {}
    for v in range(1, MAX_V + 1):
        universe = [(a, b) for a in range(v) for b in range(v)]
        for e in range(1, MAX_E + 1):
            if v > e + 1:
                continue
            t0 = time.time()
            seen = set()
            for pick in combinations(universe, e):
                pts = touched(pick)
                if len(pts) != v:
                    continue
                if not is_transitive(pick):
                    continue
                if not weakly_connected(pick, pts):
                    continue
                seen.add(canon(v, pick))
            if seen:
                D[(e, v)] = len(seen)
            print('  e=%d v=%d -> %-6d classes   (%.1fs)'
                  % (e, v, len(seen), time.time() - t0), flush=True)
    return D


def assemble(D, max_n, max_k):
    """T[n][k] from the two-variable Euler product over component classes."""
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


def shipped_table():
    rows = {}
    with open(SHIPPED) as f:
        for r in csv.DictReader(f):
            rows.setdefault(int(r['objects']), {})[int(r['arrows'])] = \
                int(r['arrow_type_semigroupoids'])
    return rows


def main():
    print('connected component classes, e <= %d, v <= %d' % (MAX_E, MAX_V), flush=True)
    t0 = time.time()
    D = connected_classes()
    print('  total %.1fs; %d types, %d classes\n'
          % (time.time() - t0, len(D), sum(D.values())), flush=True)

    T = assemble(D, max_n=11, max_k=8)

    print('ORACLE - every T(n,k) with k <= 6, n <= 10, against the shipped table')
    print('(k <= 6 is fully determined by these components: no part can exceed k arrows)')
    ship = shipped_table()
    bad = checked = 0
    for n in range(1, 11):
        for k in range(1, 7):
            want = ship.get(n, {}).get(k, 0)
            got = T[n][k]
            checked += 1
            if want != got:
                bad += 1
                print('  MISMATCH n=%d k=%d: got %d, shipped %d' % (n, k, got, want))
    print('  %d cells checked, %d mismatches -> %s\n'
          % (checked, bad, 'PASS' if bad == 0 else 'FAIL'))
    assert bad == 0, 'oracle failed; the n=11 values are not to be believed'

    print('THE THREE CELLS')
    for k in (6, 7, 8):
        print('  T(11,%d) = %d' % (k, T[11][k]))

    print('\nfull low-arrow corner, T(n,k) for k = 1..8')
    print('   n |' + ''.join('%9d' % k for k in range(1, 9)))
    for n in range(1, 12):
        mark = '  <- new' if n == 11 else ''
        print('  %2d |' % n + ''.join('%9d' % T[n][k] for k in range(1, 9)) + mark)


if __name__ == '__main__':
    main()
