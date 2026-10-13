"""Labeled transitive relations by (points, pairs) - A399095 - on the STREAMING catalogue.

Supersedes 2026-08-26-labeled-triangle-fast.py, whose twin-class optimization measured as a wash
(m=9 catalogue 158.8s against <=161s for the unoptimized trd_bivariate route) and is recorded as a
negative result. The reason it failed: trd_bivariate.canon_and_aut is only slow on highly symmetric
posets, which are rare, so a per-poset twin computation taxes the common path to speed up the tail.

trd_stream.py already had the right fix and it was in this directory the whole time:

  * McKay canonical augmentation - each isomorphism class reached exactly once, no global `seen`
    set, memory O(max_m) posets instead of O(A000112(max_m)). That is what makes m=11 attemptable:
    probe P measured 1.36 KB/poset, i.e. ~64 GB at m=11 against 8.35 GB physical.
  * _full_symmetric shortcut - when Aut(P) is the full product of symmetric groups on the colour
    classes, its ORDER is prod (class size)! with no enumeration. That is the "count, don't
    enumerate" idea, already implemented, covering the antichain case that motivated the attempt.

ONE CHANGE HERE, and it is the only reason this file exists rather than a call into trd_stream:
aut_partitions() also builds the cycle-partition table that the unlabeled Burnside sum needs, via
_set_partitions - Bell(11) = 678,570 for a symmetric class. The labeled count never looks at a
cycle. aut_order_gens() below is aut_partitions with every partition line removed; it returns the
order and the generators, which is all the walk and the labeled sum require.

usage: python 2026-08-26-labeled-stream.py [max_n]
"""
import os
import sys
import time
from math import comb, factorial

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import trd_stream as ts                                    # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
BFILE = os.path.join(HERE, '2026-08-26-a399095-bfile.txt')
BASELINE = os.path.join(HERE, '2026-08-26-labeled-triangle-n10.txt')

A000112 = [1, 1, 2, 5, 16, 63, 318, 2045, 16999, 183231, 2567284, 46749427]
A006905 = [1, 2, 13, 171, 3994, 154303, 9415189, 878222530, 122207703623,
           24890747921947, 7307450299510288, 3053521546333103057]


def aut_order_gens(m, succ, cls):
    """(|Aut(P)|, generators). trd_stream.aut_partitions with the partition work stripped."""
    if ts._full_symmetric(m, succ, cls):
        order = 1
        for members in cls:
            order *= ts.FACT[len(members)]
        ident = list(range(m))
        gens = []
        for members in cls:                     # adjacent transpositions generate
            for i in range(len(members) - 1):
                g = ident[:]
                g[members[i]], g[members[i + 1]] = members[i + 1], members[i]
                gens.append(tuple(g))
        return order, gens

    order, gens = 0, []
    ident = tuple(range(m))
    for p in ts._within_class_perms(m, cls):
        if not ts._is_aut(m, succ, p):
            continue
        order += 1
        if p != ident:
            gens.append(p)
    assert order, 'automorphism group came out empty'
    return order, gens


def compositions(n, m):
    if m == 1:
        if n >= 1:
            yield (n,)
        return
    for first in range(1, n - m + 2):
        for rest in compositions(n - first, m - 1):
            yield (first,) + rest


def run(max_n, progress=None):
    rows = {n: {} for n in range(1, max_n + 1)}
    poset_count = {m: 0 for m in range(1, max_n + 1)}
    fact = [factorial(i) for i in range(max_n + 1)]
    t0 = time.time()

    def visit(m, succ, cls, order, gens):
        poset_count[m] += 1
        if progress and m == max_n and poset_count[m] % progress == 0:
            print('    m=%d  %d posets  %.0fs' % (m, poset_count[m], time.time() - t0), flush=True)
        edges = [(v, w) for v in range(m) for w in range(m) if succ[v] >> w & 1]

        for n in range(m, max_n + 1):
            tally = {}
            nf = fact[n]
            for sizes in compositions(n, m):
                base, singles, weight = 0, 0, nf
                for s in sizes:
                    if s >= 2:
                        base += s * s
                    else:
                        singles += 1
                    weight //= fact[s]
                for v, w in edges:
                    base += sizes[v] * sizes[w]
                for j in range(singles + 1):
                    kk = base + j
                    tally[kk] = tally.get(kk, 0) + weight * comb(singles, j)
            dst = rows[n]
            for kk, v in tally.items():
                assert v % order == 0, 'non-free orbit at n=%d m=%d k=%d' % (n, m, kk)
                dst[kk] = dst.get(kk, 0) + v // order

        if m == max_n:
            return
        seen = set()
        for d in ts.down_sets(m, succ):
            if d in seen:
                continue
            orbit = ts.down_set_orbit(d, gens)
            seen |= orbit
            if d != min(orbit):
                continue
            child = ts.augment(m, succ, d)
            m1 = m + 1
            cls_c = ts.classes_of(m1, ts.refine(m1, child, ts.predecessors(m1, child)))
            order_c, gens_c = aut_order_gens(m1, child, cls_c)
            if not ts.child_is_canonical(m1, child, cls_c, gens_c):
                continue
            visit(m1, child, cls_c, order_c, gens_c)

    succ0 = (0,)
    cls0 = ts.classes_of(1, ts.refine(1, succ0, [0]))
    order0, gens0 = aut_order_gens(1, succ0, cls0)
    visit(1, succ0, cls0, order0, gens0)
    return rows, poset_count


def main():
    max_n = int(sys.argv[1]) if len(sys.argv) > 1 else 9
    progress = int(sys.argv[2]) if len(sys.argv) > 2 else 250000
    print('labeled triangle on the streaming catalogue, to n = %d' % max_n, flush=True)
    t0 = time.time()
    rows, counts = run(max_n, progress=progress)
    print('  walk finished in %.1fs' % (time.time() - t0), flush=True)

    print('\nORACLE catalogue sizes, against A000112')
    for m in sorted(counts):
        want = A000112[m] if m < len(A000112) else None
        got = counts[m]
        print('  m=%-3d %12d  A000112 %-12s %s'
              % (m, got, want, 'OK' if want == got else ('NEW' if want is None else 'MISMATCH')))
        assert want is None or got == want, 'catalogue wrong at m=%d' % m

    print('\nORACLE row sums, against A006905')
    for n in sorted(rows):
        got = sum(rows[n].values())
        want = A006905[n] if n < len(A006905) else None
        print('  n=%-3d %24d  A006905 %-24s %s'
              % (n, got, want, 'OK' if want == got else ('NEW' if want is None else 'MISMATCH')))
        assert want is None or got == want, 'row sum wrong at n=%d' % n

    if os.path.exists(BFILE):
        pub = [int(l.split()[1]) for l in open(BFILE) if l.strip() and not l.startswith('#')]
        idx, bad = 0, 0
        for n in range(0, max_n + 1):
            for kk in range(n * n + 1):
                if idx >= len(pub):
                    break
                got = 1 if n == 0 and kk == 0 else rows.get(n, {}).get(kk, 0)
                if got != pub[idx]:
                    bad += 1
                idx += 1
            if idx >= len(pub):
                break
        print('\nORACLE b-file: %d checked, %d mismatches' % (idx, bad))
        assert bad == 0, 'b-file oracle FAILED'

    if os.path.exists(BASELINE):
        txt = open(BASELINE).read()
        bad = 0
        for n in (8, 9, 10):
            tag = 'row %d (' % n
            if n > max_n or tag not in txt:
                continue
            old = [int(x) for x in txt.split(tag)[1].split('\n')[1].strip().split(',')]
            new = [rows[n].get(kk, 0) for kk in range(n * n + 1)]
            if old != new:
                bad += 1
                print('  REGRESSION at row %d' % n)
        print('ORACLE regression vs 2026-08-26 rows 8-10: %s' % ('OK' if bad == 0 else 'FAILED'))
        assert bad == 0, 'regression FAILED'

    for n in range(11, max_n + 1):
        row = [rows[n].get(kk, 0) for kk in range(n * n + 1)]
        print('\nrow %d (%d entries):' % (n, len(row)))
        print(', '.join(str(v) for v in row))


if __name__ == '__main__':
    main()
