"""Labeled transitive relations by (points, ordered pairs) - extending OEIS A399095.

A399095 (Andrew Howroyd, 19 Aug 2026) is the triangle T(n,k) = number of transitive
relations on n LABELED points with exactly k ordered pairs. Its b-file stops at row 7.

This computes it from the Brinkmann-McKay structure theorem already implemented in
trd_bivariate.py, reusing that file's poset catalogue and automorphism groups.

METHOD (pre-registered in wiki/concepts/plan-a399095-labeled-extension.md).

A transitive relation on a set is: a poset P of strong components, a size n_v >= 1 per
component, and a loop / no-loop choice per size-1 component. Arrows:

    k = sum_{n_v >= 2} n_v^2  +  sum_{(v,w) in E(P)} n_v * n_w  +  (loops chosen)

For the LABELED count, fix any labelling of P and assign actual blocks of [n]:

    T(n,k) = sum_m sum_{P on m pts} (1/|Aut(P)|) sum_{sizes,loops: arrows=k} n!/prod n_v!

The division is an orbit count of a FREE action: two block-assignments give the same
labeled relation only if some g in Aut(P) has B_v = B_{g(v)} for all v, and the blocks
are disjoint and non-empty, so g fixes every vertex. Freeness is what makes every
per-k tally divisible by |Aut(P)| - asserted below rather than assumed.

Loop choices need no enumeration: each size-1 component contributes an independent
optional +1 arrow, so a poset with s of them contributes a factor (1+x)^s.

usage: python 2026-08-26-labeled-triangle.py [max_n]
"""
import os
import sys
from math import comb, factorial

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from trd_bivariate import gen_posets            # noqa: E402

# Oracle 2: row sums. OEIS A006905, read 2026-08-26.
A006905 = [1, 2, 13, 171, 3994, 154303, 9415189, 878222530, 122207703623,
           24890747921947, 7307450299510288, 3053521546333103057]

BFILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                     '2026-08-26-a399095-bfile.txt')


def compositions(n, m):
    """All tuples of m positive integers summing to n."""
    if m == 1:
        if n >= 1:
            yield (n,)
        return
    for first in range(1, n - m + 2):
        for rest in compositions(n - first, m - 1):
            yield (first,) + rest


def labeled_rows(max_n):
    levels = gen_posets(max_n)
    rows = {0: {0: 1}}                       # the empty relation on 0 points
    for n in range(1, max_n + 1):
        row = {}
        nfact = factorial(n)
        for m in range(1, n + 1):
            edges = [(v, w) for succ in [None] for v in [] for w in []]  # placeholder
            for succ, aut in zip(levels[m], levels['aut'][m]):
                order = len(aut)
                edges = [(v, w) for v in range(m) for w in range(m)
                         if succ[v] >> w & 1]
                tally = {}
                for sizes in compositions(n, m):
                    base = sum(s * s for s in sizes if s >= 2)
                    base += sum(sizes[v] * sizes[w] for v, w in edges)
                    singles = sum(1 for s in sizes if s == 1)
                    weight = nfact
                    for s in sizes:
                        weight //= factorial(s)
                    for j in range(singles + 1):
                        kk = base + j
                        tally[kk] = tally.get(kk, 0) + weight * comb(singles, j)
                for kk, v in tally.items():
                    assert v % order == 0, (
                        'non-free orbit count at n=%d, m=%d, k=%d' % (n, m, kk))
                    row[kk] = row.get(kk, 0) + v // order
        rows[n] = row
        print('  n=%-2d done, row sum %d' % (n, sum(row.values())), flush=True)
    return rows


def main():
    max_n = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    print('labeled transitive-relation triangle to n = %d' % max_n)
    rows = labeled_rows(max_n)

    print('\nORACLE 2 - row sums against A006905')
    for n in sorted(rows):
        got = sum(rows[n].values())
        want = A006905[n] if n < len(A006905) else None
        flag = 'OK' if want == got else ('NEW' if want is None else 'MISMATCH')
        print('  n=%-2d %22d  A006905 %-22s %s' % (n, got, want, flag))
        assert want is None or got == want, 'row sum mismatch at n=%d' % n

    print('\nORACLE 1 - every published A399095 b-file value')
    if not os.path.exists(BFILE):
        print('  b-file not present at %s - SKIPPED' % BFILE)
    else:
        published = []
        with open(BFILE) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                published.append(int(line.split()[1]))
        idx, checked, bad = 0, 0, 0
        for n in range(0, max_n + 1):
            width = n * n + 1
            for kk in range(width):
                if idx >= len(published):
                    break
                got = rows[n].get(kk, 0)
                if got != published[idx]:
                    bad += 1
                    print('  MISMATCH n=%d k=%d: got %d, published %d'
                          % (n, kk, got, published[idx]))
                idx += 1
                checked += 1
            if idx >= len(published):
                break
        print('  %d of %d published values checked, %d mismatches'
              % (checked, len(published), bad))
        assert bad == 0, 'ORACLE 1 FAILED - method is void'

    print('\nNEW ROWS (not in the A399095 b-file)')
    for n in range(8, max_n + 1):
        row = [rows[n].get(kk, 0) for kk in range(n * n + 1)]
        print('  row %d (%d entries):' % (n, len(row)))
        print('   ', ', '.join(str(v) for v in row))


if __name__ == '__main__':
    main()
