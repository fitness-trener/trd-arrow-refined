"""Identification test for Egri-Nagy's arrow-type semigroupoids.

Their Table 1 caption calls them "transitively closed directed graphs with no
parallel edges up to isomorphism". This brute-forces that class directly and
compares against their published column sums, then against OEIS A091073
(transitive relations on n unlabeled points, Pfeiffer 2004 / Brinkmann-McKay 2005).

Claim under test: arrow-type semigroupoids on n objects = transitive digraphs on
n unlabelled vertices (loops allowed) with NO ISOLATED VERTICES, which equals
A091073(n) - A091073(n-1).
"""
from itertools import permutations

PUBLISHED = {1: 1, 2: 6, 3: 31, 4: 203, 5: 1653, 6: 17156, 7: 227844}
A091073 = [1, 2, 8, 39, 242, 1895, 19051, 246895, 4145108, 90325655, 2555630036,
           93810648902, 4461086120602, 274339212258846, 21775814889230580,
           2226876304576948549]


def transitive(adj, n):
    for x in range(n):
        for y in range(n):
            if adj >> (x * n + y) & 1:
                for z in range(n):
                    if (adj >> (y * n + z) & 1) and not (adj >> (x * n + z) & 1):
                        return False
    return True


def no_isolated(adj, n):
    for v in range(n):
        row = any(adj >> (v * n + y) & 1 for y in range(n))
        col = any(adj >> (x * n + v) & 1 for x in range(n))
        if not (row or col):
            return False
    return True


def canon(adj, n):
    best = None
    for p in permutations(range(n)):
        m = 0
        for x in range(n):
            for y in range(n):
                if adj >> (x * n + y) & 1:
                    m |= 1 << (p[x] * n + p[y])
        if best is None or m < best:
            best = m
    return best


def count(n, require_no_isolated):
    seen = set()
    for adj in range(1 << (n * n)):
        if not transitive(adj, n):
            continue
        if require_no_isolated and not no_isolated(adj, n):
            continue
        seen.add(canon(adj, n))
    return len(seen)


def main():
    print('%3s %14s %14s %14s' % ('n', 'brute(no iso)', 'published', 'A091073 diff'))
    for n in (1, 2, 3, 4):
        b = count(n, True)
        d = A091073[n] - A091073[n - 1]
        print('%3d %14d %14d %14d  %s'
              % (n, b, PUBLISHED[n], d, 'OK' if b == PUBLISHED[n] == d else 'MISMATCH'))
        assert b == PUBLISHED[n] == d, 'identification fails at n=%d' % n

    print('\nall transitive relations up to iso (should equal A091073):')
    for n in (1, 2, 3):
        b = count(n, False)
        print('  n=%d  brute %d  A091073 %d  %s' % (n, b, A091073[n],
                                                    'OK' if b == A091073[n] else 'MISMATCH'))
        assert b == A091073[n]

    print('\nfirst differences of A091073 = arrow-type semigroupoids by objects:')
    for n in range(1, len(A091073)):
        d = A091073[n] - A091073[n - 1]
        mark = 'published' if n in PUBLISHED else 'NOT IN THEIR TABLE'
        flag = ' <- matches' if PUBLISHED.get(n) == d else ''
        print('  %2d objects: %-22d %s%s' % (n, d, mark, flag))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
