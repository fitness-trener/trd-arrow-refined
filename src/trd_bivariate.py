"""Arrow-refined enumeration of arrow-type semigroupoids, via the poset catalogue.

Egri-Nagy & Nehaniv (arXiv:2509.00837) Table 1 counts arrow-type semigroupoids by
objects AND arrows, with the larger entries unknown. q19 established that the class
is Brinkmann-McKay's transitive relation digraphs (TRDs), whose count by objects is
A091073 (known to 15 points, 2005) - but their refinement is by strong components,
not arrows, so the arrow-refined table is not in the literature.

This computes it from their structure theorem instead of by searching semigroupoids.

STRUCTURE (Brinkmann-McKay, section 2). Contracting strong components of a TRD gives
a poset P; each component is a directed clique. Conversely a TRD on n points is:

    a poset P on m points
  + a size n_v >= 1 per point of P, with sum n_v = n
  + for each SINGLETON component, a loop / no-loop choice
    (components of size >= 2 necessarily carry all loops and both directions)

and two such data give isomorphic TRDs exactly when related by an automorphism of P.

ARROW COUNT in closed form:

    arrows = sum_v f(n_v) + sum_{(v,w) in E(P)} n_v * n_w,
    f(k) = k^2 for k >= 2,  f(1) = 1 if the loop is chosen else 0

Counting orbits per arrow count is weighted Burnside over Aut(P): a g-fixed
assignment is constant on the cycles of g, and sum_v n_v = n bounds the number of
size tuples to the compositions of n, so the inner enumeration is tiny.

ISOLATED OBJECTS. Egri-Nagy's objects are arrow-types that occur, so his table
excludes TRDs with an isolated point (a loopless singleton component isolated in P).
Both variants are computed: the unfiltered totals must reproduce A091073, and the
filtered ones his column sums.

usage: python trd_bivariate.py [max_n]
"""
import sys
from itertools import permutations

# Published values this must reproduce before any new row is believed.
A091073 = [1, 2, 8, 39, 242, 1895, 19051, 246895, 4145108, 90325655, 2555630036]
EN_BY_OBJECTS = {1: 1, 2: 6, 3: 31, 4: 203, 5: 1653, 6: 17156, 7: 227844}
# Row sums read off his Table 1. The arrows=1 entry was first transcribed as 1 from a
# mangled PDF render; it is 2 (the lone loop on one object, and x->y on two), which is
# what this program computes and what the row's two cells add to.
EN_BY_ARROWS = {1: 2, 2: 7, 3: 21, 4: 70, 5: 218, 6: 721, 7: 2360}
NUM_POSETS = {1: 1, 2: 2, 3: 5, 4: 16, 5: 63, 6: 318, 7: 2045, 8: 16999}


# ----------------------------------------------------------------- posets
def refine(m, succ, pred):
    """Colour points by iterated degree/neighbour-colour refinement."""
    col = [(bin(succ[v]).count('1'), bin(pred[v]).count('1')) for v in range(m)]
    for _ in range(m):
        new = []
        for v in range(m):
            s = tuple(sorted(col[w] for w in range(m) if succ[v] >> w & 1))
            p = tuple(sorted(col[w] for w in range(m) if pred[v] >> w & 1))
            new.append((col[v], s, p))
        codes = {c: i for i, c in enumerate(sorted(set(new)))}
        new = [codes[c] for c in new]
        if new == col:
            break
        col = new
    return col


def encode(m, succ, perm):
    """Adjacency bitstring of the relabelled poset, for lexicographic comparison."""
    out = 0
    for v in range(m):
        for w in range(m):
            if succ[v] >> w & 1:
                out |= 1 << (perm[v] * m + perm[w])
    return out


def _class_perms(m, col):
    """All permutations respecting the colour classes."""
    classes = {}
    for v in range(m):
        classes.setdefault(col[v], []).append(v)
    keys = sorted(classes)
    slots, pools = [], []
    nxt = 0
    for k in keys:
        members = classes[k]
        slots.append(members)
        pools.append(list(range(nxt, nxt + len(members))))
        nxt += len(members)

    def rec(i, perm):
        if i == len(slots):
            yield tuple(perm)
            return
        for order in permutations(pools[i]):
            for v, t in zip(slots[i], order):
                perm[v] = t
            yield from rec(i + 1, perm)

    return rec(0, [0] * m)


def _within_class_perms(m, col):
    """Permutations that permute each colour class onto ITSELF.

    Automorphisms live here. _class_perms renumbers classes into consecutive slots,
    which is what a canonical form needs but does not even contain the identity, so
    it cannot be used to detect automorphisms.
    """
    classes = {}
    for v in range(m):
        classes.setdefault(col[v], []).append(v)
    keys = sorted(classes)

    def rec(i, perm):
        if i == len(keys):
            yield tuple(perm)
            return
        members = classes[keys[i]]
        for order in permutations(members):
            for v, t in zip(members, order):
                perm[v] = t
            yield from rec(i + 1, perm)

    return rec(0, [0] * m)


def canon_and_aut(m, succ):
    """Canonical encoding and the full automorphism group (as permutation tuples)."""
    pred = [0] * m
    for v in range(m):
        for w in range(m):
            if succ[v] >> w & 1:
                pred[w] |= 1 << v
    col = refine(m, succ, pred)
    best = None
    for p in _class_perms(m, col):
        e = encode(m, succ, p)
        if best is None or e < best:
            best = e
    ident = encode(m, succ, tuple(range(m)))
    aut = [p for p in _within_class_perms(m, col) if encode(m, succ, p) == ident]
    assert aut, 'automorphism group came out empty - the identity must always be in it'
    return best, aut


def ideals(m, succ):
    """Down-sets of the poset, as bitmasks (a new maximal point may sit above any)."""
    pred = [0] * m
    for v in range(m):
        for w in range(m):
            if succ[v] >> w & 1:
                pred[w] |= 1 << v
    out = []
    for s in range(1 << m):
        ok = True
        for v in range(m):
            if s >> v & 1 and pred[v] & ~s & ((1 << m) - 1):
                ok = False
                break
        if ok:
            out.append(s)
    return out


def gen_posets(maxm):
    """Unlabelled posets on 1..maxm points, by adding one maximal point at a time."""
    levels = {0: [()]}
    for m in range(1, maxm + 1):
        seen = {}
        for succ in levels[m - 1]:
            for d in ideals(m - 1, succ):
                new = list(succ) + [0]
                for v in range(m - 1):          # new point x = m-1 sits above d
                    if d >> v & 1:
                        new[v] |= 1 << (m - 1)
                new = tuple(new)
                c, aut = canon_and_aut(m, new)
                if c not in seen:
                    seen[c] = (new, aut)
        levels[m] = [v[0] for v in seen.values()]
        levels.setdefault('aut', {})[m] = {k: v[1] for k, v in seen.items()}
        levels[m] = [seen[k][0] for k in seen]
        levels['aut'][m] = [seen[k][1] for k in seen]
    return levels


# ------------------------------------------------------- weighted Burnside
def cycles_of(perm):
    m = len(perm)
    seen = [False] * m
    out = []
    for v in range(m):
        if seen[v]:
            continue
        c = []
        x = v
        while not seen[x]:
            seen[x] = True
            c.append(x)
            x = perm[x]
        out.append(c)
    return out


def contribute(m, succ, aut, n, tally_all, tally_noiso):
    """Add this poset's TRDs on exactly n points to the two arrow-count tallies."""
    isolated = [succ[v] == 0 and all(not (succ[w] >> v & 1) for w in range(m))
                for v in range(m)]
    edges = [(v, w) for v in range(m) for w in range(m) if succ[v] >> w & 1]

    for g in aut:
        cyc = cycles_of(g)
        lens = [len(c) for c in cyc]
        cyc_of = [0] * m
        for i, c in enumerate(cyc):
            for v in c:
                cyc_of[v] = i

        def rec(i, left, sizes):
            if i == len(cyc):
                if left:
                    return
                size_of = [sizes[cyc_of[v]] for v in range(m)]
                cross = sum(size_of[v] * size_of[w] for v, w in edges)
                base = sum(size_of[v] ** 2 for v in range(m) if size_of[v] >= 2)
                # loop choices: one per cycle whose size is 1, constant on the cycle
                singles = [j for j, s in enumerate(sizes) if s == 1]
                for mask in range(1 << len(singles)):
                    arrows = base + cross
                    ok_noiso = True
                    for bit, j in enumerate(singles):
                        if mask >> bit & 1:
                            arrows += lens[j]
                        else:
                            if any(isolated[v] for v in cyc[j]):
                                ok_noiso = False
                    tally_all[arrows] = tally_all.get(arrows, 0) + 1
                    if ok_noiso:
                        tally_noiso[arrows] = tally_noiso.get(arrows, 0) + 1
                return
            # size a >= 1 for this cycle, costing lens[i]*a points
            a = 1
            while lens[i] * a <= left:
                rec(i + 1, left - lens[i] * a, sizes + [a])
                a += 1

        rec(0, n, [])


def enumerate_trds(max_n):
    levels = gen_posets(max_n)
    for m in range(1, max_n + 1):
        got, want = len(levels[m]), NUM_POSETS.get(m)
        assert want is None or got == want, 'posets on %d points: got %d, want %d' % (m, got, want)
        print('  posets on %d points: %d' % (m, got), flush=True)

    all_t, noiso_t = {}, {}
    for n in range(1, max_n + 1):
        ta, tn = {}, {}
        for m in range(1, n + 1):
            for succ, aut in zip(levels[m], levels['aut'][m]):
                raw_a, raw_n = {}, {}
                contribute(m, succ, aut, n, raw_a, raw_n)
                order = len(aut)
                for d, src in ((ta, raw_a), (tn, raw_n)):
                    for k, v in src.items():
                        assert v % order == 0, 'Burnside non-integral at n=%d' % n
                        d[k] = d.get(k, 0) + v // order
        all_t[n], noiso_t[n] = ta, tn
        print('  n=%d done (%d TRDs, %d without isolated points)'
              % (n, sum(ta.values()), sum(tn.values())), flush=True)
    return all_t, noiso_t


def main():
    max_n = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    print('generating posets and enumerating TRDs to n = %d' % max_n)
    all_t, noiso_t = enumerate_trds(max_n)

    print('\nCHECK 1 - all TRDs by points, against A091073')
    for n in sorted(all_t):
        got, want = sum(all_t[n].values()), A091073[n] if n < len(A091073) else None
        print('  n=%-3d %16d  A091073 %-16s %s'
              % (n, got, want, 'OK' if want == got else ('?' if want is None else 'MISMATCH')))
        assert want is None or got == want

    print('\nCHECK 2 - no isolated points, against Egri-Nagy Table 1 column sums')
    for n in sorted(noiso_t):
        got, want = sum(noiso_t[n].values()), EN_BY_OBJECTS.get(n)
        print('  n=%-3d %16d  published %-10s %s'
              % (n, got, want, 'OK' if want == got else ('NEW' if want is None else 'MISMATCH')))
        assert want is None or got == want

    print('\nCHECK 3 - row sums by arrow count, against his published row sums')
    by_arrows = {}
    for n in noiso_t:
        for a, c in noiso_t[n].items():
            by_arrows[a] = by_arrows.get(a, 0) + c
    for a in sorted(by_arrows)[:8]:
        want = EN_BY_ARROWS.get(a)
        # With no isolated points every object touches an arrow, so a arrows reach at
        # most 2a objects: the row is complete only once max_n >= 2a.
        complete = 2 * a <= max_n
        print('  arrows=%-3d %12d  published %-8s %s'
              % (a, by_arrows[a], want, '' if not complete else
                 ('OK' if want == by_arrows[a] else 'MISMATCH/partial')))

    print('\nTABLE - arrow-type semigroupoids by objects (rows) and arrows (columns)')
    for n in sorted(noiso_t):
        row = noiso_t[n]
        cells = ' '.join('%d:%d' % (a, row[a]) for a in sorted(row))
        print('  n=%-3d total %-12d %s' % (n, sum(row.values()), cells))

    out = 'trd-table-n%d.csv' % max_n
    with open(out, 'w', encoding='utf-8') as f:
        f.write('objects,arrows,arrow_type_semigroupoids,all_trds\n')
        for n in sorted(noiso_t):
            for a in sorted(set(noiso_t[n]) | set(all_t[n])):
                f.write('%d,%d,%d,%d\n' % (n, a, noiso_t[n].get(a, 0), all_t[n].get(a, 0)))
    print('\nwrote %s' % out)
    return 0


if __name__ == '__main__':
    sys.exit(main())
