"""Streaming arrow-refined TRD enumeration — constant memory in the poset count.

Same mathematics as trd_bivariate.py (Brinkmann-McKay structure theorem, weighted
Burnside over Aut(P)); two changes, both forced by probe P's measurements.

1. STREAMING. trd_bivariate builds every level and holds it. Probe P measured
   1.36 KB per poset, which is ~64 GB at m=11 against 8.35 GB physical. Here posets
   are visited depth-first by McKay canonical augmentation: a child is accepted only
   if the point just added lies in the Aut-orbit of the canonically-chosen maximal
   point, so each isomorphism class is reached exactly once with no global table.
   Memory is O(max_m) posets.

2. CYCLE PARTITIONS, NOT AUTOMORPHISMS. The real memory blow-up was never the poset
   count: _within_class_perms materialises all of Aut(P), and an antichain on m
   points has m! automorphisms (3.6M at m=10, 40M at m=11). But `contribute` uses g
   only through its cycles AS SETS - the fixed assignments are those constant on
   cycles, and both the arrow count and the isolated-point test read the cycle sets.
   So automorphisms are collapsed to their induced set partition, carrying the
   multiplicity prod_B (|B|-1)!. For Aut(P) = prod S_{C_i} the partitions are
   enumerated directly and Aut is never listed at all: Bell(10) = 115,975 against
   10! = 3,628,800, a 31x cut on the worst poset at m=10 and far more beyond.

The published values in trd_bivariate.py are the oracle and are re-checked here.

usage: python trd_stream.py [max_n]
"""
import sys
from itertools import permutations

A091073 = [1, 2, 8, 39, 242, 1895, 19051, 246895, 4145108, 90325655, 2555630036]
EN_BY_OBJECTS = {1: 1, 2: 6, 3: 31, 4: 203, 5: 1653, 6: 17156, 7: 227844}
A000112 = {1: 1, 2: 2, 3: 5, 4: 16, 5: 63, 6: 318, 7: 2045, 8: 16999, 9: 183231,
           10: 2567284, 11: 46749427}

FACT = [1]
for _i in range(1, 32):
    FACT.append(FACT[-1] * _i)


# ------------------------------------------------------------------ poset basics
def predecessors(m, succ):
    pred = [0] * m
    for v in range(m):
        s = succ[v]
        while s:
            b = s & -s
            pred[b.bit_length() - 1] |= 1 << v
            s ^= b
    return pred


def refine(m, succ, pred):
    """Iterated degree/neighbour-colour refinement. Same as trd_bivariate."""
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
    out = 0
    for v in range(m):
        s = succ[v]
        while s:
            b = s & -s
            w = b.bit_length() - 1
            out |= 1 << (perm[v] * m + perm[w])
            s ^= b
    return out


def classes_of(m, col):
    d = {}
    for v in range(m):
        d.setdefault(col[v], []).append(v)
    return [d[k] for k in sorted(d)]


def _class_perms(m, cls):
    """Permutations sending colour classes to consecutive label blocks."""
    slots, pools, nxt = [], [], 0
    for members in cls:
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


def _within_class_perms(m, cls):
    """Permutations mapping each colour class onto itself."""
    def rec(i, perm):
        if i == len(cls):
            yield tuple(perm)
            return
        members = cls[i]
        for order in permutations(members):
            for v, t in zip(members, order):
                perm[v] = t
            yield from rec(i + 1, perm)

    return rec(0, [0] * m)


def canon(m, succ, cls):
    """Canonical encoding and one permutation achieving it."""
    best, best_perm = None, None
    for p in _class_perms(m, cls):
        e = encode(m, succ, p)
        if best is None or e < best:
            best, best_perm = e, p
    return best, best_perm


# ------------------------------------------------- Aut as cycle partitions
def _is_aut(m, succ, perm):
    return encode(m, succ, perm) == encode(m, succ, tuple(range(m)))


def _full_symmetric(m, succ, cls):
    """Is Aut(P) the full product of symmetric groups on the colour classes?

    Refinement is isomorphism-invariant, so Aut(P) always sits inside that product.
    If in addition every within-class transposition is an automorphism, the two are
    equal - and transpositions generate, so testing pairs suffices.
    """
    ident = list(range(m))
    for members in cls:
        for i in range(len(members)):
            for j in range(i + 1, len(members)):
                p = ident[:]
                p[members[i]], p[members[j]] = members[j], members[i]
                if not _is_aut(m, succ, tuple(p)):
                    return False
    return True


def _set_partitions(items):
    """All set partitions of a list, as lists of tuples."""
    if not items:
        yield []
        return
    first, rest = items[0], items[1:]
    for part in _set_partitions(rest):
        for i in range(len(part)):
            yield part[:i] + [(first,) + part[i]] + part[i + 1:]
        yield [(first,)] + part


def _partition_key(blocks):
    return tuple(sorted(tuple(sorted(b)) for b in blocks))


def aut_partitions(m, succ, cls):
    """(order of Aut, {cycle partition: multiplicity}, generators of Aut).

    `contribute` reads g only through its cycles as sets, so the partition collapse is
    lossless: a partition into blocks B is induced by prod_B (|B|-1)! elements.

    Generators are returned separately and are NOT reconstructible from the partitions
    - a block {a,b,c} coming from a 3-cycle does not imply the transposition (a b) is
    an automorphism. Anything acting on the poset (down-set orbits, point orbits) must
    use these, not the blocks.
    """
    if _full_symmetric(m, succ, cls):
        order = 1
        for members in cls:
            order *= FACT[len(members)]
        parts = {(): 1}
        for members in cls:
            nxt = {}
            for blocks in _set_partitions(list(members)):
                mult = 1
                for b in blocks:
                    mult *= FACT[len(b) - 1]
                for key, km in parts.items():
                    k2 = key + tuple(tuple(sorted(b)) for b in blocks)
                    nxt[k2] = nxt.get(k2, 0) + km * mult
            parts = nxt
        ident = list(range(m))
        gens = []
        for members in cls:                     # adjacent transpositions generate
            for i in range(len(members) - 1):
                g = ident[:]
                g[members[i]], g[members[i + 1]] = members[i + 1], members[i]
                gens.append(tuple(g))
        return order, {tuple(sorted(k)): v for k, v in parts.items()}, gens

    parts, order, gens = {}, 0, []
    ident = tuple(range(m))
    for p in _within_class_perms(m, cls):
        if not _is_aut(m, succ, p):
            continue
        order += 1
        if p != ident:
            gens.append(p)                      # Aut is small here, keep it whole
        seen = [False] * m
        blocks = []
        for v in range(m):
            if seen[v]:
                continue
            c = []
            x = v
            while not seen[x]:
                seen[x] = True
                c.append(x)
                x = p[x]
            blocks.append(tuple(sorted(c)))
        key = _partition_key(blocks)
        parts[key] = parts.get(key, 0) + 1
    assert order, 'automorphism group came out empty'
    return order, parts, gens


def aut_orbits(m, gens):
    """Orbits of Aut on points, from real generators (union-find)."""
    par = list(range(m))

    def find(x):
        while par[x] != x:
            par[x] = par[par[x]]
            x = par[x]
        return x

    for g in gens:
        for v in range(m):
            a, b = find(v), find(g[v])
            if a != b:
                par[a] = b
    return [find(v) for v in range(m)]


def down_set_orbit(d, gens):
    """Orbit of the down-set bitmask d under the automorphism generators."""
    orbit, frontier = {d}, [d]
    while frontier:
        cur = frontier.pop()
        for g in gens:
            nxt = 0
            s = cur
            while s:
                b = s & -s
                nxt |= 1 << g[b.bit_length() - 1]
                s ^= b
            if nxt not in orbit:
                orbit.add(nxt)
                frontier.append(nxt)
    return orbit


# --------------------------------------------------------- canonical augmentation
def down_sets(m, succ):
    """Down-sets as bitmasks: a new maximal point may sit above any of them."""
    pred = predecessors(m, succ)
    full = (1 << m) - 1
    out = []
    for s in range(1 << m):
        ok = True
        for v in range(m):
            if s >> v & 1 and pred[v] & ~s & full:
                ok = False
                break
        if ok:
            out.append(s)
    return out


def augment(m, succ, d):
    """Add point m as a new maximal point sitting above the down-set d."""
    new = list(succ) + [0]
    s = d
    while s:
        b = s & -s
        new[b.bit_length() - 1] |= 1 << m
        s ^= b
    return tuple(new)


def child_is_canonical(m1, succ_child, cls_child, gens_child):
    """McKay's test: was the point just added the canonical one to have added?

    The canonical labelling picks out a distinguished maximal point (the one with the
    largest canonical label). Accept the child exactly when the newly added point -
    always the last, index m1-1 - lies in that point's Aut-orbit. Different canonical
    permutations achieving the same encoding differ by an automorphism, so the choice
    is well defined on orbits.
    """
    _, perm = canon(m1, succ_child, cls_child)
    maximal = [v for v in range(m1) if succ_child[v] == 0]
    pick = max(maximal, key=lambda v: perm[v])
    orb = aut_orbits(m1, gens_child)
    return orb[pick] == orb[m1 - 1]


# ------------------------------------------------------- weighted Burnside
def contribute(m, succ, parts, order, n, tally_all, tally_noiso):
    """Add this poset's TRDs on exactly n points to the two arrow-count tallies.

    Identical to trd_bivariate.contribute except that the outer loop runs over cycle
    partitions carrying multiplicities, rather than over automorphisms one by one.
    """
    pred = predecessors(m, succ)
    isolated = [succ[v] == 0 and pred[v] == 0 for v in range(m)]
    edges = [(v, w) for v in range(m) for w in range(m) if succ[v] >> w & 1]

    for key, mult in parts.items():
        cyc = [list(b) for b in key]
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
                base = sum(s * s for s in size_of if s >= 2)
                singles = [j for j, s in enumerate(sizes) if s == 1]
                for mask in range(1 << len(singles)):
                    arrows = base + cross
                    ok_noiso = True
                    for bit, j in enumerate(singles):
                        if mask >> bit & 1:
                            arrows += lens[j]
                        elif any(isolated[v] for v in cyc[j]):
                            ok_noiso = False
                    tally_all[arrows] = tally_all.get(arrows, 0) + mult
                    if ok_noiso:
                        tally_noiso[arrows] = tally_noiso.get(arrows, 0) + mult
                return
            a = 1
            while lens[i] * a <= left:
                rec(i + 1, left - lens[i] * a, sizes + [a])
                a += 1

        rec(0, n, [])


# ------------------------------------------------------------------ the traversal
def run(max_n, progress=None):
    all_t = {n: {} for n in range(1, max_n + 1)}
    noiso_t = {n: {} for n in range(1, max_n + 1)}
    poset_count = {m: 0 for m in range(1, max_n + 1)}

    def visit(m, succ, cls, parts, order, gens):
        poset_count[m] += 1
        if progress and poset_count[m] % progress == 0:
            print('    m=%d  %d posets' % (m, poset_count[m]), flush=True)
        for n in range(m, max_n + 1):
            ra, rn = {}, {}
            contribute(m, succ, parts, order, n, ra, rn)
            for dst, src in ((all_t[n], ra), (noiso_t[n], rn)):
                for k, v in src.items():
                    assert v % order == 0, 'Burnside non-integral at n=%d' % n
                    dst[k] = dst.get(k, 0) + v // order
        if m == max_n:
            return
        # one child per Aut(P)-orbit of down-sets, then McKay's canonical-parent test
        seen = set()
        for d in down_sets(m, succ):
            if d in seen:
                continue
            orbit = down_set_orbit(d, gens)
            seen |= orbit
            if d != min(orbit):
                continue
            child = augment(m, succ, d)
            m1 = m + 1
            cls_c = classes_of(m1, refine(m1, child, predecessors(m1, child)))
            order_c, parts_c, gens_c = aut_partitions(m1, child, cls_c)
            if not child_is_canonical(m1, child, cls_c, gens_c):
                continue
            visit(m1, child, cls_c, parts_c, order_c, gens_c)

    succ0 = (0,)
    cls0 = classes_of(1, refine(1, succ0, [0]))
    order0, parts0, gens0 = aut_partitions(1, succ0, cls0)
    visit(1, succ0, cls0, parts0, order0, gens0)
    return all_t, noiso_t, poset_count


def main():
    max_n = int(sys.argv[1]) if len(sys.argv) > 1 else 7
    print('streaming TRD enumeration to n = %d' % max_n)
    all_t, noiso_t, poset_count = run(max_n)

    print('\nposet counts (against A000112)')
    for m in sorted(poset_count):
        want = A000112.get(m)
        got = poset_count[m]
        print('  m=%-3d %12d  A000112 %-12s %s'
              % (m, got, want, 'OK' if want == got else ('?' if want is None else 'MISMATCH')))
        assert want is None or got == want

    print('\nall TRDs by points (against A091073)')
    for n in sorted(all_t):
        got = sum(all_t[n].values())
        want = A091073[n] if n < len(A091073) else None
        print('  n=%-3d %16d  A091073 %-16s %s'
              % (n, got, want, 'OK' if want == got else ('?' if want is None else 'MISMATCH')))
        assert want is None or got == want

    print('\nno isolated points (against Egri-Nagy column sums)')
    for n in sorted(noiso_t):
        got, want = sum(noiso_t[n].values()), EN_BY_OBJECTS.get(n)
        print('  n=%-3d %16d  published %-10s %s'
              % (n, got, want, 'OK' if want == got else ('NEW' if want is None else 'MISMATCH')))
        assert want is None or got == want

    out = 'trd-stream-n%d.csv' % max_n
    with open(out, 'w', encoding='utf-8') as f:
        f.write('objects,arrows,arrow_type_semigroupoids,all_trds\n')
        for n in sorted(noiso_t):
            for a in sorted(set(noiso_t[n]) | set(all_t[n])):
                f.write('%d,%d,%d,%d\n' % (n, a, noiso_t[n].get(a, 0), all_t[n].get(a, 0)))
    print('\nwrote %s' % out)
    return 0


if __name__ == '__main__':
    sys.exit(main())
