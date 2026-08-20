"""Validate the shipped table against the labeled transitive-relation literature.

Independent check, added in v1.1.0. It uses no enumeration and no code from this
repository's generators, so it is a genuine outside test of `data/trd-table-n10.csv`
rather than a self-consistency check.

The bridge
----------
OEIS carries a family counting LABELED transitive relations on an n-set by exact
number of ordered pairs, each with a closed form

    a(n) = sum_m  c_{k,m} * C(n,m)          (k = number of ordered pairs)

    k=2 A349919   k=3 A349927   k=4 A349849   k=5 A348634
    k=6 A393216   k=7 A392460   k=8 A393925   k=9 A393974

The C(n,m) factor chooses which m of the n points are active (non-isolated), so
`c_{k,m}` is the labeled count of structures whose active set is exactly those m
points. This repository's cell T(objects=m, arrows=k) counts the same structures up
to isomorphism. Same family, same slicing, same no-isolated-point convention -- one
labeled, one not.

Two consequences are tested here.

1. Bracket.  T(m,k) is a sum over isomorphism classes and c_{k,m} = sum over the
   same classes of m!/|Aut|, so

       c_{k,m} / m!   <=   T(m,k)   <=   c_{k,m}

   with the lower equality iff every class has the full symmetric group as its
   automorphism group and the upper iff every class is rigid.

2. Support.  A structure with k arcs and no isolated object spans at most 2k
   objects (k disjoint arrows), so column k of the labeled closed form must run out
   at m = 2k, and the object range of arrow-count k must agree between the two
   sources.

Coefficients below are transcribed from the %F lines of the OEIS entries, read in
text format. They are inputs to this check, not outputs of it.
"""

import csv
import os
import sys

MALA = {
    2: {2: 5, 3: 12, 4: 12},
    3: {2: 2, 3: 37, 4: 116, 5: 180, 6: 120},
    4: {2: 1, 3: 42, 4: 375, 5: 1450, 6: 2940, 7: 3360, 8: 1680},
    5: {3: 27, 4: 660, 5: 5201, 6: 21822, 7: 54600, 8: 84000, 9: 75600, 10: 30240},
    6: {3: 12, 4: 758, 5: 11600, 6: 87943, 7: 393960, 8: 1137752, 9: 2172240,
        10: 2696400, 11: 1995840, 12: 665280},
    7: {3: 6, 4: 644, 5: 18550, 6: 237210, 7: 1748587, 8: 8279224, 9: 26645976,
        10: 59729040, 11: 92862000, 12: 96465600, 13: 60540480, 14: 17297280},
    8: {4: 429, 5: 23120, 6: 478635, 7: 5507586, 8: 40073168, 9: 198343134,
        10: 696045600, 11: 1767981600, 12: 3257709840, 13: 4276752480,
        14: 3814050240, 15: 2075673600, 16: 518918400},
    9: {3: 1, 4: 216, 5: 23325, 6: 772880, 7: 13381781, 8: 143627512,
        9: 1039241647, 10: 5338438070, 11: 20105261220, 12: 56501068800,
        13: 119067989040, 14: 186575669280, 15: 211891680000, 16: 165361996800,
        17: 79394515200, 18: 17643225600},
}

OEIS_ID = {2: "A349919", 3: "A349927", 4: "A349849", 5: "A348634",
           6: "A393216", 7: "A392460", 8: "A393925", 9: "A393974"}

MAX_OBJECTS = 10          # what this repository's table covers


def factorial(n):
    f = 1
    for i in range(2, n + 1):
        f *= i
    return f


def load(path):
    cells = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            m = int(row["objects"])
            k = int(row["arrows"])
            v = int(row["arrow_type_semigroupoids"])
            if v:
                cells[(m, k)] = v
    return cells


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(here, os.pardir, "data", "trd-table-n10.csv")
    cells = load(path)
    print(f"loaded {len(cells)} non-zero cells from data/trd-table-n10.csv")
    print()

    print("1. bracket   c/m! <= T(m,k) <= c")
    checked = 0
    for k in sorted(MALA):
        row = []
        for m, c in sorted(MALA[k].items()):
            if m > MAX_OBJECTS:
                continue                      # outside this table's object range
            t = cells.get((m, k), 0)
            assert t > 0, f"table has no cell for objects={m}, arrows={k}"
            lo = c / factorial(m)
            assert lo <= t <= c, \
                f"bracket violated at objects={m}, arrows={k}: {lo} <= {t} <= {c} is false"
            row.append(m)
            checked += 1
        print(f"   arrows={k:>2} ({OEIS_ID[k]}): objects {row[0]}..{row[-1]} all inside bracket")
    print(f"   {checked} cells validated against the labeled closed forms")
    print()

    print("2. support   labeled support must end at 2k, and match this table's object range")
    for k in sorted(MALA):
        hi = max(MALA[k])
        assert hi == 2 * k, f"labeled support for k={k} ends at {hi}, expected {2*k}"
        lo_lab = min(MALA[k])
        lo_tab = min(m for (m, kk) in cells if kk == k)
        assert lo_lab == lo_tab, \
            f"support start mismatch at k={k}: labeled {lo_lab}, table {lo_tab}"
        state = "complete here" if hi <= MAX_OBJECTS else \
                f"objects {MAX_OBJECTS+1}..{hi} beyond this table"
        print(f"   arrows={k:>2}: objects {lo_lab}..{hi}   start matches; {state}")
    print()

    print("all assertions passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
