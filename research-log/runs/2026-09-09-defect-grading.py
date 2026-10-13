"""Regrade the A399030 component decomposition by DEFECT d = 2e - v.

Predictions are pre-registered in wiki/concepts/plan-defect-grading.md and are NOT restated here.

C_d = # connected arrow-type semigroupoid classes of defect d = sum_{2e-v=d} D(e,v)
D_N = # classes of defect N     = [x^N] prod_{d>0} (1-x^d)^{-C_d}

Connectivity forces v <= e+1, hence e ranges over ceil((d+1)/2) .. d+1, so with the existing
MAX_E = 6 the value C_d is COMPLETE only for d <= 5.

Reuses connected_classes() from 2026-08-30-low-arrow-cells.py verbatim - that function was
oracle-checked on 60 published cells on 2026-08-30.

usage: python 2026-09-09-defect-grading.py
"""
import importlib.util
import io
import os
import time
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))

spec = importlib.util.spec_from_file_location(
    "lowarrow", os.path.join(HERE, "2026-08-30-low-arrow-cells.py"))
lowarrow = importlib.util.module_from_spec(spec)
spec.loader.exec_module(lowarrow)

MAX_E = lowarrow.MAX_E
D_COMPLETE_TO = MAX_E - 1          # C_d complete iff d+1 <= MAX_E


def euler(C, N):
    """prod_{d>=1} (1-x^d)^{-C[d]} up to x^N."""
    P = [1] + [0] * N
    for d in range(1, N + 1):
        S = [0] * (N + 1)
        j = 0
        while d * j <= N:
            S[d * j] = comb(C.get(d, 0) + j - 1, j) if C.get(d, 0) else (1 if j == 0 else 0)
            j += 1
        Q = [0] * (N + 1)
        for a in range(N + 1):
            if not P[a]:
                continue
            for b in range(N + 1 - a):
                if S[b]:
                    Q[a + b] += P[a] * S[b]
        P = Q
    return P


def published_rows():
    """rows[n][k] from the shipped bivariate table (same file the 2026-08-30 oracle used)."""
    return lowarrow.shipped_table()


def main():
    print("connected classes, e <= %d (reused from the 2026-08-30 run)" % MAX_E, flush=True)
    t0 = time.time()
    cache = os.path.join(HERE, "2026-09-09-connected-classes.json")
    if os.path.exists(cache):
        import json
        D = {tuple(map(int, k.split(","))): v for k, v in json.load(open(cache)).items()}
        print("  [loaded from %s]" % os.path.basename(cache), flush=True)
    else:
        D = lowarrow.connected_classes()
        import json
        json.dump({"%d,%d" % k: v for k, v in D.items()}, open(cache, "w"))
    print("  %.1fs; %d types, %d classes\n" % (time.time() - t0, len(D), sum(D.values())),
          flush=True)

    print("D(e,v) table")
    for (e, v) in sorted(D):
        print("  e=%d v=%d defect=%-3d D=%d" % (e, v, 2 * e - v, D[(e, v)]))

    C = {}
    for (e, v), cnt in D.items():
        C[2 * e - v] = C.get(2 * e - v, 0) + cnt
    print("\nC_d = # connected classes of defect d")
    for d in sorted(C):
        flag = "complete" if 0 < d <= D_COMPLETE_TO else ("d<=0: not a positive defect"
                                                          if d <= 0 else "INCOMPLETE (needs e=%d)" % (d + 1))
        print("  C_%-3d = %-6d  [%s]" % (d, C[d], flag))

    # Spitz Prop. 2 analogue: the unique connected class of minimum defect is the single arrow.
    assert min(C) == 0 and C[0] == 1, "defect-0 connected classes: %r" % C.get(0)
    assert D[(1, 2)] == 1, "the defect-0 class should be e=1, v=2 (one arrow, two points)"

    Cpos = {d: c for d, c in C.items() if d <= D_COMPLETE_TO}
    DN = euler(Cpos, D_COMPLETE_TO)
    print("\nD_N = # classes of defect N   (trustworthy for N <= %d)" % D_COMPLETE_TO)
    for N, v in enumerate(DN):
        print("  D_%d = %d" % (N, v))

    rows = published_rows()
    if rows is None:
        print("\n[b-file not found; P2 check skipped]")
        return
    print("\nP2 - defect diagonal of the published A399030 b-file must be constant = D_N")
    bad = 0
    for N in range(0, D_COMPLETE_TO + 1):
        seen = []
        for k in range(max(1, 2 * N), 60):
            n = 2 * k - N
            if n < 1 or n not in rows or k not in rows[n]:
                continue
            seen.append((n, k, rows[n][k]))
        if not seen:
            print("  D_%d : no published cell with k >= 2N" % N)
            continue
        vals = {t[2] for t in seen}
        ok = (vals == {DN[N]})
        bad += 0 if ok else 1
        print("  D_%d = %-6d | b-file %s -> %s"
              % (N, DN[N], ", ".join("T(%d,%d)=%d" % t for t in seen), "OK" if ok else "MISMATCH"))
    print("\nP2 verdict: %s" % ("PASS" if bad == 0 else "FAIL (%d)" % bad))


if __name__ == "__main__":
    main()
