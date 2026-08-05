"""Probe P — where is the TRD ceiling: the poset catalogue, or the per-poset work?

The n=9 run reached 86,180,547 TRDs and matched A091073 exactly. n=10 has a published
check value too (A091073[10] = 2,555,630,036), so the oracle exists. The question A9
asks is whether n=10 is reachable, and which of the two phases is the wall.

Measures, per level m: wall time and peak RSS of gen_posets, then wall time of the
contribute phase for one target n. Extrapolates both to m=10 from the measured
per-poset constants rather than from a guess.

usage: python probe_p_poset_ceiling.py [max_m]
"""
import gc
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from trd_bivariate import canon_and_aut, contribute, ideals  # noqa: E402

# A000112 (unlabelled posets). Verified against OEIS 2026-08-02 for the terms used.
A000112 = {1: 1, 2: 2, 3: 5, 4: 16, 5: 63, 6: 318, 7: 2045, 8: 16999,
           9: 183231, 10: 2567284}


def rss_bytes():
    """Peak-ish resident size without adding a dependency."""
    try:
        import psutil
        return psutil.Process().memory_info().rss
    except ImportError:
        pass
    if sys.platform == 'win32':
        import ctypes
        from ctypes import wintypes

        class PMC(ctypes.Structure):
            _fields_ = [('cb', wintypes.DWORD),
                        ('PageFaultCount', wintypes.DWORD),
                        ('PeakWorkingSetSize', ctypes.c_size_t),
                        ('WorkingSetSize', ctypes.c_size_t),
                        ('QuotaPeakPagedPoolUsage', ctypes.c_size_t),
                        ('QuotaPagedPoolUsage', ctypes.c_size_t),
                        ('QuotaPeakNonPagedPoolUsage', ctypes.c_size_t),
                        ('QuotaNonPagedPoolUsage', ctypes.c_size_t),
                        ('PagefileUsage', ctypes.c_size_t),
                        ('PeakPagefileUsage', ctypes.c_size_t)]

        c = PMC()
        c.cb = ctypes.sizeof(c)
        ctypes.windll.psapi.GetProcessMemoryInfo(
            ctypes.windll.kernel32.GetCurrentProcess(), ctypes.byref(c), c.cb)
        return c.WorkingSetSize
    import resource
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024


def gen_level(prev_succ, m):
    """One level of gen_posets, kept separable so it can be timed on its own."""
    seen = {}
    for succ in prev_succ:
        for d in ideals(m - 1, succ):
            new = list(succ) + [0]
            for v in range(m - 1):
                if d >> v & 1:
                    new[v] |= 1 << (m - 1)
            new = tuple(new)
            c, aut = canon_and_aut(m, new)
            if c not in seen:
                seen[c] = (new, aut)
    return [v[0] for v in seen.values()], [v[1] for v in seen.values()]


def main():
    max_m = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    print('probe P — poset catalogue vs per-poset work, to m = %d' % max_m)
    print('python %s on %d cores\n' % (sys.version.split()[0], os.cpu_count()))

    base = rss_bytes()
    levels, auts = {0: [()]}, {0: [[()]]}
    gen_rows = []
    for m in range(1, max_m + 1):
        gc.collect()
        t0, r0 = time.perf_counter(), rss_bytes()
        levels[m], auts[m] = gen_level(levels[m - 1], m)
        dt, dr = time.perf_counter() - t0, rss_bytes() - r0
        want = A000112.get(m)
        assert want is None or len(levels[m]) == want, \
            'poset count at m=%d: got %d want %d' % (m, len(levels[m]), want)
        gen_rows.append((m, len(levels[m]), dt, dr, rss_bytes() - base))
        print('  gen m=%-3d posets %-10d %8.2fs  level +%7.1f MB  total %7.1f MB  '
              '%8.1f us/poset' % (m, len(levels[m]), dt, dr / 1e6,
                                  (rss_bytes() - base) / 1e6, dt / len(levels[m]) * 1e6),
              flush=True)

    print('\ncontribute phase — cost of the weighted Burnside per poset')
    con_rows = []
    for n in (max_m - 1, max_m):
        if n < 1:
            continue
        t0 = time.perf_counter()
        npos = 0
        for m in range(1, n + 1):
            for succ, aut in zip(levels[m], auts[m]):
                ra, rn = {}, {}
                contribute(m, succ, aut, n, ra, rn)
                npos += 1
        dt = time.perf_counter() - t0
        con_rows.append((n, npos, dt))
        print('  n=%-3d over %-10d posets  %8.2fs  %8.1f us/poset'
              % (n, npos, dt, dt / npos * 1e6), flush=True)

    print('\nextrapolation to the next level')
    m = max_m
    gp = gen_rows[-1]
    ratio = A000112[m + 1] / A000112[m] if m + 1 in A000112 else None
    if ratio:
        # per-poset generation cost grows because canon_and_aut is superlinear in m;
        # take the measured growth of us/poset across the last two levels as the factor.
        prev, cur = gen_rows[-2], gen_rows[-1]
        per_prev = prev[2] / prev[1]
        per_cur = cur[2] / cur[1]
        growth = per_cur / per_prev
        proj_t = A000112[m + 1] * per_cur * growth
        proj_mem = gp[4] / gp[1] * A000112[m + 1] + gp[4]
        print('  posets on %d points: %d  (%.1fx more than m=%d)'
              % (m + 1, A000112[m + 1], ratio, m))
        print('  per-poset gen cost grew %.2fx from m=%d to m=%d' % (growth, m - 1, m))
        print('  projected gen time  m=%d: %8.0f s  = %6.2f h  = %5.2f core-days'
              % (m + 1, proj_t, proj_t / 3600, proj_t / 86400))
        print('  projected peak RSS  m=%d: %8.1f GB   (holding the level in memory)'
              % (m + 1, proj_mem / 1e9))
        if con_rows:
            n, npos, dt = con_rows[-1]
            per = dt / npos
            proj_c = (sum(A000112[k] for k in range(1, m + 2) if k in A000112)) * per * growth
            print('  projected contribute m=%d: %8.0f s  = %6.2f h  = %5.2f core-days'
                  % (m + 1, proj_c, proj_c / 3600, proj_c / 86400))
            print('\n  WALL = %s' % ('poset generation' if proj_t > proj_c
                                     else 'per-poset weighted Burnside'))
    return 0


if __name__ == '__main__':
    sys.exit(main())
