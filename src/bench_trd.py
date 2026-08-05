"""Benchmark: level-storing (trd_bivariate) vs streaming (trd_stream), same n.

Both implementations do identical total work — every poset on 1..n points folded into
the arrow tallies for every n' in [m, n] — so end-to-end wall time and peak working set
are directly comparable. Run as two separate processes, sequentially, so neither
distorts the other's memory or cache.

The totals are asserted equal between the two before any timing is reported: a faster
implementation that computes something else is not a faster implementation.

usage: python bench_trd.py orig|stream <n>          # one side, prints one JSON line
       python bench_trd.py both <n>                 # runs both sequentially, compares
"""
import ctypes
import json
import subprocess
import sys
import time
from ctypes import wintypes


class PMC(ctypes.Structure):
    _fields_ = [('cb', wintypes.DWORD), ('PageFaultCount', wintypes.DWORD),
                ('PeakWorkingSetSize', ctypes.c_size_t), ('WorkingSetSize', ctypes.c_size_t),
                ('QuotaPeakPagedPoolUsage', ctypes.c_size_t), ('QuotaPagedPoolUsage', ctypes.c_size_t),
                ('QuotaPeakNonPagedPoolUsage', ctypes.c_size_t), ('QuotaNonPagedPoolUsage', ctypes.c_size_t),
                ('PagefileUsage', ctypes.c_size_t), ('PeakPagefileUsage', ctypes.c_size_t)]


def peak_bytes():
    """Peak working set. psutil first — the ctypes fallback below returned 0 in the
    first benchmark run and the failure was silent, so it is asserted here."""
    try:
        import psutil
        return psutil.Process().memory_info().peak_wset
    except ImportError:
        pass
    c = PMC()
    c.cb = ctypes.sizeof(c)
    ok = ctypes.windll.psapi.GetProcessMemoryInfo(
        ctypes.windll.kernel32.GetCurrentProcess(), ctypes.byref(c), c.cb)
    assert ok and c.PeakWorkingSetSize, 'GetProcessMemoryInfo failed (returned %r)' % ok
    return c.PeakWorkingSetSize


def one(mode, n):
    t0 = time.perf_counter()
    if mode == 'orig':
        from trd_bivariate import enumerate_trds
        all_t, noiso_t = enumerate_trds(n)
    else:
        from trd_stream import run
        all_t, noiso_t, _ = run(n)
    dt = time.perf_counter() - t0
    return {'mode': mode, 'n': n, 'seconds': dt, 'peak_bytes': peak_bytes(),
            'totals': {str(k): sum(v.values()) for k, v in sorted(all_t.items())},
            'noiso': {str(k): sum(v.values()) for k, v in sorted(noiso_t.items())}}


def main():
    mode, n = sys.argv[1], int(sys.argv[2])
    if mode in ('orig', 'stream'):
        print('@@RESULT@@' + json.dumps(one(mode, n)))
        return 0

    res = {}
    for m in ('orig', 'stream'):
        print('=== running %s at n=%d' % (m, n), flush=True)
        p = subprocess.run([sys.executable, __file__, m, str(n)],
                           capture_output=True, text=True)
        line = [l for l in p.stdout.splitlines() if l.startswith('@@RESULT@@')]
        if not line:
            print(p.stdout[-2000:]); print(p.stderr[-2000:])
            raise SystemExit('%s produced no result' % m)
        res[m] = json.loads(line[0][len('@@RESULT@@'):])

    a, b = res['orig'], res['stream']
    assert a['totals'] == b['totals'], 'TOTALS DIFFER\n%s\n%s' % (a['totals'], b['totals'])
    assert a['noiso'] == b['noiso'], 'NOISO DIFFER'
    print('\ncorrectness: totals and no-isolated counts identical for n=1..%d' % n)
    print('  n=%-3d all=%s' % (n, a['totals'][str(n)]))

    print('\n%-10s %12s %14s' % ('', 'wall (s)', 'peak RSS (MB)'))
    for m in ('orig', 'stream'):
        print('%-10s %12.1f %14.1f' % (m, res[m]['seconds'], res[m]['peak_bytes'] / 1e6))
    print('\nspeed  : %.2fx   (>1 means streaming is faster)'
          % (a['seconds'] / b['seconds']))
    print('memory : %.2fx   (>1 means streaming uses less)'
          % (a['peak_bytes'] / b['peak_bytes']))
    return 0


if __name__ == '__main__':
    sys.exit(main())
