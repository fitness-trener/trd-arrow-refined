"""Generate the OEIS submission payload for the arrow-refined TRD triangle.

Reads the n=10 CSV and emits, into this directory:
  2026-08-05-oeis-triangle-data.txt   the DATA line (truncated at OEIS's ~260 char guide)
  2026-08-05-oeis-triangle-bfile.txt  the b-file, linear index over rows
  2026-08-05-oeis-triangle-rows.txt   full rows, for the EXAMPLE field and for checking

Row structure: row n runs k = 1..n^2 inclusive, with zeros where no object exists, so
row lengths are n^2 and the triangle is `tabf`. Zeros are real: k = 8 is empty at n = 3.

usage: python make_oeis_submission.py
"""
import collections
import csv
import os

SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data',
                   'trd-table-n10.csv')
MAXN = 10
A091073 = [1, 2, 8, 39, 242, 1895, 19051, 246895, 4145108, 90325655, 2555630036]


def load():
    rows = collections.defaultdict(dict)
    with open(SRC, encoding='utf-8') as f:
        for r in csv.DictReader(f):
            n, k = int(r['objects']), int(r['arrows'])
            v = int(r['arrow_type_semigroupoids'])
            if v:
                rows[n][k] = v
    return rows


def main():
    rows = load()

    # invariant checks before anything is written for submission
    totals = [sum(rows[n].values()) for n in range(1, MAXN + 1)]
    diffs = [A091073[n] - A091073[n - 1] for n in range(1, MAXN + 1)]
    assert totals == diffs, 'row totals are not the first differences of A091073'
    for n in range(1, MAXN + 1):
        assert max(rows[n]) == n * n, 'row %d does not end at n^2' % n
        assert min(rows[n]) == -(-n // 2), 'row %d does not start at ceil(n/2)' % n

    flat, bfile, out_rows = [], [], []
    idx = 1
    for n in range(1, MAXN + 1):
        row = [rows[n].get(k, 0) for k in range(1, n * n + 1)]
        out_rows.append((n, row))
        for v in row:
            flat.append(v)
            bfile.append('%d %d' % (idx, v))
            idx += 1

    data = ', '.join(str(v) for v in flat)
    cut = data[:260].rsplit(',', 1)[0]
    with open('2026-08-05-oeis-triangle-data.txt', 'w', encoding='utf-8') as f:
        f.write(cut + '\n\n')
        f.write('# %d terms shown above (rows n=1..%d partial); full data in the b-file\n'
                % (cut.count(',') + 1, MAXN))
    with open('2026-08-05-oeis-triangle-bfile.txt', 'w', encoding='utf-8') as f:
        f.write('\n'.join(bfile) + '\n')
    with open('2026-08-05-oeis-triangle-rows.txt', 'w', encoding='utf-8') as f:
        for n, row in out_rows:
            f.write('n=%d (k=1..%d): %s\n' % (n, n * n, ', '.join(map(str, row))))

    print('flat terms: %d   b-file lines: %d' % (len(flat), len(bfile)))
    print('data line (truncated to %d chars):' % len(cut))
    print(cut)
    print()
    print('row totals  :', totals)
    print('A091073 diff:', diffs, '(identical - asserted)')
    print()
    for n, row in out_rows[:5]:
        print('n=%d: %s' % (n, ', '.join(map(str, row))))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
