# research-log

Dated scripts and their recorded outputs, exactly as they were run, for the results in two papers
about arrow-type semigroupoids (transitive relations with no isolated point) that are being prepared:

- **component decomposition, defect grading, and inversion of the published table** (the connected
  counts `D(e,v)`, `C_d`, `D_N`, the polynomials `E_N`, the three cells `T(11,6)`, `T(11,7)`, `T(11,8)`);
- **enumeration to ten points, and labelled rows 8-11 of OEIS A399095**.

Nothing in this folder was edited after it was run. Files keep their original dated names, so the order
in which things were done can be read off the names. A script that turned out to be wrong or slower than
hoped is not deleted; the outputs that show it are kept (`2026-08-26-stream-n11.txt` records a walk that
took 22,947.5 s).

## Layout

`runs/` holds everything. Scripts find their inputs next to themselves and in `../sources/`, so run them
from inside `runs/`. `trd_stream.py` and `trd_bivariate.py` are the copies that the dated scripts import;
the bivariate one differs from `../src/trd_bivariate.py` (it exposes `gen_posets` for the labelled
scripts), so both are kept as they were run.

| step | script | output | what it establishes |
|---|---|---|---|
| connected classes, `e <= 6` | `2026-09-09-defect-grading.py` (uses `2026-08-30-low-arrow-cells.py`) | `2026-09-09-defect-grading-run2.txt`, `2026-09-09-connected-classes.json` | 20 types, 293 classes; `C_0..C_5`, `D_0..D_5` |
| corner of the table | `2026-09-09-corner-table-v2.py` | `2026-09-09-corner-table-v2.txt` | the three cells of row 11, masked comparison |
| `E_N` for `N <= 5` | `2026-09-09-defect-polynomials.py` | `2026-09-09-defect-polynomials.txt` | the polynomials by enumeration |
| inversion of the b-file | `2026-09-11-bfile-inversion.py` | `.txt`, `-D.json` | `D(e,v)` for `v <= 10`; `C_0..C_8`, `D_0..D_8` |
| connected triangle | `2026-09-11-connected-triangle.py` (+ `2026-09-11-oeis-prog-triangle.py`) | `.txt`, `-bfile.txt` | OEIS A399913 and its checks |
| eventual linearity | `2026-10-13-eventual-linearity.py` | `.txt` | `A_i`, `B_i` and 28 cell checks |
| labelled rows 8-11 | `2026-08-26-labeled-triangle.py`, `2026-08-26-labeled-stream.py` | `2026-08-26-labeled-triangle-n10.txt`, `2026-08-26-stream-n11.txt`, `2026-08-26-row11-wrapped.txt` | rows 8-11 of A399095, 148/148 published terms reproduced |

`sources/2026-09-11-b399030.txt` is the published b-file of OEIS A399030 (rows 1-10) as fetched on
11 September 2026; the inversion reads it.

## Two honest notes

- `2026-10-13-eventual-linearity.py` was written **after** the pattern it checks had been seen in the
  tables. Its header says so. It tests an argument, not a prediction.
- Docstrings refer to planning pages in a private notes vault (`wiki/concepts/...`). Those pages are not
  part of this archive; the scripts do not need them.

Code: MIT. Data and outputs: CC0, as for the rest of the repository.