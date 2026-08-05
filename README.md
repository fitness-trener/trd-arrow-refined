# Arrow-type semigroupoids / transitive relation digraphs, refined by arrow count

Complete enumeration of **arrow-type semigroupoids by objects and arrows, for 1–10 objects**,
computed from a structure theorem rather than by searching semigroupoids.

Equivalently — see [Identification](#identification) — these are **transitive relation digraphs
(TRDs) on `n` unlabelled points with no isolated point, counted by the number of related pairs.**

This work extends Table 1 of Egri-Nagy & Nehaniv, *Computational Exploration of Finite Semigroupoids*
([arXiv:2509.00837](https://arxiv.org/abs/2509.00837)), which gives the refinement completely for
`n ≤ 7` and partially at `n = 8`. That paper is the source of the refinement and the motivation for
this computation.

## Results

| objects | arrow-type semigroupoids | arrow counts | peak cell |
|---|---|---|---|
| 7 | 227,844 | 4–49 | 15 arrows |
| 8 | 3,898,213 | 4–64 | 383,274 at 19 arrows |
| 9 | 86,180,547 | 5–81 | 7,548,078 at 23 arrows |
| **10** | **2,465,304,381** | **5–100** | **193,622,278 at 27 arrows** |

Full per-cell data: [`data/trd-table-n10.csv`](data/trd-table-n10.csv) (columns
`objects, arrows, arrow_type_semigroupoids, all_trds`).

### What is new here, and what is not

**New:** the refinement by arrow count for `n ≥ 8`. It is not derivable from any published sequence
and is not in OEIS at any `n`.

**Not new:** the object-wise totals. They are exactly the first differences of
[A091073](https://oeis.org/A091073) (transitive relations on `n` unlabelled points), because a
relation splits uniquely into its isolated points and a relation with none, so the totals telescope.
A091073 has been known to 15 points since Brinkmann & McKay (2005). Their agreement with this
computation is a **check on the pipeline, not a result.**

## Verification

Nothing here is worth reading if this section fails. It does not.

| check | result |
|---|---|
| unlabelled poset counts vs [A000112](https://oeis.org/A000112), `m ≤ 10` | 1, 2, 5, 16, 63, 318, 2045, 16999, 183231, 2567284 — all match |
| all TRDs by points vs [A091073](https://oeis.org/A091073), `n ≤ 10` | all match, including `n = 10`: **2,555,630,036** |
| no-isolated totals vs Egri-Nagy & Nehaniv column sums, `n ≤ 7` | 1, 6, 31, 203, 1653, 17156, 227844 — all match |
| their individual Table 1 cells | every legibly rendered cell reproduces |
| **two independent implementations** | all 185 bivariate cells at `n = 8` identical |

Run logs are in [`results/`](results/).

## Identification

Arrow-type semigroupoids are transitive relation digraphs: no multiple edges, at most one loop per
point, and `(x,y), (y,z)` present implies `(x,z)` present. `src/identify_arrowtype.py` brute-forces
the correspondence independently to `n = 4`.

The offset against A091073 is the isolated-point convention: arrow-types that *occur* leave nothing
isolated, which is why the published column sums do not appear in an OEIS search directly.

## Method

From Brinkmann & McKay §2: contracting the strong components of a TRD gives a poset `P`, each
component a directed clique. So a TRD on `n` points is exactly

- a poset `P` on `m` points,
- sizes `n_v ≥ 1` with `Σ n_v = n`,
- a loop / no-loop choice on each **singleton** component,

with two data isomorphic iff related by an automorphism of `P`. The arrow count is closed-form:

```
arrows = Σ_v f(n_v) + Σ_{(v,w) ∈ E(P)} n_v · n_w,     f(k) = k² for k ≥ 2,  f(1) = loop ? 1 : 0
```

It is `Aut(P)`-invariant, so counting orbits per arrow count is weighted Burnside within each weight
class. **The semigroupoids are never enumerated; the posets are.**

### Two implementation facts worth knowing

Both were measured, and both contradict what looked obvious beforehand.

**1. The Burnside sum is the expensive phase, not poset generation.** Per-poset generation cost is
flat-to-falling with `m` (1,089 → 820 µs from `m = 8` to `m = 9`); the per-poset Burnside more than
doubles per level (1,251 → 2,921 µs). Reaching for a faster poset generator (`nauty`, `genposetg`)
does not help — and the poset catalogue is not a limit either, since A000112 is known to 16 points.

**2. The memory wall is the automorphism groups, not the poset count.** `trd_bivariate.py` holds
every level with `Aut(P)` listed explicitly, and an antichain on `m` points has `m!` automorphisms —
3.6M at `m = 10`, 40M at `m = 11`. One poset.

`trd_stream.py` fixes both:

- **streams** the posets by McKay canonical augmentation, depth-first — a child is accepted only if
  the point just added lies in the `Aut`-orbit of the canonically-chosen maximal point, so each
  isomorphism class is reached exactly once with no global table. Memory is `O(depth)` posets.
- **collapses `Aut` to cycle partitions** with multiplicity `∏_B (|B|−1)!`, since the Burnside loop
  reads each automorphism only through its cycles *as sets*. Where `Aut(P) = ∏ S_{C_i}` the
  partitions are enumerated directly and `Aut` is never listed: Bell(10) = 115,975 against 3,628,800.

Measured at `n = 9`, both implementations end-to-end:

| | wall (s) | peak RSS (MB) |
|---|---|---|
| `trd_bivariate.py` (level-storing) | 603.1 | 275.0 |
| `trd_stream.py` (streaming) | **533.3** | **27.4** |

`n = 10` took 4.4 h and 3.25 GB with the level-storing implementation. `n = 11` is out of reach for
it (~64 GB) and looks like a few core-days for the streaming one — an estimate from a cost model
that has already been wrong by a factor of two once, so treat it as indicative.

## Layout

```
src/
  trd_bivariate.py          level-storing implementation (produced the n=10 table)
  trd_stream.py             streaming implementation (canonical augmentation + cycle partitions)
  bench_trd.py              benchmark, correctness-gated before timing
  probe_p_poset_ceiling.py  phase-by-phase cost measurement
  identify_arrowtype.py     independent brute-force check of the identification
  make_oeis_submission.py   generates the OEIS b-file and data line
data/
  trd-table-n10.csv         the full table, all cells
  oeis-b-file.txt           385 terms, linear index over rows
  triangle-rows.txt         rows as k = 1..n^2, for reading
results/
  run-n10.log, run-n9.log, benchmark-n9.txt, cost-model-n9.txt
```

## Reproducing

Python 3.11+, no dependencies (`psutil` optional, for the benchmark's memory figures).

```bash
python src/trd_bivariate.py 9     # ~10 min, writes trd-table-n9.csv
python src/trd_stream.py 9        # ~9 min, same output, ~10x less memory
python src/bench_trd.py both 9    # both, sequentially, with a correctness gate
python src/trd_bivariate.py 10    # ~4.4 h, ~3.3 GB
```

Every script asserts its published check values before reporting anything. A mismatch aborts.

## OEIS

The triangle is being submitted as an irregular triangle read by rows, `T(n,k)` = number of
transitive relations on `n` unlabelled points with no isolated point and exactly `k` related pairs,
`1 ≤ k ≤ n²`, with arXiv:2509.00837 cited as the origin of the refinement. The b-file is
[`data/oeis-b-file.txt`](data/oeis-b-file.txt). *(A-number to be added once assigned.)*

The object-wise totals are **not** being submitted — see [above](#what-is-new-here-and-what-is-not).

## References

- A. Egri-Nagy and C. L. Nehaniv, *Computational Exploration of Finite Semigroupoids*,
  [arXiv:2509.00837](https://arxiv.org/abs/2509.00837), 2025.
- G. Brinkmann and B. D. McKay, *Counting unlabelled topologies and transitive relations*,
  J. Integer Seq. **8** (2005), Article 05.2.1.
- G. Pfeiffer, *Counting transitive relations*, J. Integer Seq. **7** (2004), Article 04.3.2.
- OEIS [A091073](https://oeis.org/A091073), [A000112](https://oeis.org/A000112).

## License

Code: MIT ([`LICENSE`](LICENSE)). Data in `data/` and `results/`: CC0 — enumeration results are
facts, and they should be reusable without conditions.

## Citation

See [`CITATION.cff`](CITATION.cff). Once archived on Zenodo, cite the DOI.
