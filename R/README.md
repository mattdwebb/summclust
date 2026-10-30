# Base-R summclust

Source `summclust.R`; no installation or CRAN package is needed. This independent
implementation follows Matthew D. Webb's local Stata `summclust` 5.1.1 with the
30 September 2026 presentation revision. Frozen pre-mixed test-source hashes
are recorded in `tests/baseline/README.md`. Mathematical attribution: MacKinnon,
Nielsen, and Webb, *Leverage, Influence, and the Jackknife in Clustered Regression
Models: Reliable Inference Using summclust* (2023).

```r
source("R/summclust.R")
i <- 1:120
d <- data.frame(y = sin(i/2), x = sin(i/3), z = cos(i/4),
                group = rep(1:8, each = 15))
fit <- summclust(d, outcome = "y", focal = "x", controls = "z",
                cluster = "group", jackknife = TRUE, rankcheck = TRUE)
print(fit)
fit$cv1                         # full covariance matrix
fit$inference$cv3               # full precision for every explicit coefficient
flatten_summclust(fit, "results.csv")
```

The example above uses deterministic synthetic data.

## Arguments

`summclust(data, outcome, focal, controls, cluster, fevar=NULL, absorb=NULL,
sample=NULL, gstar=FALSE, rho=NULL, jackknife=FALSE, rankcheck=FALSE,
table=FALSE, addmeans=FALSE, regtable=FALSE, mixed=FALSE)`

- `data`: data frame or named list of equal-length vectors. Regression columns
  must be numeric. `outcome` and `focal` name single columns; `controls` names at
  least one further regressor, matching the Stata syntax.
- `cluster`: one or two column names. Numeric and string identifiers are
  supported, including nonconsecutive numeric labels. String sorting uses radix
  order rather than the session's language-dependent collation.
- `fevar`: explicit categorical effects. In one-way mode these are numeric,
  matching the Stata one-way path. Two-way mode also accepts string effects.
  With no absorption, all levels of the first effect enter and the final sorted
  level of every subsequent effect is omitted. With absorption, the final level
  is omitted for every explicit effect. There is no additional constant when
  either effect option is used. Otherwise the column order is focal, constant,
  controls. FE columns follow the controls.
- `absorb`: one effect to remove by within-group demeaning. Each absorbed group
  must lie within one cluster in every primary clustering dimension.
- `sample`: optional logical vector in original row order. Missing sample values
  are excluded, as are missing selected variables and empty string identifiers.
  Infinite numeric values are excluded because Stata has no finite-data infinity.
- `gstar` requests effective-cluster measures. `rho` in `[0.00001,1]` enables them
  too. Cluster/subcluster fixed effects make G*(1) and G*(rho) unavailable. The
  source internally uses rho=.1234567 when omitted, and returns that calculation
  when `gstar=TRUE` and it is defined; the port retains this source behavior.
- `jackknife`: additionally return one-way centered CV3J. Two-way CV3J raises an
  error, as in the authoritative package.
- `mixed`: a single nonmissing logical value, off by default. With two clustering
  dimensions, `TRUE` selects CV31, the mixed estimator `CV3_G + CV3_H - CV1_I`.
  It reuses the intersection CV1 correction `I/(I-1)*(N-1)/(N-rank-absorbed_levels)`.
  One-way use and two-way `jackknife=TRUE` are rejected. Selection never changes
  automatically with intersection count, sample size, or runtime. CV31 is a
  different estimator that may save computation; no superiority is claimed.
- `rankcheck`: calculate actual deleted-design rank losses. In one-way mode it
  also exposes the separate legacy near-zero-coefficient flags. It does not
  change the historical one-way dropped-subsample policy.
- `table`, `addmeans`, `regtable`: expose cluster tables, alternative means, and
  all-coefficient inference. Printing/graphs intentionally do not imitate Stata.
  `addmeans` has no additional effect in two-way mode.

## Results and precision

The returned `summclust` list includes Stata numerical names: `cv1`, `cv3`, optional
`cv3j`, `ng`, `lever`, `partlev`, `betanog`, optional `legacyflag`/`rankfail`, and
the documented focal inference and G* fields. Matrix dimension names give the
coefficient order. Diagnostic vector names give actual cluster identifiers.
`b` is the full coefficient row vector in both modes. `sample_mask` and
`sample_order` identify the original input rows used and their internal order.
`components` contains per-dimension matrices, complete deletion estimates and
diagnostics; `cluster_pairs` gives the actual two-way intersection identifiers.

Two-way `cv1` and `cv3` are conventional inclusion-exclusion matrices and can be
indefinite. Their component matrices and documented aliases are returned.
`cv1_max_var`/`cv3_max_var` and their inference describe only the focal
restriction. `inference` contains coefficient-specific max standard errors for
two-way clustering. There is no full max covariance matrix. All two-way
inference uses min(G,H)-1 degrees of freedom.

With `mixed=TRUE`, `cv31`/`cv31_conventional` contain the raw conventional mixed
covariance. `cv31_conventional_var` and `cv31_conventional_se` concern the focal
coefficient; negative raw variances are retained and their conventional SEs
are `NA`. `inference_conventional$cv31` contains conventional inference for
all coefficients. `cv31_max_var`/`cv31_max_se`, `cv31se`, `cv31t`, `cv31p`,
`cv31lci`, `cv31uci`, and `inference$cv31` use coefficient-specific maxima of
the mixed, G-CV3 and H-CV3 variances. No full max covariance is constructed.
`selected_estimator` is `"CV31"`; labels are `"CV31 (mixed)"` and
`"CV31 (mixed, max)"`. The degrees of freedom remain `min(G,H)-1`.

Ordinary CV3 is not calculated under mixed: its combined matrices and the
intersection CV3 matrix contain `NA`, and no ordinary combined CV3 inference
is returned. Intersection `omitted`, `betanog`, and `legacyflag` are `NA`.
Intersection sizes, leverage, partial leverage, and G-star remain available.
An explicit `rankcheck=TRUE` computes deleted-design ranks separately, without
deleted coefficient solves; this additional matrix work can reduce savings.
`table=TRUE` retains the existing primary-dimension tables. No intersection
omission table or graph is fabricated. `intersection_delete_one_solves` counts
the actual coefficient solves and is zero on both mixed paths. Each component
also reports `delete_one_solves`. The R port has no eigenvalue adjustment option.

```r
ordinary <- summclust(d, "y", "x", "z", cluster=c("firm", "year"))
mixed_fit <- summclust(d, "y", "x", "z", cluster=c("firm", "year"), mixed=TRUE)
mixed_fit$cv31_conventional
mixed_fit$inference$cv31
```

For identical partitions, CV31 is `2*CV3_G - CV1_G`, generally different from
one-way CV3. For a nested partition, it is `CV3_G + CV3_H - CV1_I` with the
actual (finer) intersection partition.

The one-way source rounds its legacy focal macros through Mata `%9.0g` before
computing t statistics, p values, and intervals. The matching top-level fields
preserve that behavior. The double-precision `b`, covariance matrices,
`inference`, `beta_full`, and `cv1se_full`/`cv3se_full` retain unrounded results.
G* fields ending in `_full` likewise retain full precision. The source's global
summary uses float-storage temporary vectors; `clustsum` reproduces that table,
while `clustsum_full` uses the underlying double-precision diagnostics. No float
conversion is applied to `ng`, `lever`, `partlev`, or `betanog`.

One-way singular deletions preserve the historical flag `any(abs(beta_minus)
<=1e-8)`, even when a zero coefficient is legitimate. `inference$cv3drop` and
`inference$cv3jdrop` keep distinct meanings. The legacy top-level `cv3sedrop`,
`cv3tdrop`, etc. reflect CV3J-drop when `jackknife=TRUE`, matching the source's
overwritten macros. Full drop matrices are exposed with `regtable=TRUE`.

`flatten_summclust` produces long CSV rows with `statistic,row,column,value,
row_label,column_label`. List metadata and nested results are excluded; all
top-level numerical matrices/vectors/scalars are included. Values missing in
Stata are represented by R `NA` and blank CSV cells.

## Numerical conventions and validation

Partial leverage residualizes the focal variable on **every other design
column**, including explicit or absorbed effects. CV1 uses the design rank plus
the number of absorbed effects in its residual degrees-of-freedom correction.
CV3 refits deterministic omitted-cluster crossproducts; it does not resample.

The generalized inverse selects independent original columns, setting omitted
rows and columns to zero, rather than using a Moore-Penrose inverse. Its
largest-current-diagonal pivots, relative singular cutoff, Schur update order
and upper-triangle symmetry were checked against Mata `invsym()` on duplicate,
scaled, zero-column and dependent-FE fixtures. Crossproducts accurately sum
rounded binary64 products, using an exact binary-lattice test to retain ordinary
crossproducts where safe. Other sums use high/low binary splitting or an
expansion of rounding residuals. This avoids a mandatory compiled extension.

Deleted-model solves reproduce the current Mata path's right-looking Cholesky
updates, reciprocal column scaling, and two-accumulator dot products in forward
and backward substitution. These arithmetic conventions matter when an
algebraically singular deleted FE design is accepted by the source's numerical
Cholesky criterion. Seventy-five exact-input empirical deletion systems were
checked in addition to the small fixtures. Extremely near-threshold or tied
systems can still be sensitive to platform arithmetic; the validation documents
the tested scope. Rank diagnostics use the separate Mata SVD rule (minimum
matrix dimension times machine epsilon times the largest singular value).
Student-t probabilities and quantiles use base R `stats`.

`clustsum_full` removes the legacy table's float-storage rounding while
retaining its original `G-1` denominator, including a selected deletion
subset. A nonzero singleton subset therefore has zero dispersion under this
diagnostic convention; a zero or nonfinite mean remains missing. This supplemental
field preserves the current baseline's singleton behavior.

Run the independent mathematical tests with:

```text
Rscript R/test_summclust.R
Rscript R/test_mixed.R
Rscript R/run_mixed_fixtures.R
```

`Rscript --vanilla R/certify.R` also runs exact-default parity checks against
the frozen pre-mixed fixture. See `validation/README.md` for the complete
Stata/R/Python comparison, including the required Stata and Python exports.

These checks use QR regression refits, a direct sandwich, exact identities,
row permutations, FE/absorb equivalence, missing values, and error cases.
The shared validation also compares full matrices, diagnostics and inference
with current Stata and the independent Python port.

The supplied fixtures are synthetic. No bootstrap calculations or simulations
are run by these scripts.
