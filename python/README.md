# Python summclust

## Explicit mixed estimator (CV31)

```python
ordinary = summclust(data, "y", "x", ["z"], ["firm", "year"])
mixed = summclust(data, "y", "x", ["z"], ["firm", "year"], mixed=True)
```

`mixed=False` is the default. The method never changes automatically with
intersection count. CV31 combines `CV3_G + CV3_H - CV1_I`, reusing the same
finite-sample-adjusted intersection CV1 as ordinary two-way CV1. This is a
different estimator that may reduce computation; no superiority is claimed.
One-way `mixed=True` and centered two-way `jackknife=True` are rejected.
`mixed` accepts Python or NumPy booleans only.

`cv31`/`cv31_conventional` retain the raw full covariance. The focal
`cv31_conventional_var` and `cv31_conventional_se` are separate from
`cv31_max_var`, `cv31_max_se`, and `inference['cv31']`: the latter uses
max(CV31 focal variance, CV3_G focal variance, CV3_H focal variance), with
`min(G,H)-1` degrees of freedom. Full-precision `cv31_max_t/p/lci/uci`
and existing-style `cv31t/p/lci/uci` aliases are available.
`inference_conventional['cv31']` supplies conventional inference; a negative
variance produces NaN SE/t/p/CI. There is no elementwise-max covariance matrix.
`selected_estimator`, `estimator_label`, `inference_label`, and
`regsum_columns` identify `CV31 (mixed)` and `CV31 (mixed, max)`.

Under mixed, `cv3`/`cv3_conventional`, `cv3_intersection`, intersection
`betanog`/`betas_deleted`/`legacyflag`, and its omit-one summary are NaN.
The G/H CV3 components remain computed. Intersection sizes, ordinary and
partial leverage, and available G-star quantities remain computed.
`rankcheck=True` computes requested intersection rank diagnostics directly,
without coefficient solves; it can still add work. `table=True` retains
the existing G/H tables. There is no explicit intersection omit-one request
in the current interface. `intersection_delete_one_solves` reports zero;
component solve counters are incremented at the actual solver call.
Pure CV3 and its deletion diagnostics require an ordinary call.
No eigenvalue adjustment option exists in this API.

Independent implementation of Matthew Webb's current Stata `summclust`, with
the cluster diagnostics and methods of MacKinnon, Nielsen, and Webb. The exact
Stata baseline is 5.1.1 with the September 30, 2026 presentation revision. Runtime
requirements: Python 3.10 or newer and NumPy 1.22 or newer. No pandas, SciPy,
statsmodels, plotting package, Stata installation, or R installation is needed
to use this module. See the repository LICENSE and source attribution.

This module contains no bootstrap procedure. The supplied tests are deterministic
and use synthetic data only.

```python
from summclust import read_csv, summclust

data = read_csv("data.csv", string_columns=["firm"])
result = summclust(
    data, outcome="y", focal="x", controls=["z", "w"], cluster="firm",
    fevar=["year"], rho=0.5, jackknife=True, rankcheck=True, table=True,
)
print(result["coefficient_names"], result["b"])
print(result["inference"]["cv3"])  # full-precision beta, SE, t, p, 95% CI, df
```

`data` is a mapping from column names to equal-length one-dimensional sequences.
Outcome, focal regressor, and controls must be numeric. Cluster identifiers may
be numeric or strings. `controls` must contain at least one name, matching the
Stata minimum-varlist requirement. No formula parser or automatic data download
is involved.

`sample` accepts a Boolean mask over the original observations. Complete-case
selection includes every regression column, cluster identifier, fixed effect,
and absorbed identifier. Numeric nonfinite values, `None`, and empty string
identifiers are missing. `read_csv` treats empty numeric cells and `.` as NaN;
declare string columns explicitly. Original arrays and their order are preserved.

`cluster=["firm", "year"]` requests native two-way inference. `fevar` accepts one
name or a sequence; `absorb` accepts one name. Both accept numeric or string
identifiers in one-way, two-way, and mixed estimation; no manual encoding is
needed. Absorbed groups must be nested inside each clustering dimension.
`gstar`, `rho`, `jackknife`, `rankcheck`,
`table`, `addmeans`, and `regtable` implement the corresponding numeric options.
`rho` must be between 0.00001 and 1, inclusive. Omitting it uses the source's
internal 0.1234567 when `gstar=True`. Graphs and Stata-specific printing are not
implemented; all graph inputs are returned.

String fixed-effect levels are sorted after sample restrictions and missing
values have been removed. Their coefficient names use the original column name
and a 1-based level position, such as `industry_1`, so spaces and punctuation in
category labels do not enter coefficient names. If this would exceed 32
characters, the column-name prefix is shortened and the suffix is
`_f<effect position>_<level position>`, keeping different long effect names
distinct. Numeric one-way fixed-effect
names retain their actual values; two-way names use level positions for both
types. Empty strings, `None`, and numeric nonfinite identifiers are missing;
nonempty strings such as `"NaN"` or `"Inf"` are ordinary category labels.

For a CSV with string columns `firm` and `industry`, declare both when reading
the file. The columns `y`, `x`, `z`, and `year` are numeric in this example:

```python
data = read_csv("data.csv", string_columns=["firm", "industry"])
oneway = summclust(data, "y", "x", ["z"], "firm", fevar="industry")
mixed = summclust(data, "y", "x", ["z"], ["firm", "year"],
                  fevar="industry", mixed=True)
absorbed = summclust(data, "y", "x", ["z"], "firm", absorb="firm")
```

If your arrays already contain strings, pass them directly; no numeric encoding
is needed. Omit `mixed=True` for ordinary two-way estimation. The last call
uses one-way clustering so each absorbed firm lies within a cluster. When also
clustering by year, firms observed in several years need `fevar="firm"`
instead of `absorb="firm"`.

The returned dictionary includes:

- `b`, `coefficient_names`, `N`, `k`, `rank`, `coeffnum`, `df_r`, `sample_mask`;
- one-way `cv1`, `cv3`, optionally `cv3j`; `ng`, `lever`, `partlev`, `betanog`;
- optional `legacyflag` and `rankfail`; flagged-cluster dropped-subsample
  covariance/inference with the same availability rules as Stata;
- optional G-star, cluster `scall` table, `clustsum` summary, `bonus` additional
  means, and `regresstab` numeric regression table;
- two-way `cv1`/`cv3` conventional inclusion-exclusion matrices, each G/H/I
  component and documented aliases, cluster-specific vectors, `dimsummary`,
  and optional `table_g`/`table_h`;
- two-way focal max variances/SEs, selected component names, and tie masks
  (1 = two-way, 2 = G, 4 = H); `regsum` contains coefficient-specific max SEs.

`components` contains full-precision per-dimension numerical arrays and actual
`cluster_ids`, including identifier pairs for the intersection. Align vectors by
these identifiers when comparing software, because string ordering may differ
across locales. `export_rows(result, specification)` produces the comparison
indexed numeric CSV rows, including coefficient labels. `flatten_summclust`
provides a compact alternative with identifiers as row labels.

There is no full max covariance matrix: two-way inference uses coefficient-wise
maxima of the conventional, G, and H variances, with `min(G,H)-1` degrees of
freedom. Negative conventional variances are preserved. Two-way CV3J is rejected.
The focal coefficient must be identified in the full and each computed
deleted-cluster two-way model (G/H only under mixed). Rank deficiency confined
to nuisance coefficients is allowed.

The design is `[focal, constant, controls]` when there are no fixed effects.
Otherwise it is `[focal, controls, indicators]`: all levels of the first explicit
fixed effect are included without absorption, and the final sorted level is
omitted for subsequent effects and every effect with absorption. Absorption
demeans the design and outcome; CV1's parameter count includes the transformed
design rank plus the number of absorbed groups. Partial leverage residualizes
the focal regressor on **all** other design columns.

The generalized inverse uses an adaptive diagonal-pivot symmetric sweep,
zeroing omitted coordinates rather than distributing coefficients across a
Moore-Penrose solution. Its relative pivot cutoff and arithmetic ordering were
checked with direct Mata `invsym()` probes, including dependent fixed effects.
Crossproducts accurately accumulate rounded double-precision products, with a
binary-lattice check retaining fast matrix multiplication for exact sums.
The Cholesky path and triangular solves preserve Mata's update and summation
order; LAPACK's equally valid ordering can select a different representation
when a deleted fixed effect makes nuisance coefficients unidentified.
Rank diagnostics use the Mata singular-value tolerance. Floating-point pivot
choices can still change in numerically marginal designs; `validation/README.md`
states the cases and tolerances actually checked. The one-way legacy flag means
that any omitted-cluster coefficient has magnitude at most 1e-8; it is deliberately
distinct from true rank loss and does not exclude estimates from primary CV3.

For one-way compatibility, `beta`, `cv1se`, `cv3se`, `cv3Jse`, G-star aliases, and
their derived t/p/CI aliases reproduce Stata's legacy `%9.0g` conversion before
inference. Full precision remains in `b`, covariance matrices, `inference`, and
the corresponding `_full` aliases. `clustsum` reproduces the float32 temporary
vectors used in Stata's summary display; `clustsum_full` retains double-precision
summary inputs. Two-way results use full precision throughout. These deliberate
compatibility conventions never round inputs to a covariance calculation.

Student-t tails and inverse tails are available as `student_t_sf(t, df)` and
`student_t_isf(p, df)`. They use an incomplete-beta continued fraction and
bracketed inversion; no normal approximation is substituted. They are checked
against analytic low-df cases, base R, and Stata across the relevant degrees of
freedom and probabilities, including tails as small as 1e-12.

Run the independent mathematical tests from this directory:

```text
python -m unittest -v
```

From the repository root, `python validation/python_certification.py` also
checks exact default parity against the frozen pre-mixed test fixture.
See `validation/README.md` for the full Stata/R/Python comparison.

Invalid specifications raise `SummclustError`, a `ValueError` subclass with a
Stata-style `code` attribute. Numerical output is returned without printing, and
diagnostic messages are collected in `result["warnings"]`.
