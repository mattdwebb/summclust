# Changelog

## Unreleased string fixed effects — 2026-11-04

- Accept string `fevar` and `absorb` identifiers directly in one-way Stata, and remove R's one-way restriction on string `fevar` columns. Both options now accept strings across all three languages and clustering modes, without a new option; R also accepts factors.
- Use sorted category ordinals for string fixed-effect coefficient names across languages; preserve existing numeric-input names and results. Long string-effect names use a shortened prefix and effect-position suffix to fit Stata's name limit.
- Preserve literal string identifiers such as `NaN` and `Inf` in Python; actual missing values and empty identifiers remain excluded.
- Add equivalent-grouping checks, numeric baseline checks, and shared cross-language comparisons.

## Unreleased Stata help and output revision — 2026-11-04

- Close example paragraphs in the help file and clarify the two-way and mixed-estimator instructions.
- Rework both Stata help files around usage, interpretation, and stored results; remove project history and implementation notes.
- State which component supplies each focal max standard error immediately below the regression table, including the cluster variable name and all components tied within the existing tolerance.
- Preserve all estimator calculations, stored results, and version numbers.

## Unreleased mixed estimator — 2026-10-30

- Added explicit, default-off `mixed` for two-way CV31 = CV3_G + CV3_H - CV1_I, with separate conventional and max-adjusted inference returns.
- Skip intersection delete-one coefficient solves and mark their uncomputed diagnostics missing; preserve ordinary defaults and version numbers.
- Include the base-R and NumPy implementations and portable mixed certification with synthetic cross-language fixtures.

## Unreleased presentation revision — 2026-09-30

- Restored the standalone twowayjack inference and diagnostic table layouts, headings, alignment, separators, cluster names and mixed numeric formats.
- Retained correct min(G,H)-1 confidence intervals and independently verified focal partial leverage; no estimator or returned-result values changed.
- Wrapped notes, warnings and input-error messages at word boundaries; preserved complete long cluster names at ordinary line widths.
- Restored column headings in optional cluster and coefficient tables and made max selectors display the actual cluster names.
- Added separate numerical and text/SMCL formatting checks, including long names, options, negative variances, ties and few-cluster warnings.
- Independently checked the recovered nlswork specification, the available minimum-wage earnings extract, and seven tsetse outcomes. See `docs/AUDIT_OUTPUT_2026-09-30.md` for evidence, qualifications and commands.

## 5.1.1 development release — 2026-09-23

- Integrated the native 5.1.0 two-way candidate into the main Stata package layout, including the thin `twowayjack` delegate and local installation metadata.
- Corrected native fixed-effect coefficient labels to match the actual omitted levels.
- Isolated native Mata work arrays with temporary names, preserving caller objects on success and error.
- Reject native two-way inference when the focal coefficient is not numerically identified in the full model or any cluster deletion, or a required variance component is missing. Nuisance-only rank loss remains supported.
- Added explicit relative-tolerance masks for tied max-variance components.
- Preserve double precision when copying the one-way dependent variable.
- Count explicit rank after absorption for one-way CV1 degrees of freedom and reject nonpositive residual degrees of freedom.
- Added independent fixed-effect refit and regression tests, an isolated local-install test, and a portable-path Windows test runner. Six suites passed on StataNow/SE 18.5.
- Updated the near-collinearity certification to require refusal when the focal coefficient cannot be numerically resolved; retained a resolved collinear case.
- Retained legacy one-way warning/drop semantics and public results for compatibility.

## 5.1.0 development candidate — 2026-08-23

- Added native `cluster(g h)` syntax while preserving the one-way path.
- Builds the two-way estimation sample and design once and reuses one internal
  component routine for G, H, and nonempty intersections.
- Returns separate CV1/CV3 component matrices and conventional
  `V_G + V_H - V_I` covariance matrices.
- Added focal max-se CV1/CV3 inference using `min(G,H)-1` degrees of freedom,
  without inventing a full max covariance matrix.
- Added dimension-specific cluster sizes, leverage, partial leverage,
  delete-cluster estimates, rank checks, coefficients of variation, and G*(0).
- Added generalized-inverse support for explicit fixed effects and a strict
  two-way nesting rule for `absorb()`.
- Rejects two-way centered CV3J because that alternative construction is not
  studied in the current two-way jackknife paper.
- Suppresses two-way graphs by policy and returns graph-ready diagnostic
  matrices instead.
- Added direct algebra, edge-case, fixed-effect, numerical-stability, negative
  variance, one-way-regression, legacy-wrapper, and performance tests.

## 5.0.0 — 2026-08-23

- Replaced per-cluster Stata/Mata crossings and dynamically growing matrices
  with a preallocated Mata core.
- Corrected partial leverage to use the Frisch-Waugh-Lovell residual of the
  focal regressor on the rest of the estimation design. This intentionally
  changes partial-leverage results relative to 4.22.
- Corrected absorbed-group processing so results no longer depend on incoming
  row order.
- Restored the caller's incoming observation order before returning.
- Added full `r(cv1)`, `r(cv3)`, and, with `jackknife`, `r(cv3j)` covariance
  matrices with coefficient names.
- Added opt-in `rankcheck`, returning the compatibility flag in
  `r(legacyflag)` and actual deleted-design rank loss in `r(rankfail)`.
- Preserved the historical warning/drop policy by default. This policy is
  coefficient-threshold based and is not itself a rank test.
- Corrected negative off-diagonal entries in the full CV3J covariance matrix.
- Added coefficient names to the full covariance matrices and corrected CV3J
  table titles.
- Corrected dropped-subsample CV3J p-values and confidence-interval degrees of
  freedom.
- Corrected `rho()` option detection for valid values near 0.12345.
- Corrected G* fixed-effect checks for string cluster identifiers and
  estimation-sample restrictions.
- Removed fixed private matrix names and fixed graph-component names introduced
  during optimization, avoiding collisions with user objects.
- Added clear errors for empty samples and samples containing fewer than two
  clusters.

## 4.22 — 2024-07-15

- Corrected a `regtable` bug.
