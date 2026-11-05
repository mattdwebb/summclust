<!-- Generated from the prose-only README.Rmd. -->

# summclust

Diagnostics and inference for linear regression with one-way or two-way clustering, implemented in Stata, base R, and Python/NumPy.
Version **5.1.1**, with an unreleased presentation revision dated 30 September 2026. Authors: James G. MacKinnon, Morten O. Nielsen, and Matthew D. Webb.

This development source includes the explicit mixed estimator (CV31) in all three languages. Stata has passed nine certification suites plus targeted mixed tests on StataNow/SE 18.5 for Windows. It declares Stata 13 compatibility; older Stata releases have not been independently tested in this integration. R uses only standard distributed packages; Python requires NumPy. This commit does not create a new versioned release.

The language implementations live in [Stata](Stata/), [R](R/README.md), and [python](python/README.md). The repository includes small synthetic verification fixtures, not empirical replication datasets.

## Install from a local copy

Extract or clone the package, then point Stata at its `Stata` subdirectory:

```stata
net install summclust, from("C:/path/to/summclust/Stata") replace
discard
which summclust
help summclust
```

`which summclust` should show version 5.1.1. An ado file in the current working directory or earlier on `adopath` can shadow an installed copy. The package contains `summclust`, its help file, and the compatibility command `twowayjack` with its help file. It uses built-in Stata/Mata functions without other ado dependencies.

## Use

The first variable is the outcome, the second is the focal regressor, and the remaining variables are controls. Supply at least one control. Factor-variable notation is not accepted in this varlist; use explicit indicators or `fevar()`.

Fixed-effect identifiers in `fevar()` and `absorb()` may be numeric or strings in both one-way and two-way models. R also accepts factors. Categories are handled automatically within the estimation sample, so no manual encoding is needed. Absorbed groups must still be nested within every clustering dimension.

```stata
summclust y treatment control, cluster(firm) nograph
summclust y treatment control, cluster(firm year) nograph
summclust y treatment control, cluster(firm year) mixed nograph
summclust y treatment control, cluster(firm year) fevar(firm year) rankcheck nograph
```

Default two-way mode computes the G, H, and nonempty intersection components of CV1 and uncentered CV3. It returns the conventional covariance `V_G + V_H - V_I`, including negative diagonal entries when they occur. Displayed focal inference uses the largest of the conventional, G, and H variances with `min(G,H)-1` degrees of freedom. This scalar rule does not define a joint covariance matrix.

The explicit `mixed` option selects CV31 = CV3_G + CV3_H - CV1_I. It is off by default and never selected automatically. It skips intersection delete-one coefficient solves; uncomputed intersection coefficients and ordinary CV3 results are missing, while inexpensive diagnostics remain available. `r(cv31_conventional)` and `r(cv31_max_se)` distinguish the raw mixed covariance from max-adjusted inference, labelled CV31 (mixed, max). CV31 is a different estimator that may save computation; no superiority is claimed. See `help summclust` for diagnostics and unsupported combinations.

Useful two-way results include `r(b)`, `r(cv1_conventional)`, `r(cv3_conventional)`, `r(cv1_max_se)`, `r(cv3_max_se)`, and the component and diagnostic matrices documented in `help summclust`. Explicit tie masks identify all components attaining the maximum within a relative tolerance of 1e-10.

Use `sample(expression)` for restrictions. `twowayjack` is now a thin delegate to the native command and also accepts `if` and `in`.

R and Python expose the same explicit choice through their existing interfaces:

```r
source("R/summclust.R")
fit <- summclust(data, "y", "x", "z", c("g", "h"), mixed = TRUE)
```

```python
from summclust import summclust  # add the python directory to the module path
fit = summclust(data, "y", "x", ["z"], ["g", "h"], mixed=True)
```

Omit `mixed`, or use `FALSE`/`False`, for ordinary CV3. See [the mixed estimator notes](docs/MIXED_CV31.md) for component normalization, unavailable intersection diagnostics, and validation scope.

## Scope and behavior

- Default CV3 is the uncentered delete-cluster estimator. The `jackknife` option requests centered CV3J in one-way mode; it is intentionally rejected in two-way mode.
- For two-way `absorb()`, every absorbed group must be nested inside a G-H intersection. Use `fevar()` for crossed fixed effects.
- Native two-way estimation rejects full or deleted models in which the focal slope is not numerically identified. Rank loss affecting only nuisance coefficients remains allowed.
- One-way mode retains the historical near-zero-coefficient flag and dropped-subsample policy for compatibility; `rankcheck` supplies actual numerical rank diagnostics. These are separate concepts.
- Partial leverage corrects a bug in 4.22. Version 5.1.1 also preserves double-precision outcomes and corrects absorbed-model CV1 degrees of freedom. Results may therefore differ from old releases.
- Native two-way calls use private temporary Mata objects. The legacy one-way path still exposes working Mata objects, as documented in the audit.

See [CHANGELOG.md](CHANGELOG.md), [the integration audit](docs/AUDIT_2026-09-23.md), and [the output audit](docs/AUDIT_OUTPUT_2026-09-30.md) for the changes, checks, and limitations.

## Reproduce the checks

On Windows, run from PowerShell:

```powershell
./tests/run_tests.ps1 -StataExe 'C:/Program Files/Stata18/StataSE-64.exe'
```

The runner copies the package and test scripts to a new temporary directory, runs all nine suites and checks plain-text and rendered SMCL output, and prints its log directory. It checks actual success markers because Stata batch exit codes alone do not reliably detect do-file errors. Local installation is tested in an isolated PLUS directory.

On other platforms, copy `Stata/*` and `tests/*.do` to a scratch directory and run the nine do-files listed in `tests/run_tests.ps1` from there. The PowerShell runner itself is Windows-specific.

The targeted three-language mixed and string fixed-effect checks and their portable manual commands are documented in [validation/README.md](validation/README.md). They test full covariance matrices, inference, unchanged defaults, actual skipped solver calls, and cluster diagnostics aligned by identifiers. Run them from the repository root; generated files go to the ignored `output/` directory.

## Background and license

MacKinnon, Nielsen, and Webb, [Leverage, Influence, and the Jackknife in Clustered Regression Models: Reliable Inference Using summclust](https://arxiv.org/abs/2205.03288).

Distributed under the [MIT License](LICENSE).
