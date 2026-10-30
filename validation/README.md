# Deterministic mixed-estimator certification

Run commands from the repository root. These tests use only the small synthetic
fixtures in `fixtures/`. They do not download data or run bootstrap/simulation
replications. Runtime dependencies are NumPy for Python, base R and its standard
packages, and Stata with built-in Mata functions.

```text
python validation/python_certification.py
Rscript --vanilla R/certify.R
```

Then run `do tests/test_stata_mixed.do` in Stata and generate R/Python comparisons:

```text
python validation/run_mixed.py
Rscript --vanilla R/run_mixed_fixtures.R
python validation/compare_mixed.py
```

On Windows, `tests/run_mixed_tests.ps1` runs the same sequence in a fresh scratch
directory, checks real Stata success markers, and reports the evidence location.
Its `-PythonExe`, `-RscriptExe`, and `-StataExe` parameters select installed
executables. The existing `tests/run_tests.ps1` separately runs the eight
ordinary Stata suites and output-format checks.

The new checks cover independent direct-score normalization, component
reconstruction, no intersection solve calls, max inference and Student-t
degrees of freedom, negative raw variances, exact default parity, identical and
nested partitions, FE/absorption, strings and missing observations, nuisance
rank deficiency, focal changes, and row/partition permutations.

`cases.csv` defines 20 existing fixture outcomes; eight successful two-way
specifications are also compared with `mixed`. Numeric covariance cells are
aligned by coefficient names and cluster diagnostics by actual identifiers.
The comparison tolerance is `2e-10 + 2e-8*max(abs(a),abs(b))`. Missing fields
must agree. Every required result is checked for presence.

Fresh outputs are written under `output/`, which is ignored by Git. Historical
empirical data, manuscript files, and release archives are not needed.
