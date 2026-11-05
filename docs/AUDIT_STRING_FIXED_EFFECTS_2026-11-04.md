# String fixed-effect identifiers — 4 November 2026

The current Stata, R, and Python sources accept numeric or string identifiers
for both `fevar` and `absorb`, in one-way, two-way, and mixed models. R also
accepts factors. No option, estimator formula, or version number was added.
Historical releases and frozen test baselines were not changed.

Stata's one-way path codes string categories after sample selection; numeric
identifiers retain their existing path. R's one-way string `fevar` rejection
was removed. Python already accepted strings, but now retains nonempty labels
such as `"NaN"` and `"Inf"` while excluding actual numeric missing values, even
in ordinary Python lists containing both strings and numbers.

String fixed-effect coefficient names use observed sorted category positions.
Long names use a shortened prefix and an effect-position suffix to fit Stata's
32-character name limit. Numeric coefficient naming is unchanged. Absorption
still requires nesting within every clustering dimension. Input data and row
order are preserved.

## Validation

Tests ran on Windows using StataNow/SE 18.5, R 4.5.3, and Python 3.14 with NumPy.
Stata 13 remains the declared compatibility level; this run did not test an
older Stata executable.

- All nine Stata certification suites passed, including the new string
  fixed-effect suite. All 14 plain-text/SMCL formatting scenarios and help
  rendering checks passed.
- Stata's focused suite passed 12 numeric/string equivalence scenarios,
  including absorption, multiple effects, mixed numeric/string effects,
  sample restrictions, missing identifiers, long names, nesting errors,
  two-way and mixed estimation, and graph generation.
- An additional before/after Stata comparison covered eight numeric-input
  specifications. Every returned matrix, scalar, macro, and matrix stripe
  was exactly unchanged with a fixed sort seed.
- R passed 36 original, 42 mixed, and 32 string/factor checks. All 20 frozen
  default fixture outcomes matched exactly.
- Python passed 33 unit tests, including eight new string fixed-effect tests.
  All 20 frozen default fixture outcomes matched exactly.
- Nine shared string fixed-effect specifications compared one-way, two-way,
  and mixed models using explicit effects, absorption, and both together.
  Across 9,564 compared cells, the maximum absolute difference between Stata
  and the other languages was `4.884981308350689e-15`, with no failures.
  Coefficients and covariance matrices were aligned by coefficient name;
  cluster diagnostics were aligned by actual cluster identifiers.
- The existing eight shared mixed specifications also passed: 7,498 cells,
  maximum absolute difference `1.0800249583553523e-10`, no failures.

Run `tests/run_tests.ps1` for the Stata suites and
`tests/run_mixed_tests.ps1` for the R/Python certifications and shared
mixed/string comparisons. See [validation/README.md](../validation/README.md)
for executable parameters and portable commands. R's certification includes
the new factor checks; Python's certification discovers the new test module.
These runners require fresh output directories and check Stata's actual
success markers rather than relying on its process exit status.
