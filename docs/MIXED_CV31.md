# Mixed estimator (CV31), unreleased development integration

The explicit Stata `mixed` option and R/Python `mixed = TRUE`/`mixed=True`
select

\[
V_{31}=V_{3,G}+V_{3,H}-V_{1,I},\qquad I=G\cap H.
\]

It is off by default. No cluster-count or runtime threshold changes the method.
The estimator may save computation; no claim of superior inference is made.

The primary CV3 components use the existing uncentered delete-cluster factor
`(J-1)/J`. The intersection component reuses ordinary CV1, including the
factor `I*(N-1)/((I-1)*(N-k_eff))`. Effective model degrees of freedom,
estimation sample, transformed design, and fitted coefficients are unchanged.

## Inference and returns

`cv31`/`cv31_conventional` are raw full covariance matrices. Conventional focal
variance and SE use `cv31_conventional_var` and `cv31_conventional_se`.
Max inference compares the conventional focal variance with the G and H
CV3 focal variances and uses `min(G,H)-1` Student-t degrees of freedom.
Explicit `cv31_max_var/se/t/p/lci/uci` returns identify this separate inference.
Existing-style `cv31se/t/p/lci/uci` aliases also refer to max inference.
The selected estimator is identified by `selected_estimator`.

There is no full elementwise-max covariance matrix. A negative conventional
variance remains negative and its SE is unavailable; max inference remains
separate. No eigenvalue-adjustment option or positive floor is introduced.
The display labels distinguish CV31 (mixed) and CV31 (mixed, max).

## Skipped work and diagnostics

The intersection path accumulates scores and inexpensive diagnostics but
does not solve deleted-intersection regressions or construct CV3_I.
The counter `intersection_delete_one_solves` is zero under mixed.
Sizes, ordinary leverage, partial leverage, and applicable G-star statistics
remain available. Intersection omit-one coefficients, related summaries, and
CV3 matrices are missing (`.` in Stata, `NA` in R, `NaN` in Python).
Ordinary CV3 return fields are never filled with CV31.

An explicit `rankcheck` calculates deleted-design ranks without coefficient
solves and may add work. Existing G/H tables and graphs retain their normal
meaning. One-way mixed and centered two-way CV3J are rejected.

## Verification

The portable checks in `validation/README.md` cover direct score normalization,
independent component construction, negative conventional variance, identical
and nested partitions, swapped dimensions, missing/string identifiers, fixed
effects, supported absorption, nuisance rank deficiency, focal changes and
row permutations. For identical dimensions the identity is
`2*CV3_G-CV1_G`, not one-way CV3. Frozen R/Python baseline test fixtures
verify exact unchanged-default behavior.

Executed validation on 2026-10-30 passed all eight preexisting Stata suites
and formatting checks, targeted Stata mixed checks, 36 original and 42 mixed
R assertions, 25 Python unit tests, and 20 default-parity outcomes in each
port. The shared synthetic comparison checked 7,498 numeric cells with a
maximum absolute difference of approximately `1.08e-10`.

Separate empirical checks covered the paper's NLS work-hours software example
and seven two-way tsetse specifications. Maximum empirical full-covariance
differences were about `2.94e-9`; focal inference differed by at most `1.99e-11`.
A small nuisance fixed-effect covariance discrepancy was reproduced in the
untouched baseline and documented in the development evidence. Empirical
datasets, manuscript attachments, historical outputs, and bootstrap draws
are not included in this repository update.

Version numbers and historical releases are unchanged. This is a development
source update; it does not create a tagged release.
