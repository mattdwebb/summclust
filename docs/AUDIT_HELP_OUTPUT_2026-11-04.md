# Stata help and max-standard-error output

James reported broken help formatting and identical CV3max and CV31max rows
in the nlswork example. Both reports were reproduced in StataNow/SE 18.5.

## Help formatting

The first three example commands and the final two-command example lacked
paragraph terminators. Stata's SMCL translator reproduced the collapsed text
and exposed directives after `keep if inrange(age,20,40)`. James's five
`{p_end}` additions fix this defect. The revised help incorporates those fixes,
closes the other paragraphs explicitly, simplifies the explanatory prose,
and documents the existing mixed component selectors and tie masks.

This behavior follows Stata's paragraph-mode rules; see the official
[SMCL manual](https://www.stata.com/manuals/psmcl.pdf).
The certification now renders the actual help file with Stata and checks
that example commands and subsequent sections remain separate, without raw
SMCL directives. The original help fails this check.

## Which standard error is reported?

The example is:

```stata
webuse nlswork, clear
keep if inrange(age,25,35)
generate vismin = inrange(race,2,3)
summclust hours vismin south, cluster(age ind_code) fevar(age birth_yr year ind_code) nograph
summclust hours vismin south, cluster(age ind_code) fevar(age birth_yr year ind_code) mixed nograph
```

There are 13,754 observations, 11 ages, 12 industries, and 132 intersections.
The focal coefficient is 1.054672 at the displayed precision.

| Focal standard error | Value |
|---|---:|
| One-way CV3, age | 0.136906654195236 |
| One-way CV3, industry | 0.521627594608576 |
| Conventional two-way CV3 | 0.497498100620775 |
| Conventional two-way CV31 | 0.498866581344820 |
| CV3max and CV31max | 0.521627594608576 |

Both maxima select H (`ind_code`), with tie mask 4, so there is no estimator
substitution error. The two-way standard errors differ; the industry CV3
component exceeds both. Inference retains min(G,H)-1 = 10 degrees of freedom.
CV1max also selects industry, with standard error 0.420219654623072.

Previously this was identified only by the terse `Max selectors` line below
the diagnostic summary. The revised display places the explanation directly
below the regression table:

```text
CV1max uses one-way CV1 clustered by ind_code (H).
CV31max uses one-way CV3 clustered by ind_code (H).
```

Conventional two-way winners are identified explicitly. Ties list all named
components within the existing relative tolerance of 1e-10. The selected
component still identifies the actual maximum; the tolerance does not change
selection or inference. In mixed mode, the one-way candidates remain CV3.

The existing `r(cv1max_component)`, `r(cv3max_component)`, and
`r(cv31max_component)` results already provide the selected component.
Their values are `two-way`, `G`, or `H`; `r(cluster_g)` and `r(cluster_h)`
give the variable names. Existing `r(cv#max_ties)` results identify ties.

## Validation and scope

The ado changes affect display text and local strings only. No Mata routines,
covariance calculations, inference calculations, stored results, dependencies,
or version numbers changed. R, Python, and historical releases were untouched.

The nlswork comparison reloads the same estimation data and fixes Stata's
sort seed before each fit. Original and revised fits have zero differences
in every returned matrix and scalar, for ordinary and mixed modes; returned
macros also agree. Fixing the sort seed avoids floating-point differences
from arbitrary ordering within tied sort keys.

Validation includes all eight existing Stata certification suites, the
targeted mixed-estimator certification, and 14 plain-text/SMCL output cases.
The output checks cover G and H winners, swapped dimensions, conventional
CV3/CV31 winners, ties, long names, negative conventional variance, and the
help-rendering regression. Source labels are compared with full-precision
component diagonals, and text output is compared with Stata-rendered SMCL.
