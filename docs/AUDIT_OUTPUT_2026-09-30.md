# Two-way output audit and revision — 30 September 2026

The revised local summclust output is ready for James to retest. It restores the original standalone twowayjack table presentation while preserving the justified numerical corrections. All eight current certification suites and eight text/SMCL formatting scenarios passed in fresh runs. Nine empirical models passed independent regression checks, and every stored scalar, macro and matrix was identical before and after the presentation change. No Monte Carlo experiment, release publication or communication to James was performed.

## Sources, versions and execution

The destination is `C:/Users/mattw/Dropbox/software/summclust`, on `reorganize-stata-r`, starting at commit `4c4409e` (summclust 5.1.1, 23 September 2026). This is an **unreleased presentation revision**, identified by a second header line dated 30 September 2026; the version number and packaged delegate version are not republished.

The genuine standalone reference is `C:/Users/mattw/Dropbox/research/multi-jack/twowayjack/twowayjack.ado`, version **1.2, 24 June 2024**, containing its own `mwsummclust` based on **4.210, 6 June 2023**. Its SHA256 is `59fd15aa4d5b0193162ee9d31388021deb603f463b3dd2de40c980371924d040`. It is byte-identical to `influence/summcluster/summcluster_twoway/twowayjack_legacy.ado`. Historical `submission_files/nlswork-twoway.out` and `mwearn-v2.out` supply the original polished output. The packaged **twowayjack 2.0.0, 23 August 2026**, delegates to native summclust and was used only for compatibility testing, not as an independent reference.

Fresh runs used **StataNow/SE 18.5 for Windows, revision 26 February 2025**. All historical sources were read or copied, not edited. Investigation, baseline and revised logs are under:

`C:/Users/mattw/OneDrive/Documents/ChatGPT/multi-jack/summclust_output_audit_20260930`

The `before` and `after` ado directories are selected explicitly with `adopath ++`; the drivers run from the audit root, which contains no ado file named summclust or twowayjack. Logs record `which ..., all` to expose shadowing. The standalone reproduction selects `reference/standalone` explicitly. `source_sha256.json` records code, data and historical-reference hashes. Built-in `webuse` failed to connect, so the nlswork data were fetched from the official `https://www.stata-press.com/data/r18/nlswork.dta` endpoint and cached under `data/`.

## Recovered nlswork model

The complete specification comes from `submission_files/nlswork-twoway.do` and is corroborated by `examples/nlswork/nlswork.do`. Its historical output matches all values quoted in the request. No specification was inferred from the reported coefficient.

```stata
webuse nlswork, clear
keep if inrange(age,25,35)
gen vismin = inrange(race,2,3)
twowayjack hours vismin south, fevar(age birth_yr year ind) cluster(age ind)
```

For the revised command, replace `twowayjack` by `summclust`. Stata expands `ind` to `ind_code`. The common complete-case sample has **13,754 observations**, **11 age clusters**, **12 industry clusters**, and **132 intersections**. All explicit row, column and intersection deletions were refitted independently with `regress hours vismin south i.age i.birth_yr i.year i.ind_code` on the same sample.

## Discrepancy inventory

| Class | Finding and resolution |
| --- | --- |
| Presentation changed | Restored `Regression Output`, `s.e.`, `Coeff`, `Sd. Err.`, `t-stat`, `P value`, `CI-lower`, `CI-upper`, and the `CV1max`/`CV3max` labels using the standalone `matlist` column/row specifications. |
| Presentation changed | Restored `Coefficients of Variation, G, and G*`, all six headings, actual cluster names and `intersect`, four decimals for the first four columns, integer G, two-decimal Gstar, and the historical spacing/separators. |
| Presentation changed | Optional cluster and coefficient tables now expose column headings. Max selectors identify actual cluster variables. Notes and warnings wrap at word boundaries. Ties are identified without changing the existing selector/tie returns. |
| Presentation changed | At width 80, a 32-character cluster name is shown completely on its own label line above aligned values; at width 120 the full name fits in the expanded row-label column. No silent abbreviation. |
| Proven numerical correction retained | Historical CV1 CI used **I−1=131** degrees of freedom, despite using **min(G,H)−1=10** for its p-value. The standalone source explicitly calls `invttail(rows(betag)-1,.025)` after its intersection run. Current CV1 and CV3 inference consistently use df=10. |
| Proven numerical correction retained | Historical partial leverage used residuals of **south**, not **vismin**. The standalone source's `reg ``varname'' ``tempcons'' CTEMP, noc` loses the undefined doubly expanded macros and treats the first control as the dependent variable. The fresh standalone run reports `e(depvar)=south`; an independent regression of south on the FE indicators reproduces all three historical values. Current FWL residuals of vismin on south and all FEs reproduce the native values. |
| Proven wrapping defect fixed | The unmodified 5.1.1 text log at width 80 splits `diagno` and `> stic` across lines. Its raw SMCL contains the unbroken word, but its text rendering hard-wraps it. Thus this is a real Stata output issue, not merely an email effect. The original email itself was not available; no claim is made about additional email reflow. |
| Missing source qualified | The submission's earnings file named `minwage.csv` is not present under that specification: the local file with that name contains `utothrs`, not `hrlyearn`. The available `minwage_1.csv` has 28,599 earnings observations and `age_12`; renaming it to `age` and applying the recorded `age=age-1` transformation reproduces the historical coefficient and both max SEs to their displayed precision. Byte identity with James's original data file is unproved. |
| Scope limit | No historical SMCL file for James's original run was found. The historical text `.out` files were preserved, and fresh standalone, baseline and revised text/SMCL logs were generated. Stata 13 and other platforms were not executed. |

No new estimator correction was necessary: this patch changes presentation and documentation, while retaining independently validated current calculations. The one-way program and all Mata numerical routines are unchanged, apart from the top-level presentation-revision comment.

## Numerical reconciliation

For nlswork, the independent slope is **1.054671801894**; max standard errors are **0.420219654623** (CV1) and **0.521627594609** (CV3). The correct 97.5th percentile is **2.228138851986273** for t(10); the historical CV1 calculation used **1.978238539230381** for t(131).

| Quantity | Historical | Retained, validated result |
| --- | --- | --- |
| CV1 95% CI | [0.223377, 1.885967] | **[0.118364, 1.990980]** |
| CV3 95% CI | [−0.107587, 2.216931] | **[−0.107587, 2.216931]** |
| Partial-leverage CV, age | 0.0927004761 | **0.1022915673** |
| Partial-leverage CV, ind_code | 1.1848831920 | **1.1919547450** |
| Partial-leverage CV, intersection | 1.1557385531 | **1.1647102492** |

The intended df rule is also stated in the local paper `drafts/MNW-twojack-v44.tex`, around lines 832–836 and 867–870. The new code checks p-values and both endpoints against full-precision t calculations; it does not fit a critical value to rounded historical endpoints. Direct `regress, vce(cluster ...)` checks each CV1 component. Repeated direct OLS fits check all CV3 deletion components. `predict, residuals` and `predict, hat` independently check partial leverage, leverage and G*(0).

## Empirical checks

The tsetse specifications use all seven outcomes and the exact controls/cluster encodings from `examples/tsetse/tsetse.do`, with `precolonial.dta`. The minimum-wage model uses the recorded submission equation and the explicitly qualified available earnings extract. These are deterministic checks, not simulation replications.

| Model | N | G/H/I | Coefficient | CV1max SE | CV3max SE | Direct references |
| --- | ---: | --- | ---: | ---: | ---: | --- |
| minwage | 28599 | 12/10/120 | 0.293407583 | 0.096488176 | 0.149218207 | Pass |
| nlswork | 13754 | 11/12/132 | 1.054671802 | 0.420219655 | 0.521627595 | Pass |
| tsetse_animal | 484 | 44/44/142 | -0.230974731 | 0.042233545 | 0.048276640 | Pass |
| tsetse_central | 467 | 44/44/140 | -0.074593753 | 0.036496703 | 0.041843541 | Pass |
| tsetse_female | 315 | 43/41/112 | 0.205707393 | 0.060367845 | 0.069576599 | Pass |
| tsetse_intensive | 485 | 44/44/142 | -0.090478344 | 0.030185371 | 0.034599821 | Pass |
| tsetse_ln_popd_murdock | 398 | 43/42/121 | -0.744617019 | 0.234392746 | 0.271882339 | Pass |
| tsetse_plow | 484 | 44/44/142 | -0.057212502 | 0.031811035 | 0.038271675 | Pass |
| tsetse_slave | 446 | 44/43/139 | 0.101328510 | 0.040475813 | 0.051913079 | Pass |

Across the nine models, every full-precision r() snapshot is byte-identical before/after. The largest direct deletion-slope discrepancy is **3.58e−12**; the largest CV3 component discrepancy is **5.79e−13**. These small differences are numerical roundoff. The numerical test tolerances were chosen independently of display rounding and no existing certification tolerance was relaxed.

## Fresh certification results

All eight suites passed: `certify_summclust`, `certify_summclust_twoway`, `certify_summclust_twoway_api`, `certify_twowayjack_delegate`, `certify_integration`, `certify_output_numerics`, `certify_output_format`, and `verify_net_install_twoway`.

The separate PowerShell formatting checks passed for ordinary output, combined options, 32-character names at widths 80/120, negative conventional variances, two clusters, identical/tied dimensions and table suppression above 52 clusters. They inspect headings, number precision, integer counts, separators, widths, complete labels and warning/option text. Stata's `smcl2log` translator renders each SMCL fixture at the matching line width; rendered and plain-text output agree after trimming line-end spaces. Certification requires actual output success markers, not process exit codes alone. Final logs are in `final_suites/`; earlier development logs are not cited as final certification evidence.

## Reproducible commands

Run the package tests from PowerShell (output directory must be new):

```powershell
$repo = 'C:/Users/mattw/Dropbox/software/summclust'
& "$repo/tests/run_tests.ps1" -StataExe 'C:/Program Files/Stata18/StataSE-64.exe'
```

To reproduce the project-specific audit using the cached reference files and data, work from the dated audit directory. The drivers write only to its dedicated before/after/reference directories; they do not write into the research source tree. To retain this audit's generated logs, copy the audit directory to a new sibling directory before rerunning.

```powershell
Set-Location 'C:/Users/mattw/OneDrive/Documents/ChatGPT/multi-jack/summclust_output_audit_20260930'
$stata = 'C:/Program Files/Stata18/StataSE-64.exe'
Start-Process $stata -ArgumentList '/e do legacy_audit.do' -WorkingDirectory $PWD -WindowStyle Hidden -Wait
Start-Process $stata -ArgumentList '/e do nlswork_audit.do before' -WorkingDirectory $PWD -WindowStyle Hidden -Wait
Start-Process $stata -ArgumentList '/e do nlswork_audit.do after' -WorkingDirectory $PWD -WindowStyle Hidden -Wait
Start-Process $stata -ArgumentList '/e do paper_examples.do before' -WorkingDirectory $PWD -WindowStyle Hidden -Wait
Start-Process $stata -ArgumentList '/e do paper_examples.do after' -WorkingDirectory $PWD -WindowStyle Hidden -Wait
python build_audit_report.py
```

## Before/after evidence: nlswork

Historical standalone output (preserved source `.out`):

```text
Regression Output

    s.e. |      Coeff   Sd. Err.   t-stat  P value    CI-lower    CI-upper
---------+----------------------------------------------------------------
  CV1max |   1.054672   0.420220   2.5098   0.0309    0.223377    1.885967
  CV3max |   1.054672   0.521628   2.0219   0.0708   -0.107587    2.216931
--------------------------------------------------------------------------

Coefficients of Variation, G, and G*

 dimension |       Ng   Leverage  Partial L.   beta no g        G      Gstar
-----------+----------------------------------------------------------------
       age |   0.0987     0.1813      0.0927      0.0431       11      10.90
  ind_code |   1.1815     0.8823      1.1849      0.1565       12       5.21
 intersect |   1.1507     0.8925      1.1557      0.0173      132      56.26
----------------------------------------------------------------------------
```

Unmodified summclust 5.1.1, freshly executed:

```text
Focal inference using MacKinnon-Nielsen-Webb max standard errors
     s.e. |Coefficient  Std. err.        t    P>|t| [95% conf.  interval]
----------+-----------------------------------------------------------------
  CV1 max |   1.054672   0.420220   2.5098   0.0309   0.118364   1.990980
  CV3 max |   1.054672   0.521628   2.0219   0.0708  -0.107587   2.216931
----------+-----------------------------------------------------------------
Max selectors: CV1=H; CV3=H
Conventional focal CV1: variance=   0.153264  s.e.=   0.391489
Conventional focal CV3: variance=   0.247504  s.e.=   0.497498
Conventional matrices use V_G + V_H - V_intersection and are not used for the ma
> x-se inference above.

Dimension-specific diagnostics

           G |    0.099     0.181     0.102     0.043    11.000    10.896
           H |    1.181     0.882     1.192     0.157    12.000     5.212
           I |    1.151     0.892     1.165     0.017   132.000    56.259
note: two-way diagnostic graphs are not drawn; use the dimension-specific diagno
> stic matrices in r().
```

Revised output, freshly executed:

```text
Regression Output

    s.e. |      Coeff   Sd. Err.   t-stat  P value    CI-lower    CI-upper
---------+----------------------------------------------------------------
  CV1max |   1.054672   0.420220   2.5098   0.0309    0.118364    1.990980
  CV3max |   1.054672   0.521628   2.0219   0.0708   -0.107587    2.216931
--------------------------------------------------------------------------

Coefficients of Variation, G, and G*

 dimension |       Ng   Leverage  Partial L.   beta no g        G      Gstar
-----------+----------------------------------------------------------------
       age |   0.0987     0.1813      0.1023      0.0431       11      10.90
  ind_code |   1.1815     0.8823      1.1920      0.1565       12       5.21
 intersect |   1.1507     0.8925      1.1647      0.0173      132      56.26
----------------------------------------------------------------------------
Inference: t(10) with df = min(G,H)-1; 95% confidence intervals.
Max selectors: CV1=ind_code; CV3=ind_code.
Conventional focal CV1: variance=   0.153264  s.e.=   0.391489
Conventional focal CV3: variance=   0.247504  s.e.=   0.497498
Conventional matrices use V_G + V_H - V_intersection. The inference table uses
max standard errors, not these matrices.
Note: two-way diagnostic graphs are not drawn; use the dimension-specific
diagnostic matrices in r().
```
