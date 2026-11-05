# summclust in R

`summclust` helps you see how much a regression result depends on particular
schools, firms, regions, or other clusters. It estimates the regression, reports
cluster-robust standard errors, and shows what happens to your coefficient of
interest when each cluster is left out.

This guide uses the same [fictional school-program data](../examples/school_program.csv)
as the [Stata](../Stata/README.md) and [Python](../python/README.md) guides.

## Get started

Download this repository using **Code → Download ZIP**, then unzip it. In R or
RStudio, set your working directory to the top-level folder containing `R`,
`Stata`, `python`, and `examples`. For example:

```r
setwd("path/to/summclust-main")
source("R/summclust.R")
```

Only base R is needed; there are no additional packages to install. You can run
the complete example with `source("examples/school_program.R")`, or follow the
steps below.

## Example: a school program and test scores

The dataset contains 1,710 observations from 24 schools over five years.
`program` is 1 when a school participates in the program in that year and 0
otherwise. We study the score difference associated with participation,
controlling for baseline scores and year effects. Standard errors are clustered
by school. These data are fictional; all results below are for illustration.

```r
d <- read.csv("examples/school_program.csv")

fit <- summclust(
  d,
  outcome = "test_score",
  focal = "program",
  controls = "baseline_score",
  cluster = "school",
  fevar = "year"
)

options(digits = 5)
print(fit)
```

Output:

```text
summclust: 1710 observations; 24 clusters; inference df = 23
      beta      se      t          p    lci    uci
cv1 3.4264 0.87067 3.9353 0.00066048 1.6252 5.2275
cv3 3.4264 1.05944 3.2341 0.00366654 1.2348 5.6180
```

The estimated program coefficient is **3.43 score points**, after accounting
for baseline scores and year effects. `lci` and `uci` are the lower and upper
limits of the 95% confidence interval.

Both rows report the same coefficient. **CV1** uses conventional clustered
standard errors. **CV3** uses a cluster jackknife calculation based on leaving
out one school at a time. In this example CV3 gives a larger standard error
and a wider interval: approximately **1.23 to 5.62 points**. The choice of
standard-error calculation does not establish whether the relationship is causal.

`focal` names the variable whose coefficient you want to examine. Supply at
least one control variable. `summclust` fits the model directly; you do not
need to run `lm()` first.

## Which schools matter most?

Partial leverage measures each school's share of the variation in `program`
remaining after the other regressors are accounted for. The following code
lists the five schools with the largest shares, alongside the estimate obtained
when each is excluded:

```r
diagnostics <- data.frame(
  school = names(fit$ng),
  N = as.numeric(fit$ng),
  partial_leverage = as.numeric(fit$partlev),
  estimate_without_school = as.numeric(fit$betanog)
)
diagnostics$change <- diagnostics$estimate_without_school - fit$beta_full
top <- diagnostics[order(diagnostics$partial_leverage, decreasing = TRUE), ]
print(head(top, 5), row.names = FALSE, digits = 4)
```

Output:

```text
    school   N partial_leverage estimate_without_school   change
 School 24 300          0.17406                   4.220  0.79384
 School 13 100          0.11050                   3.276 -0.15008
 School 04  80          0.08832                   3.469  0.04263
 School 05 100          0.06182                   3.360 -0.06642
 School 08  50          0.05549                   3.368 -0.05842
```

School 24 provides about **17.4%** of this variation. An equal share across
24 schools would be about **4.2%**; that is a comparison point, not a cutoff.
Leaving School 24 out raises the estimated coefficient from 3.43 to 4.22 points.

The share measures potential influence; the change in the coefficient shows
actual sensitivity in this regression. A large value is a reason to understand
that school's data and report sensitivity, rather than automatically exclude it.
Ordinary leverage is available in `fit$lever`; it describes a school's role in
the regression design as a whole, while partial leverage focuses on `program`.

To save the diagnostics for every school:

```r
write.csv(diagnostics, "school_diagnostics.csv", row.names = FALSE)
```

## Adapt it to your own data

Use a data frame with numeric outcome and regression variables. Cluster and
fixed-effect identifiers can be numeric, character, or factors. String labels
such as `"School 24"` do not need conversion to numeric codes.

| What you want to change | Argument |
|---|---|
| Study a different coefficient | Set `focal` to that column's name |
| Add several controls | `controls = c("baseline_score", "age")` |
| Include school and year effects | `fevar = c("school", "year")` |
| Restrict the sample | `sample = d$year >= 2021` |
| Cluster by school and year | `cluster = c("school", "year")` |
| Check whether deleting a cluster makes some coefficients unidentified | `rankcheck = TRUE` |

These are changes to the call above; the named columns must exist in your data.
Missing selected variables and empty string identifiers are excluded. Year
effects are treated as categories by `fevar`, even when `year` is numeric.

`absorb = "school"` is another way to account for school effects with one-way
clustering by school. Each absorbed group must lie within a cluster in every
clustering dimension. A school observed in several years therefore cannot be
absorbed when also clustering by year; use `fevar = "school"` for that model.
The program variable must still vary after your chosen fixed effects and controls
are accounted for.

## Two-way clustering and the mixed option

Use two clustering variables when your research design calls for them. For
example, `cluster = c("school", "year")` allows dependence within schools and
within years. This teaching dataset has only five years, so that specification
has few clusters in one dimension; it should not be treated as an assurance of
precise inference.

The printed two-way standard error is the largest of the conventional two-way
standard error and the two one-way standard errors for the coefficient of
interest. For a two-way fit, `fit$cv3max_component` identifies the selected component: `"two-way"`,
`"G"` (the first cluster variable), or `"H"` (the second).

Adding `mixed = TRUE` selects CV31, a different two-way jackknife-based
calculation that can be faster. It is optional, requires two cluster variables,
and is never selected automatically. In that case use `fit$cv31max_component`.
See the [reference](REFERENCE.md) for saved results and the precise distinctions
between these calculations.

## More help

- [Complete example script](../examples/school_program.R)
- [Arguments and saved-results reference](REFERENCE.md)
- [Shared example data and variable definitions](../examples/README.md)
- [Repository overview and research background](../README.md)
