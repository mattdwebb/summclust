# Using summclust in Python

`summclust` helps you see whether a regression result depends heavily on a few
clusters, such as schools, firms, villages, or states. It estimates your linear
regression, reports clustered standard errors, and shows what happens to the
coefficient of interest when each cluster is left out.

[Repository home](../README.md) · [Stata guide](../Stata/README.md) · [R guide](../R/README.md)

## Get started

Download and extract this repository, or clone it. Use Python 3.10 or later
with NumPy; pandas is not required. If needed, install NumPy from a terminal:

```sh
python -m pip install numpy
```

Set your working directory to the extracted repository folder—the folder
containing `Stata`, `R`, `python`, and `examples`. In a notebook, this means
opening the notebook from that folder or changing to it first. The complete
[example script](../examples/school_program.py) can also be run from a terminal:

```sh
python examples/school_program.py
```

## Example: a school program and test scores

The shared [example data](../examples/school_program.csv) are **fictional**.
There are 1,710 student observations from 24 schools over five years. Schools
differ in size, and some introduce a program during the period.

| Column | Meaning |
|---|---|
| `test_score` | Student's test score, in points |
| `program` | 1 if the school has the program that year, 0 otherwise |
| `baseline_score` | Student's baseline test score |
| `school` | School identifier, stored as a string |
| `year` | Year of observation |

We estimate the score difference associated with the program, controlling for
baseline scores and year effects. We cluster by school to allow observations
within a school to have correlated errors.

```python
import sys
import numpy as np
sys.path.insert(0, "python")
from summclust import read_csv, summclust

data = read_csv("examples/school_program.csv", string_columns=["school"])
fit = summclust(
    data, outcome="test_score", focal="program", controls=["baseline_score"],
    cluster="school", fevar="year",
)
```

`focal` names the variable whose coefficient you want to investigate. `controls`
lists the other regressors, and `fevar` adds categorical fixed effects.
`summclust` estimates the regression itself; no separate regression call is needed.

## Read the regression results

Python returns results in a dictionary. The following statements print a small
table for the coefficient on `program`:

```python
print(f"{fit['N']} observations in {fit['G']} schools")
for method in ("cv1", "cv3"):
    r = fit["inference"][method]
    print(f"{method.upper()}: estimate={r['beta']:.3f}, SE={r['se']:.3f}, "
          f"p={r['p']:.3f}, 95% CI=[{r['lci']:.3f}, {r['uci']:.3f}]")
```

Output:

```text
1710 observations in 24 schools
CV1: estimate=3.426, SE=0.871, p=0.001, 95% CI=[1.625, 5.227]
CV3: estimate=3.426, SE=1.059, p=0.004, 95% CI=[1.235, 5.618]
```

The fitted program difference is about **3.43 score points**, after accounting
for baseline scores and year effects. CV1 is the usual clustered standard-error
calculation. CV3 uses the estimates obtained by leaving out each school in turn.
They give the same coefficient but different measures of uncertainty. Here CV3
gives a larger standard error and a wider confidence interval.

These fictional results illustrate how to read the output. A causal interpretation
in your own study depends on the research design, not on the choice of standard error.

## See which schools matter most

Partial leverage measures a school's share of the variation in `program` that
remains after accounting for the other regressors. The shares sum to one. With
24 schools, an equal share would be about `1 / 24 = 0.042`; this is a useful
comparison, not a cutoff for deciding whether a school is a problem.

Show the five schools with the largest partial leverage:

```python
print("School       N   Partial leverage   Estimate without school")
for j in np.argsort(-fit["partlev"])[:5]:
    print(f"{fit['cluster_ids'][j]:<10} {int(fit['ng'][j]):3d} "
          f"{fit['partlev'][j]:18.3f} {fit['betanog'][j]:25.3f}")
```

Output:

```text
School       N   Partial leverage   Estimate without school
School 24  300              0.174                     4.220
School 13  100              0.111                     3.276
School 04   80              0.088                     3.469
School 05  100              0.062                     3.360
School 08   50              0.055                     3.368
```

School 24 contributes about 17.4% of the remaining program variation. Leaving it
out raises the estimated coefficient from 3.43 to 4.22. School 13 also has high
partial leverage, but leaving it out moves the coefficient only to 3.28. A school
can contribute a lot of identifying variation without greatly changing the result.

Use these diagnostics to understand important clusters and report sensitivity.
A large value is not, by itself, a reason to delete a school. The separate
`fit["lever"]` field measures leverage for the regression as a whole; it is not
the same as partial leverage and does not sum to one.

## Adapt the example to your study

Use your column names in `outcome`, `focal`, `controls`, and `cluster`. At least
one control is required. Regression columns must be numeric. Cluster identifiers
and fixed-effect identifiers may be numeric or strings; include string columns
in `read_csv(..., string_columns=[...])` when loading a CSV.

To restrict the sample, pass a Boolean condition over the original rows:

```python
recent = summclust(data, "test_score", "program", ["baseline_score"], "school",
                  fevar="year", sample=data["year"] >= 2021)
```

To add school effects as well as year effects:

```python
school_fe = summclust(data, "test_score", "program", ["baseline_score"], "school",
                      fevar=["school", "year"])
```

This changes the regression: the program coefficient now uses variation within
schools over time. String school names work directly. For this one-way model,
`absorb="school", fevar="year"` is another way to account for school and year
effects. An absorbed group must lie within one cluster in each clustering dimension.

For two-way clustering, supply two names, for example `cluster=["school", "year"]`.
The reported inference uses the largest of the conventional two-way standard error
and the two one-way standard errors. This example has only five year clusters,
which limits the information available for two-way inference. Choose clustering
dimensions based on the sources of dependence in your study.

The optional `mixed=True` uses CV31 in place of CV3 for two-way clustering. It can
save computation but changes the standard-error estimator. It is never selected
automatically. In mixed mode, leave-one-intersection coefficients are not computed.

Missing model values, nonfinite numeric values, and empty string identifiers are
excluded before estimation. The input data are not modified. Check
`fit["N"]`, `fit["sample_mask"]`, and `fit["warnings"]` when adapting a specification.

See the [Python reference](REFERENCE.md) for all options, saved results, and
details of two-way inference. The [paper](https://arxiv.org/abs/2205.03288) explains
the methods behind the diagnostics.
