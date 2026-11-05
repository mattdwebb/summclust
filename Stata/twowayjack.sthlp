{smcl}
{* *! version 2.0.0 23aug2026}{...}
{vieweralsosee "summclust" "help summclust"}{...}
{title:Title}

{pstd}
{cmd:twowayjack} {hline 2} Two-way cluster-robust inference and diagnostics{p_end}

{title:Syntax}

{phang}
{cmd:twowayjack} {it:depvar focalvar controls} {ifin},
{opt cluster(var1 var2)} [{it:options}]{p_end}

{pstd}
The first variable is the outcome, the second is the regressor of interest,
and the remaining variables are controls. At least one control and exactly
two clustering variables are required. Factor-variable notation and weights
are not supported.{p_end}

{title:Description}

{pstd}
{cmd:twowayjack} reports the same two-way CV1 and CV3 inference and cluster
diagnostics as {help summclust}. Use {cmd:if} or {cmd:in} to restrict the
sample. If you also specify {opt sample()}, observations must satisfy all
of the restrictions.{p_end}

{pstd}
The displayed {cmd:CV1max} and {cmd:CV3max} standard errors use the largest
of the conventional two-way variance and the two one-way variances for the
coefficient of interest. Inference uses min(G,H)-1 degrees of freedom.
The output names the component supplying each maximum.{p_end}

{title:Options}

{pstd}
{opt cluster(var1 var2)} specifies two numeric or string clustering variables.
The first is labelled G and the second H in the returned results.{p_end}

{pstd}
The options {opt fevar()}, {opt absorb()}, {opt gstar}, {opt rho()},
{opt sample()}, {opt table}, {opt regtable}, and {opt rankcheck} have the same
meaning as in {help summclust}. Both {opt fevar()} and {opt absorb()} accept
numeric or string identifiers; no {cmd:encode} step is needed. Empty string
identifiers and numeric missing values are excluded from the estimation sample.
Fixed effects specified with {opt absorb()}
must be nested within intersections of the two clustering dimensions.
Use {opt fevar()} for fixed effects that cross either dimension.{p_end}

{pstd}
No graphs are drawn. {opt nograph} suppresses the note about graphs.
{opt addmeans} is accepted but adds no summary statistics.
{opt jackknife} is not available with two-way clustering.{p_end}

{pstd}
To request the mixed estimator CV31, use {cmd:summclust} with {opt mixed}.
{cmd:twowayjack} does not accept {opt mixed}. Use {opt sample()} to express
any sample restriction when calling {cmd:summclust}.{p_end}

{title:Examples}

{phang2}{cmd:. webuse nlswork, clear}{p_end}
{phang2}{cmd:. generate vismin = inrange(race,2,3)}{p_end}
{phang2}{cmd:. twowayjack hours vismin south if inrange(age,25,35), cluster(age ind_code) fevar(age birth_yr year ind_code)}{p_end}

{pstd}Use the same sample and fixed effects with the mixed estimator:{p_end}
{phang2}{cmd:. summclust hours vismin south, cluster(age ind_code) fevar(age birth_yr year ind_code) sample(inrange(age,25,35)) mixed nograph}{p_end}

{pstd}
For your own data, a string variable {cmd:industry} containing categories
such as "Manufacturing" and "Retail" can enter directly as a fixed effect:
{p_end}
{phang2}{cmd:. twowayjack y treatment controls, cluster(firm year) fevar(industry) nograph}{p_end}

{title:Stored results}

{pstd}
{cmd:twowayjack} stores the two-way results documented in {help summclust}
in {cmd:r()}. Use {cmd:return list} to inspect them. The conventional covariance
matrices and the displayed max standard errors are stored separately.{p_end}
