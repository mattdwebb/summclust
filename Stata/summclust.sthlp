{smcl}
{* *! version 5.1.1 23sep2026}{...}
{vieweralsosee "summclust" "help summclust"}{...}
{hline}

{title:Title}

{pstd}
{cmd:summclust} {hline 2} One- and two-way cluster-robust variance estimation,
leverage, influence, and cluster-jackknife diagnostics{p_end}

{title:Syntax}

{phang}
{cmd:summclust} {it:depvar focalvar controls},
{opt cluster(clustvar)} [{it:options}]{p_end}

{phang}
{cmd:summclust} {it:depvar focalvar controls},
{opt cluster(clustvar1 clustvar2)} [{it:options}]{p_end}

{pstd}
The varlist must contain at least three variables. The first is the dependent
variable, the second is the regressor of interest, and all remaining
variables are controls. Factor-variable notation is not accepted in this
varlist; use {opt fevar()} for categorical fixed effects.{p_end}

{pstd}
Use {opt sample()} to restrict the estimation sample. Standard Stata
{cmd:if}, {cmd:in}, and weight syntax are not supported.{p_end}

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{opt cluster(varlist)}}one or two numeric or string cluster identifiers; required{p_end}
{synopt:{opt fevar(varlist)}}include categorical fixed effects with numeric or string identifiers{p_end}
{synopt:{opt absorb(varname)}}partial out one numeric or string fixed effect; see the nesting conditions below{p_end}
{synopt:{opt jack:knife}}also calculate one-way CV3J; not available with two-way clustering{p_end}
{synopt:{opt mixed}}use CV31 instead of CV3 with two-way clustering{p_end}
{synopt:{opt gstar}}calculate effective-cluster measures G*(0) and, when available, G*(1){p_end}
{synopt:{opt rho(#)}}also calculate G*(rho); 0.00001 <= # <= 1{p_end}
{synopt:{opt sam:ple(exp)}}restrict the estimation sample using a logical expression{p_end}
{synopt:{opt tab:le}}display cluster-by-cluster diagnostics; in two-way mode, display separate G and H tables{p_end}
{synopt:{opt add:means}}display harmonic, geometric, and quadratic means and ratios; one-way only{p_end}
{synopt:{opt reg:table}}display results for every explicit regression coefficient{p_end}
{synopt:{opt nog:raph}}suppress the one-way four-panel graph; no graphs are drawn with two-way clustering{p_end}
{synopt:{opt rank:check}}check whether deleting a cluster reduces the rank of the regression design{p_end}
{synoptline}

{title:Description}

{pstd}
{cmd:summclust} reports cluster sizes, leverage, partial leverage,
omit-one-cluster coefficient estimates, CV1 and CV3 inference, and optional
CV3J and effective-cluster calculations. By default it draws four scatter
plots combining cluster size, leverage, partial leverage, and the
omit-one-cluster coefficient for one-way clustering. With two cluster
variables it reports a focal-coefficient table and separate
G, H, and nonempty-intersection diagnostics.{p_end}

{pstd}
The focal coefficient is the coefficient on the second variable in the
varlist. Put the regressor you want to study in that position. CV3 uses
deviations of the omit-one-cluster estimates from the full-sample estimate.
With one-way clustering, {opt jackknife} also reports CV3J, which centers
those estimates on their own mean. The main one-way p values and confidence
intervals use G-1 degrees of freedom, where G is the number of clusters.{p_end}

{pstd}
With one-way clustering, {opt regtable} displays the full coefficient table
using CV3 standard errors, or CV3J standard errors when {opt jackknife}
is also specified.{p_end}

{pstd}
Stata 13 or later is required.{p_end}

{title:Two-way clustering}

{pstd}
Specify two variables in {opt cluster()}, for example {cmd:cluster(firm year)}.
G denotes the first clustering dimension, H the second, and I their nonempty
intersections. The default estimators are CV1 and CV3. For either estimator,
the conventional two-way covariance is{p_end}

{phang2}
V_two-way = V_G + V_H - V_I.{p_end}

{pstd}
Each CV1 component uses a finite-sample correction based on its own cluster
count. The conventional covariance can have negative diagonal entries.
These are retained, with missing conventional standard errors; no eigenvalue
adjustment is made.{p_end}

{pstd}
The focal-coefficient table reports {cmd:CV1max} and {cmd:CV3max}. Each uses
the largest of the conventional two-way variance and the two one-way
variances for the focal coefficient:{p_end}

{phang2}
V_max = max(V_two-way[1,1], V_G[1,1], V_H[1,1]).{p_end}

{pstd}
All two-way t tests and 95% confidence intervals use min(G,H)-1 degrees of
freedom, including when a one-way component supplies the maximum. The max
rule applies to individual coefficients. It does not define a covariance
matrix for joint hypothesis tests. {opt regtable} reports max standard
errors separately for each coefficient.{p_end}

{title:The mixed estimator}

{pstd}
{opt mixed} replaces CV3 by CV31, using the intersection CV1 component in
place of intersection CV3:{p_end}

{phang2}
V_CV31 = V_CV3,G + V_CV3,H - V_CV1,I.{p_end}

{pstd}
CV31 avoids calculating delete-intersection estimates and can therefore be
faster. It is used only when you specify {opt mixed}, requires two cluster
variables, and cannot be combined with {opt jackknife}.{p_end}

{pstd}
With {opt mixed}, the focal inference row is labelled {cmd:CV31max}.
Its variance is max(V_CV31[1,1], V_CV3,G[1,1], V_CV3,H[1,1]). The two one-way
candidates remain CV3, and the degrees of freedom remain min(G,H)-1.
The conventional CV31 variance and standard error are reported separately.{p_end}

{pstd}
Intersection sizes, leverage, partial leverage, and applicable G* measures
remain available with {opt mixed}. Intersection omit-one coefficients,
their coefficient of variation, and intersection CV3 are missing because
they are not calculated. {opt rankcheck} still checks the rank after deleting
each intersection, which can reduce the time saved by {opt mixed}.{p_end}

{phang2}{cmd:summclust y x z, cluster(g h) nograph}{p_end}
{phang2}{cmd:summclust y x z, cluster(g h) mixed nograph}{p_end}

{title:Reading the two-way output}

{pstd}
{it:Regression Output} reports {cmd:CV1max} and either {cmd:CV3max} or,
with {opt mixed}, {cmd:CV31max}. Each row gives the focal coefficient,
standard error, t statistic, p value, and 95% confidence interval.{p_end}

{pstd}
Immediately below this table, the output identifies which component supplies
each max standard error: the conventional two-way estimator or the one-way
estimator for a named cluster variable. Ties are listed using the relative
tolerance in {cmd:r(max_tie_tolerance)}. With {opt mixed}, the two one-way
candidates are still CV3. Thus {cmd:CV3max} and {cmd:CV31max} can be identical
even when conventional CV3 and CV31 differ: the same one-way component may
be largest in both comparisons. The conventional focal variances and standard
errors are printed separately below the diagnostic summary.{p_end}

{pstd}
{it:Coefficients of Variation, G, and G*} identifies the actual cluster variables
and {cmd:intersect}. Its first four columns ({it:Ng}, {it:Leverage},
{it:Partial L.}, {it:beta no g}) are coefficients of variation displayed to four
decimals. G is an integer and Gstar is G*(0), displayed to two decimals.
The tables are designed for line widths of at least 80. Long cluster names
appear on a separate line when needed.{p_end}

{pstd}
Returned results retain full precision. {cmd:r(dimsummary)} has row names
G H I and column names
{cmd:cv_N cv_leverage cv_partlev cv_betanog clusters Gstar0}; their order maps to
{cmd:r(cluster_g)}, {cmd:r(cluster_h)}, and their intersection.{p_end}

{pstd}
{opt table} displays separate cluster-level diagnostic tables for G and H
when the corresponding dimension has at most 52 clusters, and returns both
tables regardless of size. It does not print the potentially much larger
intersection table. In all dimension-specific diagnostic matrices, row j
corresponds to the jth cluster in sorted nonmissing identifier order; for
intersections, sorting is by the (G,H) identifier pair. This rule applies to
numeric and string identifiers.{p_end}

{pstd}
{opt jackknife} is available only with one-way clustering. Specifying it with
two cluster variables gives error 198.{p_end}

{pstd}
Two-way estimation requires at least two nonempty clusters in each primary
dimension. Exactly two are allowed, but inference then uses one degree of
freedom and the command warns that it may be unreliable. Identical and nested
cluster dimensions are allowed. Without {opt mixed}, identical dimensions
reduce algebraically to one-way clustering. With {opt mixed}, they give
CV31 = 2 CV3_G - CV1_G instead.{p_end}

{title:Partial leverage}

{pstd}
For the focal regressor x, cluster g's partial leverage is{p_end}

{phang2}
sum(i in g) xtilde_i^2 / sum(i) xtilde_i^2,{p_end}

{pstd}
where xtilde is the Frisch-Waugh-Lovell residual from projecting x on every
other regressor, including the constant and fixed effects when present.
Partial leverage measures a cluster's share of the variation that identifies
the focal coefficient; the cluster shares sum to one.{p_end}

{title:Fixed effects and G*}

{pstd}
{opt fevar()} includes indicators for one or more categorical variables.
Both {opt fevar()} and {opt absorb()} accept numeric or string identifiers
with one-way or two-way clustering. String categories are grouped internally;
you do not need to create numeric codes or use {cmd:encode}. This also applies
to two-way models with {opt mixed}. Empty strings and numeric missing
values are excluded from the estimation sample.
Indicators included by {opt fevar()} enter the leverage and covariance calculations.
{opt absorb()} removes one fixed effect by subtracting group means. Each
absorbed group must lie wholly within a cluster; with two-way clustering it
must lie wholly within a G-H intersection. Use {opt fevar()} for fixed effects
that cross either clustering dimension. If the nesting condition fails,
the command stops with error 198.{p_end}

{pstd}
For string fixed effects, coefficient names use {it:fevar_j}, where j is the
category's position in the sorted nonempty levels in the estimation sample.
Levels outside that sample do not affect the names or omitted category.
If a name would exceed 32 characters, the variable-name prefix is shortened
and the suffix is {it:_fp_j}, where p is the position in {opt fevar()}.
For numeric fixed effects, one-way names use the numeric level and two-way
names use its sorted position.
All levels of the first fixed effect are included when there is no absorption;
the last level is omitted for subsequent fixed effects and for every fixed
effect used with {opt absorb()}.{p_end}

{pstd}
With one-way clustering, request G* using {opt gstar} or {opt rho()}.
If fixed effects occur at the cluster or subcluster level, only G*(0) is
available. Supplying {opt rho()} turns on the G* calculations even without
{opt gstar}.{p_end}

{pstd}
In two-way mode, the dimension summary always reports G*(0) separately for G,
H, and their intersections. {opt gstar} additionally returns G*(1), when
defined, by dimension; {opt rho()} likewise returns G*(rho). Availability is
checked separately for each dimension: G*(1) and G*(rho) are unavailable for
a dimension if fixed-effect groups are nested within its clusters. They
may remain available for the other dimensions. With {opt absorb()}, only
G*(0) is available in all three dimensions.{p_end}

{pstd}
Specify {opt rho()} to choose the value of rho used for G*(rho).
With {opt gstar} alone, the stored G*(rho) results use rho = 0.1234567.
This applies to both one-way and two-way clustering.{p_end}

{title:Deleted-cluster diagnostics}

{pstd}
With one-way clustering, an omit-one-cluster estimate is flagged if
any coefficient has absolute value at most 1e-8. This is not a mathematical
rank test: a coefficient that is truly zero can trigger it. The main CV3
result uses all omit-one estimates. When estimates are flagged, additional
results exclude those estimates.{p_end}

{pstd}
Specify {opt rankcheck} to return both the coefficient flags in
{cmd:r(legacyflag)} and actual numerical rank-loss indicators in
{cmd:r(rankfail)}. This option is diagnostic only and does not change warnings,
the estimates used in those additional results, or the reported inference.{p_end}

{pstd}
For two-way clustering, {opt rankcheck} checks G, H, and I separately and
returns {cmd:r(rankfail_g)}, {cmd:r(rankfail_h)}, and
{cmd:r(rankfail_intersection)}. An entry is 1 if deleting that cluster reduces
the rank relative to the full model, and 0 otherwise. Without {opt rankcheck},
these vectors are missing. Two-way rank warnings use this rank comparison,
rather than the one-way coefficient flag.{p_end}

{pstd}
Regardless of {opt rankcheck}, two-way estimation exits with error 498
if the focal coefficient is not numerically identified in the full model or
any delete-G or delete-H model. Without {opt mixed}, the same check applies
to every delete-intersection model. With {opt mixed}, those models are not
estimated. Rank loss confined to other coefficients is allowed. If the focal regressor
cannot be distinguished numerically from the other regressors, check for
collinearity and consider rescaling variables with very different magnitudes.
A missing required variance component also causes error 498.{p_end}

{title:Stored results: one-way clustering}

{pstd}
Let G be the number of clusters and K the number of explicit columns in the
estimation design. K includes the constant only when neither {opt fevar()} nor
{opt absorb()} is specified; absorbed indicators are not explicit columns.{p_end}

{p2colset 5 23 25 2}{...}
{p2col 5 23 25 2:Matrices}{p_end}
{synopt:{cmd:r(ng)}}G x 1 cluster sizes{p_end}
{synopt:{cmd:r(lever)}}G x 1 cluster leverages{p_end}
{synopt:{cmd:r(partlev)}}G x 1 partial leverages{p_end}
{synopt:{cmd:r(betanog)}}G x 1 focal-coefficient estimates after deleting each cluster{p_end}
{synopt:{cmd:r(cv1)}}K x K CV1 covariance matrix{p_end}
{synopt:{cmd:r(cv3)}}K x K CV3 covariance matrix{p_end}
{synopt:{cmd:r(cv3j)}}K x K CV3J covariance matrix; with {opt jackknife}{p_end}
{synopt:{cmd:r(legacyflag)}}G x 1 coefficient flags described above; with {opt rankcheck}{p_end}
{synopt:{cmd:r(rankfail)}}G x 1 deleted-design rank-loss flags; with {opt rankcheck}{p_end}
{synopt:{cmd:r(scall)}}G x 4 cluster table; with {opt table} when G <= 52{p_end}
{synopt:{cmd:r(cv3drop)}}CV3 covariance using unflagged deleted samples; with {opt regtable} and without {opt jackknife}, when coefficient flags occur{p_end}
{synopt:{cmd:r(cv3jdrop)}}CV3J covariance using unflagged deleted samples; available with {opt regtable jackknife} when coefficient flags occur{p_end}

{p2col 5 23 25 2:Numeric results stored as macros}{p_end}
{synopt:{cmd:r(beta)}}focal coefficient estimate{p_end}
{synopt:{cmd:r(cv1se)}}CV1 standard error{p_end}
{synopt:{cmd:r(cv1t)}}CV1 t statistic{p_end}
{synopt:{cmd:r(cv1p)}}CV1 two-sided p value using G-1 degrees of freedom{p_end}
{synopt:{cmd:r(cv1lci)}}CV1 95% confidence-interval lower endpoint{p_end}
{synopt:{cmd:r(cv1uci)}}CV1 95% confidence-interval upper endpoint{p_end}
{synopt:{cmd:r(cv3se)}}CV3 standard error{p_end}
{synopt:{cmd:r(cv3t)}}CV3 t statistic{p_end}
{synopt:{cmd:r(cv3p)}}CV3 two-sided p value using G-1 degrees of freedom{p_end}
{synopt:{cmd:r(cv3lci)}}CV3 95% confidence-interval lower endpoint{p_end}
{synopt:{cmd:r(cv3uci)}}CV3 95% confidence-interval upper endpoint{p_end}
{synopt:{cmd:r(cv3Jse)}}CV3J standard error; with {opt jackknife}{p_end}
{synopt:{cmd:r(cv3Jt)}}CV3J t statistic; with {opt jackknife}{p_end}
{synopt:{cmd:r(cv3Jp)}}CV3J two-sided p value; with {opt jackknife}{p_end}
{synopt:{cmd:r(cv3Jlci)}}CV3J 95% confidence-interval lower endpoint; with {opt jackknife}{p_end}
{synopt:{cmd:r(cv3Juci)}}CV3J 95% confidence-interval upper endpoint; with {opt jackknife}{p_end}
{synopt:{cmd:r(gstarzero)}}G*(0); with {opt gstar} or {opt rho()}{p_end}
{synopt:{cmd:r(gstarone)}}G*(1), when available{p_end}
{synopt:{cmd:r(gstarrho)}}G*(rho); with {opt gstar} or {opt rho()}, when available{p_end}
{p2colreset}{...}

{pstd}
The numeric results in the second group are stored as macros. To save one
as a scalar, expand its macro value, for example:{p_end}

{phang2}{cmd:. scalar focal_se = `r(cv3se)'}{p_end}

{pstd}
Use {cmd:return list} to inspect the available results. Save results you need
before running another command that replaces {cmd:r()}.{p_end}

{title:Stored results: two-way clustering}

{pstd}
K is the number of explicit regression columns. The suffixes {cmd:_g},
{cmd:_h}, and {cmd:_intersection} identify the two clustering dimensions and
their intersections. Unless stated otherwise, results below apply to both
ordinary and mixed estimation.{p_end}

{pstd}
With {opt mixed}, use {cmd:cv31} in place of {cmd:cv3} for the conventional
focal variance and standard error, max inference, and component selectors.
The {cmd:cv3} scalar and macro results for these quantities are not returned.
The {cmd:cv3} matrices listed below have the stated exceptions.{p_end}

{p2colset 5 34 36 2}{...}
{synopt:{cmd:r(b)}}1 x K coefficient vector{p_end}
{synopt:{cmd:r(cv1_conventional)}}K x K conventional CV1 matrix V_G+V_H-V_intersection{p_end}
{synopt:{cmd:r(cv3_conventional)}}K x K conventional CV3 matrix V_G+V_H-V_intersection; missing with {opt mixed}{p_end}
{synopt:{cmd:r(cv1)}, {cmd:r(cv3)}}aliases for the two conventional matrices{p_end}
{synopt:{cmd:r(cv1_g)}, {cmd:r(cv1_h)}, {cmd:r(cv1_intersection)}}K x K CV1 components{p_end}
{synopt:{cmd:r(cv3_g)}, {cmd:r(cv3_h)}, {cmd:r(cv3_intersection)}}K x K CV3 components; intersection component missing with {opt mixed}{p_end}
{synopt:{cmd:r(ng_g)}, {cmd:r(ng_h)}, {cmd:r(ng_intersection)}}dimension-specific cluster sizes{p_end}
{synopt:{cmd:r(lever_g)}, {cmd:r(lever_h)}, {cmd:r(lever_intersection)}}dimension-specific leverages{p_end}
{synopt:{cmd:r(partlev_g)}, {cmd:r(partlev_h)}, {cmd:r(partlev_intersection)}}dimension-specific partial leverages{p_end}
{synopt:{cmd:r(betanog_g)}, {cmd:r(betanog_h)}, {cmd:r(betanog_intersection)}}dimension-specific focal delete-cluster estimates; intersection estimates missing with {opt mixed}{p_end}
{synopt:{cmd:r(rankfail_g)}, {cmd:r(rankfail_h)}, {cmd:r(rankfail_intersection)}}dimension-specific rank-loss flags; with {opt rankcheck}, otherwise missing{p_end}
{synopt:{cmd:r(table_g)}, {cmd:r(table_h)}}cluster-level N, leverage, partial leverage, and delete-cluster coefficient; with {opt table}{p_end}
{synopt:{cmd:r(dimsummary)}}3 x 6 coefficients of variation, counts, and G*(0){p_end}
{synopt:{cmd:r(N)}}number of observations{p_end}
{synopt:{cmd:r(k)}, {cmd:r(rank)}}number of explicit regression columns and their numerical rank{p_end}
{synopt:{cmd:r(cluster_g)}, {cmd:r(cluster_h)}}names of the clustering variables{p_end}
{synopt:{cmd:r(N_clust_g)}, {cmd:r(N_clust_h)}, {cmd:r(N_clust_intersection)}}cluster counts{p_end}
{synopt:{cmd:r(G)}, {cmd:r(H)}, {cmd:r(I)}}short aliases for the cluster counts{p_end}
{synopt:{cmd:r(df_r)}}min(G,H)-1 inference degrees of freedom{p_end}
{synopt:{cmd:r(cv1_conventional_var)}, {cmd:r(cv3_conventional_var)}}focal conventional variances{p_end}
{synopt:{cmd:r(cv1_conventional_se)}, {cmd:r(cv3_conventional_se)}}focal conventional standard errors; missing when the variance is negative{p_end}
{synopt:{cmd:r(cv1_max_var)}, {cmd:r(cv3_max_var)}}focal max variances{p_end}
{synopt:{cmd:r(cv1_max_se)}, {cmd:r(cv3_max_se)}}focal max standard errors{p_end}
{synopt:{cmd:r(cv1max_component)}, {cmd:r(cv3max_component)}}selected component: two-way, G, or H{p_end}
{synopt:{cmd:r(cv1max_ties)}, {cmd:r(cv3max_ties)}}tie masks: add 1 for conventional, 2 for G, 4 for H; 7 means all three{p_end}
{synopt:{cmd:r(cv31)}, {cmd:r(cv31_conventional)}}K x K conventional CV31 covariance; with {opt mixed}{p_end}
{synopt:{cmd:r(cv31_conventional_var)}, {cmd:r(cv31_conventional_se)}}focal conventional CV31 variance and standard error; with {opt mixed}{p_end}
{synopt:{cmd:r(cv31_max_var)}, {cmd:r(cv31_max_se)}}focal max variance and standard error; with {opt mixed}{p_end}
{synopt:{cmd:r(cv31_max_t)}, {cmd:r(cv31_max_p)}}focal max t statistic and p value; numeric scalars with {opt mixed}{p_end}
{synopt:{cmd:r(cv31_max_lci)}, {cmd:r(cv31_max_uci)}}focal max 95% confidence limits; numeric scalars with {opt mixed}{p_end}
{synopt:{cmd:r(cv31max_component)}}with {opt mixed}, selected component: two-way (CV31), G (CV3), or H (CV3){p_end}
{synopt:{cmd:r(cv31max_ties)}}with {opt mixed}, tie mask with the same coding as above{p_end}
{synopt:{cmd:r(max_tie_tolerance)}}relative tie tolerance, 1e-10, using max(abs(selected variance),1e-300) as scale{p_end}
{synopt:{cmd:r(selected_estimator)}}{cmd:CV3} by default or {cmd:CV31} with {opt mixed}{p_end}
{synopt:{cmd:r(mixed)}}1 with {opt mixed}, 0 otherwise{p_end}
{synopt:{cmd:r(gstarzero_g)}, {cmd:r(gstarzero_h)}, {cmd:r(gstarzero_intersection)}}dimension-specific G*(0){p_end}
{synopt:{cmd:r(gstarone_g)}, {cmd:r(gstarone_h)}, {cmd:r(gstarone_intersection)}}dimension-specific G*(1), when available{p_end}
{synopt:{cmd:r(gstarrho_g)}, {cmd:r(gstarrho_h)}, {cmd:r(gstarrho_intersection)}}dimension-specific G*(rho); with {opt gstar} or {opt rho()}, when available{p_end}
{p2colreset}{...}

{pstd}
Conventional variances are retained when negative; their standard errors
are then missing. With {opt mixed}, use the {cmd:cv31} results for the
conventional covariance and max inference; the {cmd:cv3} and
{cmd:cv3_conventional} matrices are missing. The component selectors and
tie masks refer to the focal coefficient, not every coefficient in
{opt regtable}.{p_end}

{pstd}
Alternative names include {cmd:r(V_cv#_g)}, {cmd:r(V_cv#_h)},
and {cmd:r(V_cv#_i)}, where # is 1 or 3. Diagnostic names ending in {cmd:_i}
are aliases for the corresponding {cmd:_intersection} results.
{cmd:r(cv#var_max)} and {cmd:r(cv#se_max)} are aliases for the max variance
and standard error; # is 1 and either 3 or, with {opt mixed}, 31.{p_end}

{pstd}
The focal macros {cmd:r(beta)}, {cmd:r(cv1se)}, {cmd:r(cv1t)},
{cmd:r(cv1p)}, and the confidence limits are also returned in two-way mode.
The corresponding macros have the prefix {cmd:cv3} by default or
{cmd:cv31} with {opt mixed}. These standard errors and inference macros use
the max-se rule; the full covariance matrices are conventional two-way
matrices.{p_end}

{title:Examples}

{phang2}{cmd:. webuse nlswork, clear}{p_end}
{phang2}{cmd:. keep if inrange(age,20,40)}{p_end}
{phang2}{cmd:. summclust ln_wage msp union race, cluster(ind_code) fevar(grade age birth_yr)}{p_end}

{pstd}Use a sample restriction and suppress the graph:{p_end}
{phang2}{cmd:. summclust ln_wage msp union race, cluster(ind_code) fevar(grade age birth_yr) sample(south==1) nograph}{p_end}

{pstd}Absorb a fixed effect nested within clusters:{p_end}
{phang2}{cmd:. summclust ln_wage msp union race, cluster(ind_code) absorb(ind_code) jackknife nograph}{p_end}

{pstd}Request G*(rho), a full regression table, and rank diagnostics:{p_end}
{phang2}{cmd:. summclust ln_wage msp union race, cluster(ind_code) fevar(grade age birth_yr) rho(.5) jackknife regtable rankcheck nograph}{p_end}

{pstd}Compare conventional CV3 and mixed CV31 with clustering by age and industry:{p_end}
{phang2}{cmd:. webuse nlswork, clear}{p_end}
{phang2}{cmd:. keep if inrange(age,25,35)}{p_end}
{phang2}{cmd:. generate vismin = inrange(race,2,3)}{p_end}
{phang2}{cmd:. summclust hours vismin south, cluster(age ind_code) fevar(age birth_yr year ind_code) nograph}{p_end}
{phang2}{cmd:. summclust hours vismin south, cluster(age ind_code) fevar(age birth_yr year ind_code) mixed nograph}{p_end}

{pstd}
In this example, {cmd:CV3max} and {cmd:CV31max} both use the one-way CV3
component for {cmd:ind_code}. The conventional CV3 and CV31 standard errors differ.{p_end}

{pstd}For a dataset with outcome {cmd:y}, regressor {cmd:treatment}, controls, and firm and year identifiers:{p_end}
{phang2}{cmd:. summclust y treatment controls, cluster(firm year) rankcheck nograph}{p_end}

{pstd}Display separate firm- and year-level diagnostic tables:{p_end}
{phang2}{cmd:. summclust y treatment controls, cluster(firm year) table nograph}{p_end}

{pstd}Two-way fixed effects via explicit indicators:{p_end}
{phang2}{cmd:. summclust y treatment controls, cluster(firm year) fevar(firm year) rankcheck nograph}{p_end}

{pstd}Absorb intersection-level effects (validly nested in both dimensions):{p_end}
{phang2}{cmd:. egen firmyear = group(firm year)}{p_end}
{phang2}{cmd:. summclust y treatment controls, cluster(firm year) absorb(firmyear) nograph}{p_end}

{pstd}
If {cmd:industry} contains string categories such as "Manufacturing" and
"Retail", use it directly in {opt fevar()}. The cluster variables may also
be numeric or strings:{p_end}
{phang2}{cmd:. summclust y treatment controls, cluster(firm) fevar(industry) nograph}{p_end}
{phang2}{cmd:. summclust y treatment controls, cluster(firm year) fevar(industry) mixed nograph}{p_end}

{pstd}
A string firm identifier can also define an absorbed effect. With one-way
clustering by firm, each absorbed group is a cluster:{p_end}
{phang2}{cmd:. summclust y treatment controls, cluster(firm) absorb(firm) nograph}{p_end}

{pstd}
With two-way clustering, absorbing {cmd:firm} is only valid if each firm
belongs to one year in the estimation sample. For firms observed in several
years, use {opt fevar(firm)} to include firm fixed effects. If your model
calls for firm-year fixed effects, absorb the intersections as shown above.
Observations with missing model variables or empty string identifiers are
excluded from these examples.{p_end}

{title:Author}

{p 4}Matthew D. Webb{p_end}
{p 4}matt.webb@carleton.ca{p_end}

{title:Citation}

{p 4 8 2}
James G. MacKinnon, Morten O. Nielsen, and Matthew D. Webb. 2023.
"Leverage, Influence, and the Jackknife in Clustered Regression Models:
Reliable Inference Using summclust."{p_end}
