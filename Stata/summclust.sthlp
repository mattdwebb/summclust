{smcl}
{* *! version 5.1.1 23sep2026}{...}
{vieweralsosee "summclust" "help summclust"}{...}
{hline}

{title:Title}

{pstd}
{cmd:summclust} {hline 2} One- and two-way cluster-robust variance estimation,
leverage, influence, and cluster-jackknife diagnostics

{title:Syntax}

{phang}
{cmd:summclust} {it:depvar focalvar controls},
{opt cluster(clustvar)} [{it:options}]

{phang}
{cmd:summclust} {it:depvar focalvar controls},
{opt cluster(clustvar1 clustvar2)} [{it:options}]

{pstd}
The varlist must contain at least three variables. The first is the dependent
variable, the second is the coefficient of interest, and all remaining
variables are controls. Factor-variable notation is not accepted in this
varlist; use {opt fevar()} for categorical fixed effects.

{pstd}
The command retains its historical {opt sample()} option. Standard Stata
{cmd:if}, {cmd:in}, and weight syntax are not supported.

{synoptset 28 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{opt cluster(varlist)}}one or two numeric or string cluster identifiers; required{p_end}
{synopt:{opt fevar(varlist)}}include indicators for numeric categorical variables{p_end}
{synopt:{opt absorb(varname)}}partial out one fixed effect; see the nesting conditions below{p_end}
{synopt:{opt jack:knife}}also calculate one-way CV3J; not available with two-way clustering{p_end}
{synopt:{opt gstar}}calculate effective-cluster measures G*(0) and, when available, G*(1){p_end}
{synopt:{opt rho(#)}}also calculate G*(rho); 0.00001 <= # <= 1{p_end}
{synopt:{opt sam:ple(exp)}}restrict the estimation sample using a logical expression{p_end}
{synopt:{opt tab:le}}display cluster-by-cluster diagnostics; in two-way mode, display separate G and H tables{p_end}
{synopt:{opt add:means}}display additional measures based on alternative means{p_end}
{synopt:{opt reg:table}}display results for every explicit regression coefficient{p_end}
{synopt:{opt nog:raph}}suppress the one-way four-panel graph; two-way graphs are suppressed by policy{p_end}
{synopt:{opt rank:check}}return an actual deleted-design rank diagnostic in addition to the legacy flag{p_end}
{synoptline}

{title:Description}

{pstd}
{cmd:summclust} reports cluster sizes, leverage, partial leverage,
omit-one-cluster coefficient estimates, CV1 and CV3 inference, and optional
CV3J and effective-cluster calculations. By default it draws four scatter
plots combining cluster size, leverage, partial leverage, and the
omit-one-cluster coefficient for one-way clustering. With two cluster
variables it instead reports a concise focal-coefficient table and separate
G, H, and nonempty-intersection diagnostics.

{pstd}
In two-way output, G is the first variable named in {opt cluster()}, H is the
second, and I denotes their nonempty intersections. The output header repeats
the actual variable names and reports the inference degrees of freedom as
min(G,H)-1.

{pstd}
Version 5.1.1 uses an optimized Mata implementation while retaining
{cmd:version 13}. Coefficients, CV1, CV3, CV3J standard errors, leverage,
omit-one estimates, and G* agree numerically with version 4.22 in validated
well-posed specifications, up to floating-point tolerance.
The dependent-variable copy now preserves double precision. Absorbed-model
CV1 degrees of freedom use the rank after demeaning; these corrections can
also change results relative to older versions. This release was tested on
StataNow/SE 18.5; the declared Stata 13 minimum has not been independently tested.

{title:Native two-way clustering}

{pstd}
Specify two variables in {opt cluster()}, for example
{cmd:cluster(firm year)}. The estimation sample and design are constructed
once. The command computes separate one-way CV1 and CV3 component matrices for
firm (G), year (H), and their nonempty intersections (I), then returns the
conventional inclusion-exclusion matrices

{phang2}
V^(3) = V_G + V_H - V_I.

{pstd}
Each CV1 component uses its own finite-sample correction based on G, H, or I.
The conventional two-way matrix can be indefinite in finite samples and is
returned without eigenvalue repair in {cmd:r(cv1_conventional)} or
{cmd:r(cv3_conventional)}; {cmd:r(cv1)} and {cmd:r(cv3)} are retained aliases.
If the focal diagonal is negative, the command preserves and reports that
variance, returns a missing conventional standard error, warns explicitly,
and continues to report the max-se inference.

{pstd}
For the focal coefficient, displayed inference uses the max-se rule of
MacKinnon, Nielsen, and Webb:

{phang2}
V_max = max{V^(3)[1,1], V_G[1,1], V_H[1,1]}.

{pstd}
The t distribution has min(G,H)-1 degrees of freedom. Max-se is a
single-restriction rule. The command deliberately does not construct or return
a full "max covariance matrix". {opt regtable} displays coefficient-wise max
standard errors, which likewise must not be used as a joint covariance matrix.

{pstd}
{opt table} displays separate cluster-level diagnostic tables for G and H
when the corresponding dimension has at most 52 clusters, and returns both
tables regardless of size. It does not print the potentially much larger
intersection table. In all dimension-specific diagnostic matrices, row j
corresponds to the jth cluster in sorted nonmissing identifier order; for
intersections, sorting is by the (G,H) identifier pair. This rule applies to
numeric and string identifiers.

{pstd}
The default two-way CV3 components use deviations of each delete-cluster
estimate from the full-sample estimate, as defined in the paper. The
alternative centered CV3J requested by {opt jackknife} is not studied there;
two-way use of that option exits with error 198 rather than exposing an
unvalidated estimator. One-way {opt jackknife} remains supported.

{pstd}
Two-way estimation requires at least two nonempty clusters in each primary
dimension. Exactly two are allowed, but inference then uses one degree of
freedom and the command warns that it may be unreliable. Identical cluster
variables are allowed and reduce algebraically to one-way clustering;
nested, nonidentical cluster dimensions are also allowed.

{title:Partial leverage}

{pstd}
For the focal regressor x, version 5.0 defines cluster g's partial leverage as

{phang2}
sum(i in g) xtilde_i^2 / sum(i) xtilde_i^2,

{pstd}
where xtilde is the Frisch-Waugh-Lovell residual from projecting x on every
other column of the estimation design. This corrects a version-4.22 macro
indirection bug. Consequently, {cmd:r(partlev)}, the partial-leverage column in
tables, and graph panels using partial leverage may intentionally differ from
version 4.22. Coefficients and covariance estimators are unaffected by this
correction.

{title:Fixed effects and G*}

{pstd}
{opt fevar()} explicitly adds indicator columns and therefore contributes to
leverage and full covariance calculations. {opt absorb()} uses within-group
demeaning and is intended only when every absorbed group belongs to one
cluster for one-way clustering. For two-way clustering, every absorbed group
must be nested in both primary dimensions, and therefore in one intersection.
This stronger rule is needed because one shared within-transformed design must
remain valid for delete-G, delete-H, and delete-I calculations. Use
{opt fevar()} for fixed effects that cross either dimension. The command exits
with error 198 if the applicable nesting condition fails.

{pstd}
In two-way returned coefficient names, {it:fevar_j} denotes the jth sorted
nonmissing level in the estimation sample, not necessarily the numeric value j.
All levels of the first fixed effect are included when there is no absorption;
the last level is omitted for subsequent fixed effects and for every fixed
effect used with {opt absorb()}.

{pstd}
When fixed effects occur at the cluster or subcluster level, G*(1) and G*(rho)
are not available; G*(0) is still reported. Supplying {opt rho()} turns on the
G* calculations even without {opt gstar}.

{pstd}
In two-way mode, the dimension summary always reports G*(0) separately for G,
H, and their intersections. {opt gstar} additionally returns G*(1), when
defined, by dimension; {opt rho()} likewise returns G*(rho).

{title:Deleted-cluster diagnostics}

{pstd}
For one-way backward compatibility, the default warning and dropped-subsample policy
uses the historical flag: a deleted-cluster coefficient vector is flagged if
any coefficient has absolute value at most 1e-8. This is not a mathematical
rank test; a legitimately zero coefficient can trigger it. The primary CV3
calculation continues to use all deleted-cluster estimates, with a generalized
inverse when needed.

{pstd}
Specify {opt rankcheck} to return both the historical flags in
{cmd:r(legacyflag)} and actual numerical rank-loss indicators in
{cmd:r(rankfail)}. This option is diagnostic only and does not change warnings,
dropped-subsample calculations, or displayed recommendations.

{pstd}
For two-way clustering, {opt rankcheck} is evaluated independently for G, H,
and I and returns {cmd:r(rankfail_g)}, {cmd:r(rankfail_h)}, and
{cmd:r(rankfail_i)}. Two-way warnings use these true rank checks; the historical
near-zero-coefficient proxy is not used by the native two-way path.

{pstd}
Regardless of {opt rankcheck}, native two-way estimation exits with error 498
if the focal coefficient is not numerically identified in the full model or
any delete-G, delete-H, or delete-intersection model. Rank loss confined to
nuisance coefficients is allowed. The numerical check requires the focal row
of invsym(X'X)*(X'X) to agree with its identity row within 1e-7. Very ill-conditioned
designs can therefore be rejected; rescale or reconsider the specification.
Missing variance components also cause error 498 rather than partial inference.

{title:Stored results}

{pstd}
Let G be the number of clusters and K the number of explicit columns in the
estimation design. K includes the constant only when neither {opt fevar()} nor
{opt absorb()} is specified; absorbed indicators are not explicit columns.

{p2colset 5 23 25 2}{...}
{p2col 5 23 25 2:Matrices}{p_end}
{synopt:{cmd:r(ng)}}G x 1 cluster sizes{p_end}
{synopt:{cmd:r(lever)}}G x 1 cluster leverages{p_end}
{synopt:{cmd:r(partlev)}}G x 1 partial leverages{p_end}
{synopt:{cmd:r(betanog)}}G x 1 focal-coefficient estimates after deleting each cluster{p_end}
{synopt:{cmd:r(cv1)}}K x K CV1 covariance matrix{p_end}
{synopt:{cmd:r(cv3)}}K x K CV3 covariance matrix{p_end}
{synopt:{cmd:r(cv3j)}}K x K CV3J covariance matrix; with {opt jackknife}{p_end}
{synopt:{cmd:r(legacyflag)}}G x 1 historical flags; with {opt rankcheck}{p_end}
{synopt:{cmd:r(rankfail)}}G x 1 deleted-design rank-loss flags; with {opt rankcheck}{p_end}
{synopt:{cmd:r(scall)}}G x 4 cluster table; with {opt table} when G <= 52{p_end}
{synopt:{cmd:r(cv3drop)}}CV3 covariance using unflagged deleted samples; available with {opt regtable} when historical flags occur{p_end}
{synopt:{cmd:r(cv3jdrop)}}CV3J covariance using unflagged deleted samples; available with {opt regtable jackknife} when historical flags occur{p_end}

{p2col 5 23 25 2:Numeric macros (legacy interface)}{p_end}
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
{synopt:{cmd:r(gstarrho)}}G*(rho), when {opt rho()} is supplied and the measure is available{p_end}
{p2colreset}{...}

{pstd}
Numeric results in the second group remain returned macros, rather than
returned scalars, to preserve the version-4.22 interface. New code should use
the documented {cmd:r()} results. Legacy named Stata matrices and Mata working
objects used by earlier releases remain available where the corresponding
option is requested, but they are not recommended as a stable programming
interface.

{title:Additional stored results with two-way clustering}

{p2colset 5 34 36 2}{...}
{synopt:{cmd:r(b)}}1 x K coefficient vector{p_end}
{synopt:{cmd:r(cv1_conventional)}}K x K conventional CV1 matrix V_G+V_H-V_intersection{p_end}
{synopt:{cmd:r(cv3_conventional)}}K x K conventional CV3 matrix V_G+V_H-V_intersection{p_end}
{synopt:{cmd:r(cv1)}, {cmd:r(cv3)}}aliases for the two conventional matrices{p_end}
{synopt:{cmd:r(cv1_g)}, {cmd:r(cv1_h)}, {cmd:r(cv1_intersection)}}K x K CV1 components{p_end}
{synopt:{cmd:r(cv3_g)}, {cmd:r(cv3_h)}, {cmd:r(cv3_intersection)}}K x K CV3 components{p_end}
{synopt:{cmd:r(ng_g)}, {cmd:r(ng_h)}, {cmd:r(ng_intersection)}}dimension-specific cluster sizes{p_end}
{synopt:{cmd:r(lever_g)}, {cmd:r(lever_h)}, {cmd:r(lever_intersection)}}dimension-specific leverages{p_end}
{synopt:{cmd:r(partlev_g)}, {cmd:r(partlev_h)}, {cmd:r(partlev_intersection)}}dimension-specific partial leverages{p_end}
{synopt:{cmd:r(betanog_g)}, {cmd:r(betanog_h)}, {cmd:r(betanog_intersection)}}dimension-specific focal delete-cluster estimates{p_end}
{synopt:{cmd:r(rankfail_g)}, {cmd:r(rankfail_h)}, {cmd:r(rankfail_intersection)}}dimension-specific rank-loss flags{p_end}
{synopt:{cmd:r(table_g)}, {cmd:r(table_h)}}cluster-level N, leverage, partial leverage, and delete-cluster coefficient; with {opt table}{p_end}
{synopt:{cmd:r(dimsummary)}}3 x 6 coefficients of variation, counts, and G*(0){p_end}
{synopt:{cmd:r(N_clust_g)}, {cmd:r(N_clust_h)}, {cmd:r(N_clust_intersection)}}cluster counts{p_end}
{synopt:{cmd:r(G)}, {cmd:r(H)}, {cmd:r(I)}}short aliases for the cluster counts{p_end}
{synopt:{cmd:r(df_r)}}min(G,H)-1 inference degrees of freedom{p_end}
{synopt:{cmd:r(cv1_conventional_var)}, {cmd:r(cv3_conventional_var)}}focal conventional variances{p_end}
{synopt:{cmd:r(cv1_conventional_se)}, {cmd:r(cv3_conventional_se)}}focal conventional standard errors; missing when the variance is negative{p_end}
{synopt:{cmd:r(cv1_max_var)}, {cmd:r(cv3_max_var)}}focal MNW max variances{p_end}
{synopt:{cmd:r(cv1_max_se)}, {cmd:r(cv3_max_se)}}focal MNW max standard errors{p_end}
{synopt:{cmd:r(cv1max_component)}, {cmd:r(cv3max_component)}}selected component: two-way, G, or H{p_end}
{synopt:{cmd:r(cv1max_ties)}, {cmd:r(cv3max_ties)}}tie masks: add 1 for conventional, 2 for G, 4 for H; 7 means all three{p_end}
{synopt:{cmd:r(max_tie_tolerance)}}relative tie tolerance, 1e-10, using max(abs(selected variance),1e-300) as scale{p_end}
{synopt:{cmd:r(gstarzero_g)}, {cmd:r(gstarzero_h)}, {cmd:r(gstarzero_intersection)}}dimension-specific G*(0){p_end}
{synopt:{cmd:r(gstarone_g)}, {cmd:r(gstarone_h)}, {cmd:r(gstarone_intersection)}}dimension-specific G*(1), when available{p_end}
{synopt:{cmd:r(gstarrho_g)}, {cmd:r(gstarrho_h)}, {cmd:r(gstarrho_intersection)}}dimension-specific G*(rho), with {opt rho()} when available{p_end}
{p2colreset}{...}

{pstd}
For compatibility with the pre-release 5.1 test interface,
{cmd:r(V_cv#_g)}, {cmd:r(V_cv#_h)}, {cmd:r(V_cv#_i)}, diagnostic names ending
in {cmd:_i}, {cmd:r(cv#var_max)}, and {cmd:r(cv#se_max)} remain aliases for the
corresponding results above.

{pstd}
The legacy focal macros {cmd:r(beta)}, {cmd:r(cv1se)}, {cmd:r(cv1t)},
{cmd:r(cv1p)}, confidence limits, and their CV3 analogues are also returned in
two-way mode. In that mode the standard errors and inference macros refer to
the max-se rule, while the full {cmd:r(cv1)} and {cmd:r(cv3)} matrices are the
conventional inclusion-exclusion matrices.

{title:Examples}

{phang2}{cmd:. webuse nlswork, clear}
{phang2}{cmd:. keep if inrange(age,20,40)}
{phang2}{cmd:. summclust ln_wage msp union race, cluster(ind) fevar(grade age birth_yr)}

{pstd}Use a sample restriction and suppress the graph:{p_end}
{phang2}{cmd:. summclust ln_wage msp union race, cluster(ind) fevar(grade age birth_yr) sample(south==1) nograph}

{pstd}Absorb a fixed effect nested within clusters:{p_end}
{phang2}{cmd:. summclust ln_wage msp union race, cluster(ind) absorb(ind) jackknife nograph}

{pstd}Request G*(rho), a full regression table, and rank diagnostics:{p_end}
{phang2}{cmd:. summclust ln_wage msp union race, cluster(ind) fevar(grade age birth_yr) rho(.5) jackknife regtable rankcheck nograph}

{pstd}Native two-way clustering by firm and year:{p_end}
{phang2}{cmd:. summclust y treatment controls, cluster(firm year) rankcheck nograph}

{pstd}Display separate firm- and year-level diagnostic tables:{p_end}
{phang2}{cmd:. summclust y treatment controls, cluster(firm year) table nograph}

{pstd}Two-way fixed effects via explicit indicators:{p_end}
{phang2}{cmd:. summclust y treatment controls, cluster(firm year) fevar(firm year) rankcheck nograph}

{pstd}Absorb intersection-level effects (validly nested in both dimensions):{p_end}
{phang2}{cmd:. egen firmyear = group(firm year)}
{phang2}{cmd:. summclust y treatment controls, cluster(firm year) absorb(firmyear) nograph}

{title:Author}

{p 4}Matthew D. Webb{p_end}
{p 4}matt.webb@carleton.ca{p_end}

{title:Citation}

{p 4 8 2}
James G. MacKinnon, Morten O. Nielsen, and Matthew D. Webb. 2023.
"Leverage, Influence, and the Jackknife in Clustered Regression Models:
Reliable Inference Using summclust."{p_end}
