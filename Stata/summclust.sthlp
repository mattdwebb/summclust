{smcl}
{* *! version 4.004 March 2023}{...}
{help summclust:summclust}
{hline}{...}

{title:Title}

{pstd}
Cluster specific summary statistics{p_end}

{title:Syntax}

{phang}
{cmd:summclust} {it:varlist}, {it:cluster(varname)}  [  {it:options}]

{phang}
{it:varlist}the dependent variable, the independent variable of
interest, and other (binary or continuous) independent variables,

{phang}
{it:cluster}the clustering variable. 

{synoptset 45 tabbed}{...}
{synopthdr}
{synoptline}

{synopt:{opt fevar(varlist)}}creates fixed effects for the included
variables, using {cmd:i.varname}.{p_end}

{synopt:{opt absorb(varname)}}partials out this variable from the
regression before computing other statistics. This should only be used
for variables that are nested within the specified clusters. This
option is computationally faster than using {cmd: fevar}. {p_end}

{synopt:{opt jack:knife}}calculates the jackknife variance estimator
CV_3J in addition to CV_3.{p_end}

{synopt:{opt add:means}}displays additional summary statistics for
cluster variability based on alternative means.{p_end}

{synopt:{opt gstar}}calculates the effective number of clusters G*()
and G*(1).{p_end}

{synopt:{opt rho(scalar)}}calculates the effective number of
clusters, G*(rho), in addition to G*(0) and G*(1).  This option can be
used with out without the {cmd: gstar} option.  The value of rho must
be between 0 and 1.{p_end}

{synopt:{opt tab:le}}displays the cluster by cluster statistics.{p_end}

{synopt:{opt sam:ple}}allows for sample restrictions. For instance,
in order to restrict the analysis to individuals 25 years of age or
older based on a variable "age", sample(age>=25) should be added to
the list of options. {p_end}

{synopt:{opt nog:raph}}suppresses creation of the figure, which is
otherwise shown by default. {p_end}

{synopt:{opt reg:table}}displays the full regression table similar to
Stata's {cmd: regress} table, but with CV_3 standard errors. {p_end}

{title:Updates}

{marker description}{...}
{title:Description}

{pstd}{cmd:summclust} is a stand-alone command that summarizes
cluster variability and calculates a cluster jackknife variance
estimator. MacKinnon, Nielsen, and Webb (2023) documents it more fully
than this help file. The command calculates measures of cluster-level
influence and leverage. It can optionally calculate the effective
number of clusters. By default, it reports CV_1 and CV_3 standard
errors, and it can optionally report a CV_3J standard error. It
also, optionally, calculates additional measures of cluster-level
heterogeneity. By default, it produces a figure which can help
identify the source of cluster-level heterogeneity. Finally, it can
produce a full table of regression results.

{pstd}{cmd:summclust} by default calculates the CV_3 standard error.
With well-behaved samples, this should match the standard error
calculated using either Stata's native {cmd: jackknife: reg y x,
cluster(group)} or {cmd: reg y x, cluster(group) vce(jackknife)}
commands.  However, many samples are not well behaved. Specifically,
some of the omit-one-cluster subsamples may be singular. When they
are, {cmd:summclust} calculates two standard errors. One 
drops the singular subsamples, as the native Stata routines do. The
other uses a generalized inverse. {cmd:summclust} provides guidance as
to which standard error is likely to be more reliable. When {cmd:
regtable} is specified, and singular subsamples are present, two
versions of the regression table are displayed. Similarly, if 
{cmd:jackknife} is specified and there are singular subsamples, four
different standard errors are shown, either CV_3 or CV_3J, combined
with either the generalized inverse or one computed after
dropping the singular subsamples. 

{pstd}{cmd: nograph} suppresses creation of the figure, which is
otherwise shown by default. The figure shows four scatter plots:
leverage against observations per cluster, partial leverage against
observations per cluster, leverage against omit-one-cluster
coefficients, and partial leverage against omit-one-cluster
coefficients. This figure can be quite informative, but it is
computationally costly to produce. We therefore recommend invoking
this option after you have inspected the figure.

{pstd}
{cmd: regtable} when {cmd:jackknife} is specified, regtable uses the
CV_3J estimates to produce the regression table. Otherwise, CV_3
estimates are used.

{title:Stored results}

{pstd}
{cmd:summclust} stores the following in {cmd:r()}:

{p2col 5 20 24 2: Matrices}{p_end}
 
{synopt:{cmd:r(ng)}}The number of observations per cluster.{p_end}

{synopt:{cmd:r(lever)}}The cluster-specific leverage.{p_end}

{synopt:{cmd:r(part)}}The cluster-specific partial leverage.{p_end}

{synopt:{cmd:r(betanog)}}The estimate of beta when the g_th
cluster is omitted.{p_end}

{p2col 5 20 24 2: Scalars}{p_end}

{synopt:{cmd:r(gstarzero)}}The effective number of clusters for the 
coefficient of interest using rho=0.{p_end}

{synopt:{cmd:r(gstarrho)}}The effective number of clusters for the 
coefficient of interest using the scalar rho from the {cmd:rho}
option.{p_end}

{synopt:{cmd:r(gstarone)}}The effective number of clusters for the 
coefficient of interest using rho=1.{p_end}

{synopt:{cmd:r(beta)}}The estimated beta for the coefficient of
interest.{p_end}

{synopt:{cmd:r(cv1se)}}The CV_1 standard error for the coefficient of
interest.{p_end}

{synopt:{cmd:r(cv1t)}}The CV_1 t-statistic for the coefficient of
interest.{p_end}

{synopt:{cmd:r(cv1p)}}The P value for the null hypothesis that beta=0
for the coefficient of interest using the CV_1 standard error.{p_end}

{synopt:{cmd:r(cv1lci)}}The lower bound of the 95% confidence interval
for beta using the CV_1 standard error.{p_end}

{synopt:{cmd:r(cv1uci)}}The upper bound of the 95% confidence 
interval for beta using the CV_1 standard error.{p_end}

{synopt:{cmd:r(cv3se)}}The CV_3 standard error for the coefficient of
interest.{p_end}

{synopt:{cmd:r(cv3t)}}The CV_3 t-statistic for the coefficient of
interest.{p_end}

{synopt:{cmd:r(cv3p)}}The P value for the null hypothesis that beta=0
for the coefficient of interest using the CV_3 standard error.{p_end}

{synopt:{cmd:r(cv3lci)}}The lower bound of the 95% confidence interval
for beta using the CV_3 standard error.{p_end}

{synopt:{cmd:r(cv3uci)}}The upper bound of the 95% confidence  interval
for beta using the CV_3 standard error.{p_end}

{synopt:{cmd:r(cv3Jse)}}The CV_3J standard error for the coefficient
of interest.{p_end}

{synopt:{cmd:r(cv3Jt)}}The CV_3J t-statistic for the coefficient of
interest.{p_end}

{synopt:{cmd:r(cv3Jp)}}The P value for the null hypothesis that beta=0
for the coefficient of interest using the CV_3J standard error.{p_end}

{synopt:{cmd:r(cv3Jlci)}}The lower bound of the 95% confidence interval
for beta using the CV_3J standard error.{p_end}

{synopt:{cmd:r(cv3Juci)}}The upper bound of the 95% confidence 
interval for beta using the CV_3J standard error.{p_end}

{pstd}
{cmd:summclust} stores the following in {cmd:mata}:

{p2col 5 20 24 2: Matrices}{p_end}

{synopt:{cmd:cvstuff}}The matrix with the standard errors, t-statistics,
etc.{p_end}

{synopt:{cmd:clustsum}}The matrix with the measures of cluster
variability.{p_end}

{synopt:{cmd:bonus}}The matrix with additional measures of cluster
variability. Only calculated when the option {cmd:addmeans} is
specified.{p_end}

{synopt:{cmd:scall}}The matrix with the cluster-by-cluster statistics.
Only calculated when the option {cmd:table} is specified.{p_end}

{synopt:{cmd:cnames}}The string matrix with the cluster names, to
match with elements in scall. Only calculated when the option
{cmd:table} is specified.{p_end}

{synopt:{cmd:regresstab}}The matrix that is displayed when the
{cmd:regresstab} option is specified.{p_end}

{synopt:{cmd:regresstab}}The additional matrix that is displayed when
the {cmd:regresstab} option is specified and there are singular
clusters.{p_end}

{title:Examples}

{hline}

{pstd} nlswork -- using {cmd:regress}

{phang2}{cmd:. webuse nlswork, clear}

{phang2}{cmd:. keep if inrange(age,20,40)}

{phang2}{cmd:. reg ln_wage i.grade i.age i.birth_yr union race msp,
cluster(ind)}

{pstd} nlswork -- using {cmd:summclust}

{phang2}{cmd:. summclust ln_wage msp union race, fevar(grade age
birth_yr) cluster(ind) }

{pstd} adding industry fixed effects using {cmd:absorb}

{phang2}{cmd:. summclust ln_wage msp union race, fevar(grade age
birth_yr) absorb(ind) cluster(ind)}

{pstd} sample restrictions - using {cmd:sample}

{phang2}{cmd:. summclust ln_wage msp union race, fevar(grade age
birth_yr) sample(south==1) cluster(ind)}

{pstd} Effective Number of Clusters using {cmd:gstar} or {cmd:rho}.

{phang2}{cmd:. summclust ln_wage msp union race, fevar(grade age
birth_yr) cluster(ind) gstar}

{phang2}{cmd:. summclust ln_wage msp union race, fevar(grade age
birth_yr) cluster(ind) rho(0.5)}

{pstd} All Output.

{phang2}{cmd:. summclust ln_wage msp union race, fevar(grade age
birth_yr) absorb(ind) cluster(ind) table addmeans jack rho(0.5) regtable}


{title:Author}

{p 4}Matthew D. Webb{p_end}
{p 4}matt.webb@carleton.ca{p_end}

{title:Citation}

{p 4 8 2}{cmd:summclust} is not an official Stata command. It is a
free contribution to the research community.

Please cite:
{p 8 8 2} James G. MacKinnon,  Morten Ø. Nielsen, and Matthew D. Webb.
2023. Leverage, Influence, and the Jackknife in Clustered Regression
Models: Reliable Inference Using summclust.{p_end}
