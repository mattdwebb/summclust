{smcl}
{* *! version 3.2.3 1 June 2021}{...}
{help summclust:summclust}
{hline}{...}

{title:Title}

{pstd}
Cluster specific summary statistics{p_end}

{title:Syntax}

{phang}
{cmd:summclust} {it:depvar}, {it:yvar(varname)} {it:xvar(varlist)}  {it:cluster(varname)}  [  {it:options}]

{phang}
 {it:depvar} is the dependendent variable of interest,

{phang}
 {it:xvar} are other, binary or continuous, dependendent variables of interest,

{phang}
 {it:yvar} is the independendent variable of interest, 

{phang}
 {it:cluster} is the clustering variable. 




{synoptset 45 tabbed}{...}
{synopthdr}
{synoptline}

{synopt:{opt fevar(varlist)}} creates fixed effects for the included varaibles, using {cmd:i.varname}.{p_end}

{synopt:{opt absorb(varname)}} creates fixed effects for the included varaibles, using {cmd:i.varname}.  This should
only be used for variables that are nested within the specified clusters. {p_end}

{synopt:{opt jack:knife}} calculates the jackknife variance estimator cv_3J in addition to cv_3.{p_end}


{synopt:{opt sv:ars}} displays additional summary statistics  for cluster variability based on alternative means.{p_end}

{synopt:{opt gstar}} calcutaes the effective number of clusters G*() and G*(1).{p_end}

{synopt:{opt rho(scalar)}} calcutaes the effective number of clusters, G*(rho), in addition to G*() and G*(1).  This option can be used with out without the {cmd: gstar} option.  The value of rho must be between 0 and 1.{p_end}

{synopt:{opt tab:le}} displays the cluster by cluster statistics.{p_end}

{synopt:{opt sam:ple}} allows for sample restrictions.  For instance, in order to restrict the
analysis to individuals 25 years of age or older based on a variable "age", sample(age>=25)
should be added to the list of options. {p_end}


{title:Updates}

{marker description}{...}
{title:Description}

{pstd}
{cmd:summclust} is a stand alone command that summarizes cluster variability and calculates a cluster jackknife variance estimator. MacKinnon, Nielsen, and Webb (2022) documents it more fully than this help file.  The command 
will calculate measures of cluster level influence and leverage. It will optionally calculate the effective number of clusters. By default it estimates  CV_1 and CV_3 standard errors, it will optionally estimate a CV_3J standard error. Finally, it will optionally calculate additional measures of cluster level heterogeneity.



{title:Stored results}

{pstd}
{cmd:summclust} stores the following in {cmd:r()}:

{p2col 5 20 24 2: Matrices}{p_end}
 

{synopt:{cmd:r(ng)}}The number of observations per cluster.{p_end}

{synopt:{cmd:r(lever)}}The cluster specific leverage.{p_end}

{synopt:{cmd:r(part)}}The cluster specific partial leverage.{p_end}

{synopt:{cmd:r(betanog)}}The estimate of beta when omitting the g_th cluster.{p_end}

{p2col 5 20 24 2: Scalars}{p_end}


{synopt:{cmd:r(gstarzero)}}The effective number of clusters for the 
coefficient of interest using rho=0.{p_end}

{synopt:{cmd:r(gstarrho)}}The effective number of clusters for the 
coefficient of interest using the scalar rho from the {cmd:rho} option.{p_end}

{synopt:{cmd:r(gstarone)}}The effective number of clusters for the 
coefficient of interest using rho=1.{p_end}


{synopt:{cmd:r(beta)}}The estimated beta for the coefficient of interest.{p_end}


{synopt:{cmd:r(cv1se)}}The CV_1 standard error for the 
coefficient of interest.{p_end}


{synopt:{cmd:r(cv1t)}}The CV_1 t-statistic for the coefficient
of interest.{p_end}


{synopt:{cmd:r(cv1p)}}The P value for the null hypothesis that 
beta = 0 for the coefficient of interest using the CV_1
standard error.{p_end}

{synopt:{cmd:r(cv1lci)}}The lower bound of the 95% confidence interval
for beta using the CV_1 standard error.{p_end}


{synopt:{cmd:r(cv1uci)}}The upper bound of the 95% confidence  interval
for beta using the CV_1 standard error.{p_end}

{synopt:{cmd:r(cv3se)}}The CV_3 standard error for the 
coefficient of interest.{p_end}


{synopt:{cmd:r(cv3t)}}The CV_3 t-statistic for the coefficient
of interest.{p_end}


{synopt:{cmd:r(cv3p)}}The P value for the null hypothesis that 
beta = 0 for the coefficient of interest using the CV_3
standard error.{p_end}

{synopt:{cmd:r(cv3lci)}}The lower bound of the 95% confidence interval
for beta using the CV_3 standard error.{p_end}


{synopt:{cmd:r(cv3uci)}}The upper bound of the 95% confidence  interval
for beta using the CV_3 standard error.{p_end}


{synopt:{cmd:r(cv3Jse)}}The CV_3J standard error for the 
coefficient of interest.{p_end}


{synopt:{cmd:r(cv3Jt)}}The CV_3J t-statistic for the coefficient
of interest.{p_end}


{synopt:{cmd:r(cv3Jp)}}The P value for the null hypothesis that 
beta = 0 for the coefficient of interest using the CV_3J
standard error.{p_end}

{synopt:{cmd:r(cv3Jlci)}}The lower bound of the 95% confidence interval
for beta using the CV_3J standard error.{p_end}


{synopt:{cmd:r(cv3Juci)}}The upper bound of the 95% confidence  interval
for beta using the CV_3J standard error.{p_end}

{pstd}
{cmd:summclust} stores the following in {cmd:mata}:

{p2col 5 20 24 2: Matrices}{p_end}

{synopt:{cmd:cvstuff}}The matrix with the standard errors, t-statistics, etc.{p_end}

{synopt:{cmd:clustsum}}The matrix with the measures of cluster variability.{p_end}

{synopt:{cmd:bonus}}The matrix with additional measures of cluster variability.  
Only calculated when the option {cmd:svars} is specified.{p_end}

{synopt:{cmd:scall}}The matrix with the cluster by cluster statistics.  
Only calculated when the option {cmd:table} is specified.{p_end}


{title:Examples}

{hline}

{pstd} nlswork - using {cmd:regress}


{phang2}{cmd:. webuse nlswork, clear}


{phang2}{cmd:. reg ln_wage  i.grade i.age  i.birth_yr union race msp, cluster(ind)}

{pstd} nlswork - using {cmd:summclust}

{phang2}{cmd:. summclust msp, yvar(ln_wage) xvar(union race) fevar(grade age birth_yr) cluster(ind)}

{pstd}  adding industry fixed effects using {cmd:absorb}

{phang2}{cmd:. summclust msp, yvar(ln_wage) xvar(union race) fevar(grade age birth_yr) absorb(ind) cluster(ind)}

{pstd} sample restrictions - using {cmd:sample}

{phang2}{cmd:. summclust msp, yvar(ln_wage) xvar(union race i.grade) fevar(age birth_yr) sample(age>=25) cluster(ind)}


{pstd} Effective Number of Clusters using {cmd:gstar} or {cmd:rho}.

{phang2}{cmd:. summclust msp, yvar(ln_wage) xvar(union race) fevar(grade age birth_yr) cluster(ind) gstar}

{phang2}{cmd:. summclust msp, yvar(ln_wage) xvar(union race) fevar(grade age birth_yr) cluster(ind) rho(0.5)}

{pstd} All Output.

{phang2}{cmd:. summclust msp, yvar(ln_wage) xvar(union race) fevar(grade age birth_yr) absorb(dubind) cluster(ind) table svars jack rho(0.5)}



{title:Author}

{p 4}Matthew D. Webb{p_end}
{p 4}matt.webb@carleton.ca{p_end}


{title:Citation}

{p 4 8 2}{cmd:summclust} is not an official Stata command. It is a free contribution to the research community.
Please cite: {p_end}
{p 8 8 2} James G. MacKinnon,  Morten Ø. Nielsen, and Matthew D. Webb. 2022. Leverage, Influence, and the Jackknife in Clustered Regression Models: Reliable Inference Using summclust.{p_end}
