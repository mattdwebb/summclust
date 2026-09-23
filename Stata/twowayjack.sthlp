{smcl}
{* *! version 2.0.0 23aug2026}{...}
{title:Title}

{pstd}
{cmd:twowayjack} {hline 2} Compatibility front end for native two-way
{cmd:summclust}

{title:Syntax}

{phang}
{cmd:twowayjack} {it:depvar focalvar controls} {ifin},
{opt cluster(var1 var2)} [{it:summclust options}]

{title:Description}

{pstd}
Version 2 is a thin r-class delegate. It marks {cmd:if}/{cmd:in}, combines that
restriction with optional {opt sample()}, and calls
{cmd:summclust ..., cluster(var1 var2)}. It makes no separate G, H, or
intersection command calls and returns the native two-way {cmd:summclust}
results unchanged.

{pstd}
See {help summclust} for estimator definitions, stored results, fixed-effect
requirements, the max-se rule, and the reason two-way {opt jackknife} is not
available.
