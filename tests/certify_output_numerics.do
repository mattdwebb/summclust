version 13
clear all
set more off
capture log close _all
log using certify_output_numerics_results.log, text replace
run summclust.ado
do output_audit_helpers.do
set seed 30926
set obs 240
generate int firm=ceil(_n/24)
generate int year=mod(_n-1,6)+1
generate double x=rnormal()+.1*firm
generate double z=rnormal()
generate double y=.7*x-.2*z+firm/10-year/5+rnormal()
audit_model y x z, cluster(firm year) id(numeric_plain) out(.)
audit_model y x z, cluster(firm year) fevar(firm year) id(numeric_fe) out(.)

/* Independent table moments, max components and unrounded inference.
   Repeated fits may change summation order at machine precision. */
quietly summclust y x z, cluster(firm year) nograph
tempname original1 original3 summary b
matrix `original1'=r(cv1)
matrix `original3'=r(cv3)
matrix `summary'=r(dimsummary)
matrix `b'=r(b)
quietly summclust y x z, cluster(firm year) table regtable addmeans nograph
assert mreldif(r(cv1),`original1')<1e-14
assert mreldif(r(cv3),`original3')<1e-14
assert mreldif(r(dimsummary),`summary')<1e-14
assert mreldif(r(b),`b')<1e-14
local row=0
foreach d in g h intersection {
    local ++row
    local column=0
    foreach stat in ng lever partlev betanog {
        local ++column
        mata: __q=st_matrix("r(`stat'_`d')"); st_numscalar("__cv",sqrt(variance(__q))/abs(mean(__q)))
        assert abs(__cv-el(r(dimsummary),`row',`column'))<1e-12
    }
}
foreach method in cv1 cv3 {
    assert r(`method'_max_var)==max(el(r(`method'_conventional),1,1),el(r(`method'_g),1,1),el(r(`method'_h),1,1))
    assert r(`method'_max_se)==sqrt(r(`method'_max_var))
    assert abs(real("`r(`method't)'")-el(r(b),1,1)/r(`method'_max_se))<1e-12
    assert abs(real("`r(`method'p)'")-2*ttail(r(df_r),abs(real("`r(`method't)'"))))<1e-12
    assert abs(real("`r(`method'lci)'")-(el(r(b),1,1)-invttail(r(df_r),.025)*r(`method'_max_se)))<1e-12
    assert abs(real("`r(`method'uci)'")-(el(r(b),1,1)+invttail(r(df_r),.025)*r(`method'_max_se)))<1e-12
}
display "SUMMCLUST OUTPUT NUMERICS CERTIFICATION PASSED"
log close
exit, clear
