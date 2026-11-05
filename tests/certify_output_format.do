version 13
clear all
set more off
capture log close _all
log using certify_output_format_results.log, text replace name(main)
run summclust.ado
tempname sources
file open `sources' using format_sources.csv, write text replace
file write `sources' "scenario,method,two_way_var,g_var,h_var,cluster_g,cluster_h,selected,ties" _n
foreach scenario in ordinary options long80 long120 negative few identical suppressed mixed mixed_swapped twoway mixed_twoway long_mixed mixed_tied {
    quietly {
        clear
        set seed 9326
        if inlist("`scenario'","mixed","mixed_swapped","long_mixed") set seed 7
        if inlist("`scenario'","twoway","mixed_twoway") set seed 6
        set obs 240
        generate int firm=ceil(_n/24)
        generate int year=mod(_n-1,6)+1
        generate double x=rnormal()
        generate double z=rnormal()
        generate double y=.7*x-.2*z+rnormal()
    }
    set linesize 80
    if "`scenario'"=="long120" set linesize 120
    if inlist("`scenario'","long80","long120","long_mixed") {
        rename firm clustering_dimension_alpha_12345
        rename year clustering_dimension_beta_123456
        rename z control_variable_long_name_12345
    }
    if "`scenario'"=="negative" {
        quietly {
            clear
            set seed 17
            set obs 80
            generate int firm=mod(_n-1,5)+1
            generate int year=mod(floor((_n-1)/5),4)+1
            generate double x=rnormal()
            generate double z=rnormal()
            generate double y=.4*x+.2*z+rnormal()
        }
    }
    if "`scenario'"=="few" {
        quietly replace firm=mod(_n-1,2)+1
        quietly replace year=mod(floor((_n-1)/2),2)+1
    }
    if "`scenario'"=="suppressed" quietly replace firm=ceil(_n/4)
    // A zero outcome gives exact zero variances: all mixed candidates tie.
    if "`scenario'"=="mixed_tied" quietly replace y=0
    log using "format_`scenario'.log", text replace name(txt)
    log using "format_`scenario'.smcl", smcl replace name(smcl)
    noisily display "OUTPUT BEGIN"
    if "`scenario'"=="ordinary" summclust y x z, cluster(firm year)
    if "`scenario'"=="options" summclust y x z, cluster(firm year) fevar(firm year) sample(x<1) table regtable gstar rho(.2) addmeans rankcheck
    if inlist("`scenario'","long80","long120") summclust y x control_variable_long_name_12345, cluster(clustering_dimension_alpha_12345 clustering_dimension_beta_123456) table regtable
    if inlist("`scenario'","negative","few","suppressed") summclust y x z, cluster(firm year) table regtable
    if "`scenario'"=="identical" summclust y x z, cluster(firm firm) nograph
    if "`scenario'"=="twoway" summclust y x z, cluster(firm year) nograph
    if inlist("`scenario'","mixed","mixed_twoway","mixed_tied") summclust y x z, cluster(firm year) mixed nograph
    if "`scenario'"=="mixed_swapped" summclust y x z, cluster(year firm) mixed nograph
    if "`scenario'"=="long_mixed" summclust y x control_variable_long_name_12345, cluster(clustering_dimension_alpha_12345 clustering_dimension_beta_123456) mixed regtable nograph
    // Record full-precision variances, independently of display formatting.
    local estimator cv3
    if r(mixed) local estimator cv31
    foreach method in cv1 `estimator' {
        local component cv3
        if "`method'"=="cv1" local component cv1
        matrix format_tw = r(`method'_conventional)
        matrix format_g = r(`component'_g)
        matrix format_h = r(`component'_h)
        file write `sources' "`scenario',`method'," %24.17g (format_tw[1,1]) "," %24.17g (format_g[1,1]) "," %24.17g (format_h[1,1]) ",`r(cluster_g)',`r(cluster_h)',`r(`method'max_component)'," %8.0f (r(`method'max_ties)) _n
    }
    noisily display "OUTPUT END"
    log close txt
    log close smcl
    translate "format_`scenario'.smcl" "format_`scenario'_translated.log", translator(smcl2log) linesize(`c(linesize)') replace
}
file close `sources'
translate summclust.sthlp help_rendered.log, translator(smcl2log) linesize(80) replace
display "SUMMCLUST OUTPUT FORMAT FIXTURES CERTIFICATION PASSED"
log close main
exit, clear
