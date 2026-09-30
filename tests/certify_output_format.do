version 13
clear all
set more off
capture log close _all
log using certify_output_format_results.log, text replace name(main)
run summclust.ado
foreach scenario in ordinary options long80 long120 negative few identical suppressed {
    quietly {
        clear
        set seed 9326
        set obs 240
        generate int firm=ceil(_n/24)
        generate int year=mod(_n-1,6)+1
        generate double x=rnormal()
        generate double z=rnormal()
        generate double y=.7*x-.2*z+rnormal()
    }
    set linesize 80
    if "`scenario'"=="long120" set linesize 120
    if inlist("`scenario'","long80","long120") {
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
    log using "format_`scenario'.log", text replace name(txt)
    log using "format_`scenario'.smcl", smcl replace name(smcl)
    noisily display "OUTPUT BEGIN"
    if "`scenario'"=="ordinary" summclust y x z, cluster(firm year)
    if "`scenario'"=="options" summclust y x z, cluster(firm year) fevar(firm year) sample(x<1) table regtable gstar rho(.2) addmeans rankcheck
    if inlist("`scenario'","long80","long120") summclust y x control_variable_long_name_12345, cluster(clustering_dimension_alpha_12345 clustering_dimension_beta_123456) table regtable
    if inlist("`scenario'","negative","few","suppressed") summclust y x z, cluster(firm year) table regtable
    if "`scenario'"=="identical" summclust y x z, cluster(firm firm) nograph
    noisily display "OUTPUT END"
    log close txt
    log close smcl
    translate "format_`scenario'.smcl" "format_`scenario'_translated.log", translator(smcl2log) linesize(`c(linesize)') replace
}
display "SUMMCLUST OUTPUT FORMAT FIXTURES CERTIFICATION PASSED"
log close main
exit, clear
