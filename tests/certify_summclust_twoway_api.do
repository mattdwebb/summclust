version 13
clear all
set more off
capture log close _all
log using certify_summclust_twoway_api_results.log, replace text
run summclust.ado
run twowayjack.ado

display as text "=== TWO-WAY PUBLIC API CERTIFICATION ==="
set seed 9151
set obs 240
generate int firm = ceil(_n/24)
generate int year = mod(_n-1,6)+1
generate double x = rnormal()
generate double z = rnormal()
generate double y = .7*x-.2*z+rnormal()

* Default output and the preferred versus compatibility result names.
summclust y x z, cluster(firm year)
tempname apiinc1 apiinc3 apinative1 apinative3
matrix `apiinc1' = r(cv1_g)+r(cv1_h)-r(cv1_intersection)
matrix `apiinc3' = r(cv3_g)+r(cv3_h)-r(cv3_intersection)
assert mreldif(r(cv1),`apiinc1')<1e-14
assert mreldif(r(cv3),`apiinc3')<1e-14
scalar __api_cv1max = max(el(r(cv1),1,1),el(r(cv1_g),1,1),el(r(cv1_h),1,1))
scalar __api_cv3max = max(el(r(cv3),1,1),el(r(cv3_g),1,1),el(r(cv3_h),1,1))
assert reldif(r(cv1_max_var),scalar(__api_cv1max))<1e-14
assert reldif(r(cv3_max_var),scalar(__api_cv3max))<1e-14
assert reldif(r(cv1_max_se),sqrt(scalar(__api_cv1max)))<1e-14
assert reldif(r(cv3_max_se),sqrt(scalar(__api_cv3max)))<1e-14
local api_cv1_names : colnames r(cv1)
local api_cv3_names : rownames r(cv3)
assert "`api_cv1_names'"=="x _cons z"
assert "`api_cv3_names'"=="x _cons z"
assert r(N_clust_g)==10 & r(N_clust_h)==6 & r(N_clust_intersection)==60
assert r(G)==r(N_clust_g) & r(H)==r(N_clust_h) & r(I)==r(N_clust_intersection)
assert r(df_r)==min(r(G),r(H))-1
assert mreldif(r(cv1),r(cv1_conventional))<1e-15
assert mreldif(r(cv3),r(cv3_conventional))<1e-15
assert mreldif(r(cv1_g),r(V_cv1_g))<1e-15
assert mreldif(r(cv1_h),r(V_cv1_h))<1e-15
assert mreldif(r(cv1_intersection),r(V_cv1_i))<1e-15
assert mreldif(r(cv3_g),r(V_cv3_g))<1e-15
assert mreldif(r(cv3_h),r(V_cv3_h))<1e-15
assert mreldif(r(cv3_intersection),r(V_cv3_i))<1e-15
assert r(cv1_max_var)==r(cv1var_max) & r(cv1_max_se)==r(cv1se_max)
assert r(cv3_max_var)==r(cv3var_max) & r(cv3_max_se)==r(cv3se_max)
assert reldif(r(cv1_conventional_var),el(r(cv1_conventional),1,1))<1e-14
assert reldif(r(cv3_conventional_var),el(r(cv3_conventional),1,1))<1e-14
assert "`r(cluster_g)'"=="firm" & "`r(cluster_h)'"=="year"

* The compatibility command delegates and returns the native results.
matrix `apinative1' = r(cv1)
matrix `apinative3' = r(cv3)
quietly twowayjack y x z, cluster(firm year) nograph
assert mreldif(r(cv1),`apinative1')<1e-14
assert mreldif(r(cv3),`apinative3')<1e-14

* Option behavior: G*/rho, rank checks, tables, regression table, and graph policy.
summclust y x z, cluster(firm year) gstar rho(.5) rankcheck table regtable
confirm scalar r(gstarzero_g)
confirm scalar r(gstarzero_h)
confirm scalar r(gstarzero_intersection)
confirm scalar r(gstarone_g)
confirm scalar r(gstarrho_intersection)
confirm matrix r(rankfail_g)
confirm matrix r(rankfail_intersection)
confirm matrix r(table_g)
confirm matrix r(table_h)
assert rowsof(r(table_g))==r(G) & colsof(r(table_g))==4
assert rowsof(r(table_h))==r(H) & colsof(r(table_h))==4
summclust y x z, cluster(firm year) nograph

* Explicit and absorbed fixed effects on their valid two-way paths.
summclust y x z, cluster(firm year) fevar(firm year) rankcheck nograph
egen long firmyear = group(firm year)
summclust y x z, cluster(firm year) absorb(firmyear) nograph

* Identical dimensions reduce to one-way; nested dimensions remain supported.
summclust y x z, cluster(firm) nograph
local onecv1se "`r(cv1se)'"
local onecv3se "`r(cv3se)'"
summclust y x z, cluster(firm firm) nograph
assert reldif(`onecv1se',r(cv1_conventional_se))<1e-6
assert reldif(`onecv3se',r(cv3_conventional_se))<1e-6
generate int region = ceil(firm/2)
summclust y x z, cluster(region firm) nograph
assert r(G)==5 & r(H)==10 & r(I)==10 & r(df_r)==4

* String identifiers follow the same sorted-row contract.
generate str8 firms = "F"+string(firm)
generate str8 years = "Y"+string(year)
summclust y x z, cluster(firms years) table rankcheck nograph
assert rowsof(r(ng_g))==10 & rowsof(r(ng_h))==6 & rowsof(r(ng_intersection))==60
assert rowsof(r(table_g))==10 & rowsof(r(table_h))==6

* A deterministic negative conventional diagonal is preserved and flagged.
clear
set seed 17
set obs 80
generate int g = mod(_n-1,5)+1
generate int h = mod(floor((_n-1)/5),4)+1
generate double x = rnormal()
generate double z = rnormal()
generate double y = .4*x+.2*z+rnormal()
summclust y x z, cluster(g h) regtable
assert r(cv1_conventional_var)<0 & r(cv3_conventional_var)<0
assert missing(r(cv1_conventional_se)) & missing(r(cv3_conventional_se))
assert r(cv1_max_var)>0 & r(cv3_max_var)>0
assert r(cv1_max_se)<. & r(cv3_max_se)<.

* Too few clusters fail specifically; two-by-two is allowed with df=1.
replace g = 1
capture noisily summclust y x z, cluster(g h) nograph
assert _rc==498
replace g = mod(_n-1,2)+1
replace h = mod(floor((_n-1)/2),2)+1
summclust y x z, cluster(g h) nograph
assert r(df_r)==1

* Centered CV3J is one-way only; the two-way request is rejected deliberately.
capture noisily summclust y x z, cluster(g h) jackknife nograph
assert _rc==198
summclust y x z, cluster(g) jackknife nograph
confirm matrix r(cv3j)

display as result "TWO-WAY API CERTIFICATION PASSED"
log close
exit, clear
