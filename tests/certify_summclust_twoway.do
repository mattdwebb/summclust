version 13
clear all
set more off
set seed 5102026
capture log close _all
log using certify_summclust_twoway_results.log, replace text
run summclust.ado

capture program drop _tw_matclose
program define _tw_matclose
    args A B tol
    mata: st_numscalar("__tw_diff",max(abs(st_matrix("`A'"):-st_matrix("`B'"))))
    assert scalar(__tw_diff) < `tol'
end

/* 1. Balanced design: direct CV1 and explicit delete-cluster CV3 references. */
clear
set obs 180
generate int g = ceil(_n/30)
generate int h = mod(_n-1,5)+1
generate double x = rnormal()+.15*g-.08*h
generate double z = .25*x+rnormal()
bysort g: generate double ag = rnormal() if _n==1
bysort g: replace ag = ag[1]
bysort h: generate double ah = rnormal() if _n==1
bysort h: replace ah = ah[1]
generate double y = .8*x-.3*z+ag+ah+rnormal()
egen long i = group(g h)

tempname b0 V1g0 V1h0 V1i0 V3g0 V3h0 V3i0 BG BH BI raw
quietly regress y x z
matrix `b0' = (_b[x] \ _b[_cons] \ _b[z])
quietly regress y x z, vce(cluster g)
matrix `raw' = e(V)
mata: __v=st_matrix("`raw'"); st_matrix("`V1g0'",__v[(1,3,2),(1,3,2)])
quietly regress y x z, vce(cluster h)
matrix `raw' = e(V)
mata: __v=st_matrix("`raw'"); st_matrix("`V1h0'",__v[(1,3,2),(1,3,2)])
quietly regress y x z, vce(cluster i)
matrix `raw' = e(V)
mata: __v=st_matrix("`raw'"); st_matrix("`V1i0'",__v[(1,3,2),(1,3,2)])
matrix `BG' = J(3,6,.)
forvalues j=1/6 {
    quietly regress y x z if g!=`j'
    matrix `BG'[1,`j'] = _b[x]
    matrix `BG'[2,`j'] = _b[_cons]
    matrix `BG'[3,`j'] = _b[z]
}
matrix `BH' = J(3,5,.)
forvalues j=1/5 {
    quietly regress y x z if h!=`j'
    matrix `BH'[1,`j'] = _b[x]
    matrix `BH'[2,`j'] = _b[_cons]
    matrix `BH'[3,`j'] = _b[z]
}
matrix `BI' = J(3,30,.)
forvalues j=1/30 {
    quietly regress y x z if i!=`j'
    matrix `BI'[1,`j'] = _b[x]
    matrix `BI'[2,`j'] = _b[_cons]
    matrix `BI'[3,`j'] = _b[z]
}
mata: __B=st_matrix("`BG'"); __b=st_matrix("`b0'"); __D=__b*J(1,6,1)-__B; st_matrix("`V3g0'",(5/6)*(__D*__D'))
mata: __B=st_matrix("`BH'"); __b=st_matrix("`b0'"); __D=__b*J(1,5,1)-__B; st_matrix("`V3h0'",(4/5)*(__D*__D'))
mata: __B=st_matrix("`BI'"); __b=st_matrix("`b0'"); __D=__b*J(1,30,1)-__B; st_matrix("`V3i0'",(29/30)*(__D*__D'))

generate double shuffle = runiform()
sort shuffle
generate long incoming = _n
quietly summclust y x z, cluster(g h) nograph rankcheck gstar rho(.2) regtable
assert incoming==_n
assert r(G)==6 & r(H)==5 & r(I)==30 & r(df_r)==4
assert rowsof(r(ng_g))==6 & rowsof(r(ng_h))==5 & rowsof(r(ng_i))==30
assert rowsof(r(dimsummary))==3 & colsof(r(dimsummary))==6
assert r(gstarzero_g)>0 & r(gstarzero_g)<=6
assert r(gstarzero_h)>0 & r(gstarzero_h)<=5
assert r(gstarzero_i)>0 & r(gstarzero_i)<=30
tempname A1g A1h A1i A3g A3h A3i A1 A3 C1 C3
matrix `A1g'=r(V_cv1_g)
matrix `A1h'=r(V_cv1_h)
matrix `A1i'=r(V_cv1_i)
matrix `A3g'=r(V_cv3_g)
matrix `A3h'=r(V_cv3_h)
matrix `A3i'=r(V_cv3_i)
matrix `A1'=r(cv1)
matrix `A3'=r(cv3)
_tw_matclose `A1g' `V1g0' 1e-8
_tw_matclose `A1h' `V1h0' 1e-8
_tw_matclose `A1i' `V1i0' 1e-8
_tw_matclose `A3g' `V3g0' 1e-8
_tw_matclose `A3h' `V3h0' 1e-8
_tw_matclose `A3i' `V3i0' 1e-8
matrix `C1'=`V1g0'+`V1h0'-`V1i0'
matrix `C3'=`V3g0'+`V3h0'-`V3i0'
_tw_matclose `A1' `C1' 1e-8
_tw_matclose `A3' `C3' 1e-8
scalar mx1=max(el(`C1',1,1),el(`V1g0',1,1),el(`V1h0',1,1))
scalar mx3=max(el(`C3',1,1),el(`V3g0',1,1),el(`V3h0',1,1))
assert abs(r(cv1var_max)-mx1)<1e-10
assert abs(r(cv3var_max)-mx3)<1e-10
capture mata: TWX
assert _rc!=0
capture quietly summclust y x z, cluster(g h) jackknife nograph
assert _rc==198

/* 2. Coincident dimensions reduce exactly to one-way clustering. */
quietly summclust y x z, cluster(g) nograph
tempname OG1 OG3
matrix `OG1'=r(cv1)
matrix `OG3'=r(cv3)
quietly summclust y x z, cluster(g g) nograph
assert r(G)==r(H) & r(H)==r(I)
tempname TG1 TG3
matrix `TG1'=r(cv1)
matrix `TG3'=r(cv3)
_tw_matclose `TG1' `OG1' 1e-8
_tw_matclose `TG3' `OG3' 1e-8

/* 3. Nested dimensions reduce to the coarser component. */
generate int coarse = ceil(g/2)
quietly summclust y x z, cluster(coarse) nograph
tempname OC1 OC3
matrix `OC1'=r(cv1)
matrix `OC3'=r(cv3)
quietly summclust y x z, cluster(coarse g) nograph
assert r(G)==3 & r(H)==6 & r(I)==6
tempname TC1 TC3
matrix `TC1'=r(cv1)
matrix `TC3'=r(cv3)
_tw_matclose `TC1' `OC1' 1e-8
_tw_matclose `TC3' `OC3' 1e-8

/* 4. String identifiers, missing identifiers, and sample() restrictions. */
generate str8 gs = "G"+string(g)
generate str8 hs = "H"+string(h)
quietly summclust y x z, cluster(gs hs) nograph
tempname S1 S3
matrix `S1'=r(cv1)
matrix `S3'=r(cv3)
quietly summclust y x z, cluster(g h) nograph
tempname N1 N3
matrix `N1'=r(cv1)
matrix `N3'=r(cv3)
_tw_matclose `S1' `N1' 1e-8
_tw_matclose `S3' `N3' 1e-8
generate byte keep = mod(incoming,4)!=0
replace gs = "" if mod(incoming,17)==0
quietly count if keep & gs!=""
local expectedN=r(N)
quietly summclust y x z, cluster(gs hs) sample(keep==1) nograph
assert r(N)==`expectedN'
assert incoming==_n

/* 5. Unbalanced layout with empty intersections. */
clear
set obs 63
generate int g = ceil(_n/7)
generate int h = mod(_n-1,7)+1
drop if mod(g+2*h,4)==0
generate long cell = _n
expand 1+mod(cell,4)
generate double x = rnormal()+.1*g
generate double z = rnormal()
generate double y = .6*x-.2*z+rnormal()
egen long inter = group(g h)
quietly summarize g, meanonly
local UG=r(max)
quietly summarize h, meanonly
local UH=r(max)
quietly summarize inter, meanonly
local UI=r(max)
assert `UI' < `UG'*`UH'
foreach d in g h inter {
    quietly summclust y x z, cluster(`d') nograph
    tempname Q1`d' Q3`d'
    matrix `Q1`d''=r(cv1)
    matrix `Q3`d''=r(cv3)
}
quietly summclust y x z, cluster(g h) nograph
assert r(G)==`UG' & r(H)==`UH' & r(I)==`UI'
tempname U1g U1h U1i U3g U3h U3i
matrix `U1g'=r(V_cv1_g)
matrix `U1h'=r(V_cv1_h)
matrix `U1i'=r(V_cv1_i)
matrix `U3g'=r(V_cv3_g)
matrix `U3h'=r(V_cv3_h)
matrix `U3i'=r(V_cv3_i)
_tw_matclose `U1g' `Q1g' 1e-8
_tw_matclose `U1h' `Q1h' 1e-8
_tw_matclose `U1i' `Q1inter' 1e-8
_tw_matclose `U3g' `Q3g' 1e-8
_tw_matclose `U3h' `Q3h' 1e-8
_tw_matclose `U3i' `Q3inter' 1e-8

/* 6. Two-way fixed effects via generalized inverse, including rank flags. */
clear
set obs 192
generate int g = ceil(_n/24)
generate int h = mod(_n-1,6)+1
generate double x = rnormal()+.2*g-.1*h
generate double z = rnormal()
generate double y = .5*x+.25*z+g/4-h/5+rnormal()
egen long inter = group(g h)
foreach d in g h inter {
    quietly summclust y x z, cluster(`d') fevar(g h) rankcheck nograph
    tempname F1`d' F3`d'
    matrix `F1`d''=r(cv1)
    matrix `F3`d''=r(cv3)
}
quietly summclust y x z, cluster(g h) fevar(g h) rankcheck nograph
assert el(r(rankfail_g),1,1)==1
mata: st_numscalar("__rfg",sum(st_matrix("r(rankfail_g)"))); st_numscalar("__rfh",sum(st_matrix("r(rankfail_h)"))); st_numscalar("__rfi",sum(st_matrix("r(rankfail_i)")))
assert scalar(__rfg)==8 & scalar(__rfh)==6 & scalar(__rfi)==0
assert abs(el(r(V_cv1_g),1,1)-el(`F1g',1,1))<1e-8
assert abs(el(r(V_cv1_h),1,1)-el(`F1h',1,1))<1e-8
assert abs(el(r(V_cv1_i),1,1)-el(`F1inter',1,1))<1e-8
assert abs(el(r(V_cv3_g),1,1)-el(`F3g',1,1))<1e-8
assert abs(el(r(V_cv3_h),1,1)-el(`F3h',1,1))<1e-8
assert abs(el(r(V_cv3_i),1,1)-el(`F3inter',1,1))<1e-8

/* 7. Valid intersection absorb; invalid crossed absorb fails without reordering. */
foreach d in g h inter {
    quietly summclust y x z, cluster(`d') absorb(inter) nograph
    tempname W1`d' W3`d'
    matrix `W1`d''=r(cv1)
    matrix `W3`d''=r(cv3)
}
quietly summclust y x z, cluster(g h) absorb(inter) nograph
assert abs(el(r(V_cv1_g),1,1)-el(`W1g',1,1))<1e-8
assert abs(el(r(V_cv1_h),1,1)-el(`W1h',1,1))<1e-8
assert abs(el(r(V_cv1_i),1,1)-el(`W1inter',1,1))<1e-8
assert abs(el(r(V_cv3_g),1,1)-el(`W3g',1,1))<1e-8
assert abs(el(r(V_cv3_h),1,1)-el(`W3h',1,1))<1e-8
assert abs(el(r(V_cv3_i),1,1)-el(`W3inter',1,1))<1e-8
generate double shuf=runiform()
sort shuf
generate long incoming=_n
capture quietly summclust y x z, cluster(g h) absorb(g) nograph
assert _rc==198
assert incoming==_n

/* 8. True rank deficiency and a highly collinear design. */
generate double rare = (g==1)*rnormal()
quietly summclust y x z rare, cluster(g h) rankcheck nograph
mata: st_numscalar("__rare",sum(st_matrix("r(rankfail_g)")))
assert scalar(__rare)>=1
generate double znear = x + 1e-8*rnormal()
capture quietly summclust y x znear, cluster(g h) rankcheck nograph
assert _rc==498
replace znear = x + .1*rnormal()
quietly summclust y x znear, cluster(g h) rankcheck nograph
assert !missing(r(cv1var_max)) & !missing(r(cv3var_max))

/* 9. Negative conventional focal variances remain exposed; max-se stays valid. */
clear
set seed 17
set obs 80
generate int g = mod(_n-1,5)+1
generate int h = mod(floor((_n-1)/5),4)+1
generate double x = rnormal()
generate double z = rnormal()
generate double y = .4*x+.2*z+rnormal()
quietly summclust y x z, cluster(g h) nograph
assert el(r(cv1),1,1)<0 & el(r(cv3),1,1)<0
assert r(cv1var_max)>0 & r(cv3var_max)>0
assert real("`r(cv1se)'")>0 & real("`r(cv3se)'")>0
assert inlist("`r(cv1max_component)'","G","H")
assert inlist("`r(cv3max_component)'","G","H")

display as result _newline "SUMMCLUST 5.1.1 TWO-WAY CERTIFICATION PASSED"
log close
exit, clear
