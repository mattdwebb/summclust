version 13
clear all
set more off
set seed 23082026
capture log close _all
log using certify_summclust_results.log, text replace

run summclust.ado

/* Core specification and independent CV1/CV3/CV3J calculations. */
set obs 600
generate int g = ceil(_n/50)
generate double x = rnormal()
generate double z = .4*x + rnormal()
generate double w = rnormal()
bysort g: generate double ag = rnormal() if _n==1
bysort g: replace ag = ag[1]
generate double y = .7*x - .35*z + .2*w + ag + rnormal()

quietly regress y x z w, vce(cluster g)
tempname bref V1ref
matrix `bref' = (_b[x] \ _b[_cons] \ _b[z] \ _b[w])
matrix `V1ref' = e(V)
mata: __vr=st_matrix("`V1ref'"); st_matrix("`V1ref'",__vr[(1,4,2,3),(1,4,2,3)])

quietly summclust y x z w, cluster(g) jackknife gstar rankcheck regtable nograph
tempname V1 V3 V3J NG PL BN RF LF RT
matrix `V1' = r(cv1)
matrix `V3' = r(cv3)
matrix `V3J' = r(cv3j)
matrix `NG' = r(ng)
matrix `PL' = r(partlev)
matrix `BN' = r(betanog)
matrix `RF' = r(rankfail)
matrix `LF' = r(legacyflag)
matrix `RT' = regresstab
local beta = real("`r(beta)'")
local cv1se = real("`r(cv1se)'")
local cv3se = real("`r(cv3se)'")
local cv3jse = real("`r(cv3Jse)'")

/* Legacy numeric macros use display precision; matrices retain full precision. */
assert abs(`beta' - el(`bref',1,1)) < 1e-6
assert abs(`cv1se' - sqrt(el(`V1',1,1))) < 1e-6
assert abs(`cv3se' - sqrt(el(`V3',1,1))) < 1e-6
assert abs(`cv3jse' - sqrt(el(`V3J',1,1))) < 1e-6
mata: st_numscalar("__dV1",max(abs(st_matrix("`V1'"):-st_matrix("`V1ref'"))))
display as text "CV1 max difference from regress: " as result %12.4e scalar(__dV1)
assert scalar(__dV1) < 1e-8
assert rowsof(`V1')==4 & colsof(`V1')==4
assert rowsof(`V3')==4 & colsof(`V3')==4
assert rowsof(`V3J')==4 & colsof(`V3J')==4
local cvnames : rownames `V1'
assert "`cvnames'" == "x _cons z w"
local rtnames : rownames `RT'
assert "`rtnames'" == "x _cons z w"
forvalues j=1/4 {
    assert abs(el(`RT',`j',1)-el(`bref',`j',1)) < 1e-8
}
mata: st_numscalar("__symj",max(abs(st_matrix("`V3J'"):-st_matrix("`V3J'")')))
mata: st_numscalar("__missj",sum(missing(st_matrix("`V3J'"))))
assert scalar(__symj) < 1e-12
assert scalar(__missj)==0
mata: st_numscalar("__nrf",sum(st_matrix("`RF'")))
assert scalar(__nrf)==0

tempname B V3ref V3Jref
matrix `B' = J(4,12,.)
forvalues j=1/12 {
    quietly regress y x z w if g!=`j'
    matrix `B'[1,`j'] = _b[x]
    matrix `B'[2,`j'] = _b[_cons]
    matrix `B'[3,`j'] = _b[z]
    matrix `B'[4,`j'] = _b[w]
}
mata: __B=st_matrix("`B'"); __b=st_matrix("`bref'"); __G=cols(__B)
mata: __D=__b*J(1,__G,1)-__B
mata: st_matrix("`V3ref'",((__G-1)/__G)*(__D*__D'))
mata: __bm=mean(__B'); __DJ=__B'-J(__G,1,1)*__bm
mata: st_matrix("`V3Jref'",((__G-1)/__G)*quadcross(__DJ,__DJ))
mata: st_numscalar("__dV3",max(abs(st_matrix("`V3'"):-st_matrix("`V3ref'"))))
mata: st_numscalar("__dV3J",max(abs(st_matrix("`V3J'"):-st_matrix("`V3Jref'"))))
assert scalar(__dV3) < 1e-8
assert scalar(__dV3J) < 1e-8

/* Independent Frisch-Waugh-Lovell partial leverage. */
quietly regress x z w
tempvar xr xr2 gnum gden p tag
quietly predict double `xr', residuals
quietly generate double `xr2' = `xr'^2
quietly bysort g: egen double `gnum' = total(`xr2')
quietly egen double `gden' = total(`xr2')
quietly generate double `p' = `gnum'/`gden'
quietly egen byte `tag' = tag(g)
quietly sort g
tempname PLref
quietly mkmat `p' if `tag', matrix(`PLref')
mata: st_numscalar("__dPL",max(abs(st_matrix("`PL'"):-st_matrix("`PLref'"))))
assert scalar(__dPL) < 1e-10
mata: st_numscalar("__sumNG",sum(st_matrix("`NG'")))
assert scalar(__sumNG)==600

/* Every valid rho value must be distinguishable from an omitted rho(). */
quietly summclust y x z w, cluster(g) rho(.12345) nograph
assert "`r(gstarzero)'" != ""
assert "`r(gstarrho)'" != ""

/* The optimized command must not destroy similarly named user objects. */
matrix __summclust_cv1 = (42)
quietly summclust y x z w, cluster(g) nograph
assert el(__summclust_cv1,1,1)==42

/* Numeric fixed effects with a string cluster identifier and G*. */
generate int nested = (g-1)*2 + mod(_n,2) + 1
generate str8 gs = "c" + string(g)
capture quietly summclust y x z w, cluster(gs) fevar(nested) gstar nograph
assert _rc==0
assert "`r(gstarzero)'" != ""
assert "`r(gstarone)'" == ""

/* Absorption is invariant to incoming row order and matches areg beta. */
generate int absorbed = (g-1)*5 + floor(mod(_n-1,50)/10) + 1
generate double shuffle = runiform()
sort shuffle
generate long incoming_order = _n
quietly areg y x z w, absorb(absorbed)
local bareg = _b[x]
quietly summclust y x z w, cluster(g) absorb(absorbed) jackknife nograph
assert abs(real("`r(beta)'")-`bareg') < 1e-6
assert incoming_order==_n

/* sample(), missing values, alternate focal regressor, and storage types. */
generate byte useit = mod(_n,4)!=0
replace y = . if mod(_n,97)==0
recast float x
quietly regress y z x w if useit==1
local bz = _b[z]
quietly summclust y z x w, cluster(g) sample(useit==1) nograph
assert abs(real("`r(beta)'")-`bz') < 1e-6

/* Dropped-subsample CV3J table uses the retained-cluster degrees of freedom. */
clear
set seed 23082031
set obs 1200
generate int g = ceil(_n/100)
generate double x = rnormal()
generate double z = cond(g==1,rnormal(),0)
generate double w = rnormal()
generate double y = .6*x + z + .25*w + rnormal()
quietly summclust y x z w, cluster(g) rankcheck jackknife regtable nograph
tempname LFD VJD RTD
matrix `LFD' = r(legacyflag)
matrix `VJD' = r(cv3jdrop)
matrix `RTD' = regresstabdrop
mata: st_numscalar("__nflagd",sum(st_matrix("`LFD'")))
local gkeep = 12-scalar(__nflagd)
assert scalar(__nflagd)==1
assert abs(el(`RTD',1,2)-sqrt(el(`VJD',1,1))) < 1e-10
local pdrop = 2*min(ttail(`gkeep'-1,el(`RTD',1,3)),1-ttail(`gkeep'-1,el(`RTD',1,3)))
assert abs(el(`RTD',1,4)-`pdrop') < 1e-10
local cidrop = invttail(`gkeep'-1,.025)*el(`RTD',1,2)
assert abs(el(`RTD',1,5)-(el(`RTD',1,1)-`cidrop')) < 1e-10
assert abs(el(`RTD',1,6)-(el(`RTD',1,1)+`cidrop')) < 1e-10

/* Exact deleted-design rank failure. */
clear
set seed 23082027
set obs 1200
generate int g = ceil(_n/100)
generate double x1 = rnormal()
generate double noise = rnormal()
generate double x2 = x1 + cond(g==1,.1*noise,0)
generate double z = rnormal()
generate double y = .5*x1-.2*x2+.1*z+rnormal()
quietly summclust y x1 x2 z, cluster(g) rankcheck jackknife nograph
tempname RF2
matrix `RF2' = r(rankfail)
assert el(`RF2',1,1)==1

/* A zero coefficient can trigger the compatibility flag without rank loss. */
clear
set seed 23082028
set obs 800
generate int g = ceil(_n/100)
generate double x = rnormal()
generate double z = rnormal()
generate double y = 1 + 2*x
quietly summclust y x z, cluster(g) rankcheck nograph
tempname LF2 RF3
matrix `LF2' = r(legacyflag)
matrix `RF3' = r(rankfail)
mata: st_numscalar("__nlf2",sum(st_matrix("`LF2'")))
mata: st_numscalar("__nrf3",sum(st_matrix("`RF3'")))
assert scalar(__nlf2)>0
assert scalar(__nrf3)==0

/* Graph path and documented error codes. */
clear
set seed 23082029
set obs 400
generate int g = ceil(_n/40)
generate double x = rnormal()
generate double z = rnormal()
generate double y = .5*x+rnormal()
capture graph drop _all
capture noisily summclust y x z, cluster(g)
assert _rc==0
capture graph describe Graph
assert _rc==0
capture graph drop _all

capture quietly summclust y x z, cluster(g) rho(1.5) nograph
assert _rc==198
generate int badabs = mod(_n-1,7)+1
capture quietly summclust y x z, cluster(g) absorb(badabs) nograph
assert _rc==198
capture quietly summclust y x z, cluster(g) sample(0) nograph
assert _rc==2000
replace g = 1
capture quietly summclust y x z, cluster(g) nograph
assert _rc==498

display as result _newline "SUMMCLUST 5.0.0 CERTIFICATION PASSED"
log close
