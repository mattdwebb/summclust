version 13
clear all
set more off
set seed 23092026
set type float
capture log close _all
log using certify_integration_results.log, replace text
run summclust.ado

/* Exact pilot design: independent refits, including every FE deletion. */
foreach m in 1 25 {
    clear
    set obs 64
    generate int g=ceil(_n/8)
    generate int h=mod(_n-1,8)+1
    generate int cell=_n
    expand `m'
    generate double x=(1-2*(mod(g-1,4)>=2))*(1-2*(mod(h-1,4)>=2))
    generate double y=.123456789*x+.2*g-.3*h+rnormal()
    forvalues j=2/8 {
        generate byte rg`j'=g==`j'
        generate byte ch`j'=h==`j'
    }
    quietly regress y x rg* ch*
    scalar fullb=_b[x]
    foreach d in g h cell {
        quietly regress y x rg* ch*, vce(cluster `d')
        scalar ref1_`d'=_se[x]^2
        quietly summarize `d', meanonly
        local J=r(max)
        matrix refb_`d'=J(`J',1,.)
        scalar ref3_`d'=0
        forvalues j=1/`J' {
            quietly regress y x rg* ch* if `d'!=`j'
            matrix refb_`d'[`j',1]=_b[x]
            scalar ref3_`d'=ref3_`d'+(fullb-_b[x])^2*(`J'-1)/`J'
        }
    }
    quietly summclust y x rg* ch*, cluster(g h) nograph rankcheck
    assert abs(el(r(b),1,1)-fullb)<1e-12
    foreach pair in "g g" "h h" "cell intersection" {
        gettoken d tag : pair
        local tag=strtrim("`tag'")
        assert abs(el(r(cv1_`tag'),1,1)-ref1_`d')<1e-10
        assert abs(el(r(cv3_`tag'),1,1)-ref3_`d')<1e-10
        mata: st_numscalar("err",max(abs(st_matrix("r(betanog_`tag')")-st_matrix("refb_`d'"))))
        assert err<1e-10
    }
    mata: st_numscalar("err",max(abs(st_matrix("r(partlev_g)"):-1/8)))
    assert err<1e-12
    mata: st_numscalar("err",max(abs(st_matrix("r(partlev_h)"):-1/8)))
    assert err<1e-12
    quietly summclust y x rg* ch*, cluster(g) jackknife nograph
    assert !missing(el(r(cv3j),1,1))
    assert abs(el(r(cv1),1,1)-ref1_g)<1e-10
    assert abs(el(r(cv3),1,1)-ref3_g)<1e-10
    mata: st_numscalar("err",max(abs(st_matrix("r(partlev)"):-1/8)))
    assert err<1e-12
}

/* FE column names must identify the columns actually estimated. */
clear
set obs 240
generate int g=ceil(_n/24)
generate int h=mod(_n-1,6)+1
egen long cell=group(g h)
generate double x=rnormal()
generate double z=rnormal()
generate double y=.5*x+.2*z+g/3-h/5+rnormal()
forvalues j=1/10 {
    generate byte gd`j'=g==`j'
}
forvalues j=1/5 {
    generate byte hd`j'=h==`j'
}
quietly regress y x z gd* hd*, noconstant
matrix ref=e(b)
quietly summclust y x z, cluster(g h) fevar(g h) nograph
matrix actual=r(b)
local names : colnames actual
assert "`names'"=="x z g_1 g_2 g_3 g_4 g_5 g_6 g_7 g_8 g_9 g_10 h_1 h_2 h_3 h_4 h_5"
mata: st_numscalar("err",max(abs(st_matrix("actual")-st_matrix("ref"))))
assert err<1e-10

/* Temporary native Mata matrices cannot destroy the caller's objects. */
mata: TWX=42; TWY=43; TWX0=44; TWXC=45; TWG=46; TWH=47; TWI=48; TWFE=49; TWAB=50; TWABINFO=51
quietly summclust y x z, cluster(g g) nograph
assert r(cv1max_ties)==7 & r(cv3max_ties)==7
assert r(max_tie_tolerance)==1e-10
mata: assert((TWX,TWY,TWX0,TWXC,TWG,TWH,TWI,TWFE,TWAB,TWABINFO)==(42,43,44,45,46,47,48,49,50,51))
generate double shuffle=runiform()
sort shuffle
generate long incoming=_n
capture quietly summclust y x z, cluster(g h) absorb(g) nograph
assert _rc==198
assert incoming==_n
mata: assert(TWX==42 & TWY==43)

/* Nuisance-only rank loss is allowed; loss of focal identification is not. */
generate double rare=(g==1)*rnormal()
quietly summclust y x z rare, cluster(g h) rankcheck nograph
assert el(r(rankfail_g),1,1)==1
foreach option in "" "rankcheck" {
    capture quietly summclust y rare z, cluster(g h) `option' nograph
    assert _rc==498
    assert incoming==_n
}
generate double xcopy=x
capture quietly summclust y x xcopy, cluster(g h) nograph
assert _rc==498
assert incoming==_n
mata: assert(TWX==42 & TWY==43)

/* Absorbed degrees of freedom count surviving explicit columns only. */
generate double one=1
quietly regress y x z i.cell, vce(cluster g)
scalar refv=_se[x]^2
quietly summclust y x z one, cluster(g) absorb(cell) nograph
assert abs(el(r(cv1),1,1)-refv)<1e-10
quietly summclust y x z one, cluster(g h) absorb(cell) nograph
assert abs(el(r(cv1_g),1,1)-refv)<1e-10

/* Saturated one-way models must fail before dividing by residual df=0. */
clear
set obs 4
generate int g=_n
generate double x=(_n==1)
generate double z=(_n==2)
generate double w=(_n==3)
generate double y=rnormal()
capture quietly summclust y x z w, cluster(g) nograph
assert _rc==2001

display as result _newline "SUMMCLUST INTEGRATION CERTIFICATION PASSED"
log close
exit, clear
