version 18
clear all
set more off
capture log close _all
capture mkdir output
log using output/stata_mixed_tests.log, text replace
run Stata/summclust.ado
mata:
void mx_write(string scalar spec) {
    real scalar f,j,k
    real matrix a
    string rowvector vars
    string scalar name,rl,cl
    string matrix rn,cn
    f=fopen("output/stata_mixed.csv","a")
    vars=tokens(st_local("matrices"))
    for(k=1;k<=cols(vars);k++) {
        name=vars[k]; a=st_matrix("r("+name+")")
        rn=st_matrixrowstripe("r("+name+")"); cn=st_matrixcolstripe("r("+name+")")
        for(j=1;j<=rows(a);j++) {
            real scalar c
            for(c=1;c<=cols(a);c++) {
                rl=rows(rn)>=j ? rn[j,2] : ""
                cl=rows(cn)>=c ? cn[c,2] : ""
                fput(f,spec+","+name+","+strofreal(j)+","+strofreal(c)+","+strofreal(a[j,c],"%24.17g")+","+rl+","+cl)
            }
        }
    }
    vars=tokens(st_local("scalars"))
    for(k=1;k<=cols(vars);k++) {
        name=vars[k]
        fput(f,spec+","+name+",1,1,"+strofreal(st_numscalar("r("+name+")"),"%24.17g")+",,")
    }
    fclose(f)
}
end
file open f using output/stata_mixed.csv, write replace
file write f "specification,statistic,row,column,value,row_label,column_label" _n
file close f
file open f using output/stata_mixed_keys.csv, write replace
file write f "specification,dimension,index,identifier1,identifier2" _n
file close f

foreach spec in two_full two_strings two_fe two_absorb two_identical two_nested two_singular_nuisance two_negative {
    local dataset synthetic
    if "`spec'"=="two_negative" local dataset negative
    import delimited validation/fixtures/`dataset'.csv, clear asdouble
    local controls "z w"
    local clusters "g h"
    local options ""
    if "`spec'"=="two_strings" {
        local clusters "gs hs"
        keep if keep==1
    }
    if "`spec'"=="two_fe" local options "fevar(g h)"
    if "`spec'"=="two_absorb" local options "absorb(a)"
    if "`spec'"=="two_identical" local clusters "g g"
    if "`spec'"=="two_nested" local clusters "region g"
    if "`spec'"=="two_singular_nuisance" local controls "z zdup"
    if "`spec'"=="two_negative" local controls "z"
    foreach v in y x `controls' `clusters' {
        quietly drop if missing(`v')
    }
    quietly summclust y x `controls', cluster(`clusters') `options' rho(.5) nograph rankcheck
    matrix expected=r(cv3_g)+r(cv3_h)-r(cv1_intersection)
    matrix oldb=r(b)
    matrix oldcv1=r(cv1)
    matrix oldg=r(cv3_g)
    matrix oldh=r(cv3_h)
    matrix oldi=r(cv1_intersection)
    local icount=r(I)
    assert r(intersection_delete_one_solves)==`icount'
    summclust y x `controls', cluster(`clusters') `options' rho(.5) nograph rankcheck table regtable mixed
    assert r(intersection_delete_one_solves)==0
    assert "`r(selected_estimator)'"=="CV31"
    assert mreldif(r(cv31_conventional),expected)<1e-13
    assert mreldif(r(b),oldb)<1e-13
    assert mreldif(r(cv1),oldcv1)<1e-13
    assert mreldif(r(cv3_g),oldg)<1e-13
    assert mreldif(r(cv3_h),oldh)<1e-13
    assert mreldif(r(cv1_intersection),oldi)<1e-13
    mata: assert(all(st_matrix("r(cv3)"):==.))
    mata: assert(all(st_matrix("r(cv3_intersection)"):==.))
    mata: assert(all(st_matrix("r(betanog_intersection)"):==.))
    assert r(df_r)==min(r(G),r(H))-1
    assert reldif(r(cv31_max_var),max(el(expected,1,1),el(oldg,1,1),el(oldh,1,1)))<1e-13
    assert reldif(r(cv31_max_p),2*ttail(r(df_r),abs(r(cv31_max_t))))<1e-13
    assert reldif(r(cv31_max_lci),el(oldb,1,1)-invttail(r(df_r),.025)*r(cv31_max_se))<1e-13
    if "`spec'"=="two_negative" {
        assert r(cv31_conventional_var)<0
        assert missing(r(cv31_conventional_se))
        assert r(cv31_max_se)>0 & r(cv31_max_se)<.
    }
    local matrices : r(matrices)
    local scalars : r(scalars)
    mata: mx_write("`spec'")
    local cg : word 1 of `clusters'
    local ch : word 2 of `clusters'
    foreach dim in g h intersection {
        preserve
        local keyvars "`c`dim''"
        if "`dim'"=="intersection" local keyvars "`cg' `ch'"
        local keyvars : list uniq keyvars
        keep `keyvars'
        duplicates drop
        sort `keyvars'
        file open f using output/stata_mixed_keys.csv, write append
        forvalues j=1/`=_N' {
            local id1=""
            local id2=""
            local n=0
            foreach v of local keyvars {
                local ++n
                local id`n'=`v'[`j']
            }
            if "`dim'"=="intersection" & "`cg'"=="`ch'" local id2="`id1'"
            file write f "`spec',`dim',`j',`id1',`id2'" _n
        }
        file close f
        restore
    }
}

import delimited validation/fixtures/synthetic.csv, clear asdouble
capture summclust y x z, cluster(g) mixed nograph
assert _rc==198
capture summclust y x z, cluster(g h) mixed jackknife nograph
assert _rc==198
quietly summclust y x z w, cluster(g h) mixed nograph
matrix swapref=r(cv31)
matrix scoretarget=r(cv1_intersection)
mata:
    xx=(st_data(.,"x"),J(st_nobs(),1,1),st_data(.,("z","w")))
    yy=st_data(.,"y")
    ai=invsym(quadcross(xx,xx))
    uu=yy-xx*(ai*quadcross(xx,yy))
    codes=1000*st_data(.,"g")+st_data(.,"h")
    oo=order(codes,1)
    inf=panelsetup(codes[oo],1)
    ss=panelsum((xx:*uu)[oo,.],inf)
    ii=rows(inf); nn=rows(xx); kk=cols(xx)
    manual=(ii/(ii-1))*((nn-1)/(nn-kk))*(ai*quadcross(ss,ss)*ai)
    assert(max(abs(manual-st_matrix("scoretarget")))<1e-12)
end
assert r(intersection_delete_one_solves)==0
mata: assert(all(st_matrix("r(rankfail_intersection)"):==.))
quietly summclust y x z w, cluster(h g) mixed nograph
assert mreldif(r(cv31),swapref)<1e-13
gsort -g -h
quietly summclust y x z w, cluster(g h) mixed nograph
assert mreldif(r(cv31),swapref)<1e-11
quietly summclust y z x w, cluster(g h) nograph
matrix otherfocal=r(cv3_g)+r(cv3_h)-r(cv1_intersection)
quietly summclust y z x w, cluster(g h) mixed nograph
assert mreldif(r(cv31),otherfocal)<1e-13
quietly summclust y x z w, cluster(g g) mixed nograph
matrix identity=2*r(cv3_g)-r(cv1_g)
assert mreldif(r(cv31),identity)<1e-13
display "STATA MIXED CERTIFICATION PASSED"
log close
exit, clear
