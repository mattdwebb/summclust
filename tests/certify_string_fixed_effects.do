version 13
clear all
set more off
capture log close _all
log using certify_string_fixed_effects_results.log, text replace
run summclust.ado

/* Snapshot every public numeric result. Coefficient names are checked
   separately because the numeric and string input columns have different names. */
capture program drop sf_save
program define sf_save
    global SF_matrices : r(matrices)
    global SF_scalars : r(scalars)
    global SF_macros : r(macros)
    local j=0
    foreach name of global SF_matrices {
        local ++j
        matrix SF_m`j'=r(`name')
    }
    local j=0
    foreach name of global SF_scalars {
        local ++j
        scalar SF_s`j'=r(`name')
    }
    local j=0
    foreach name of global SF_macros {
        local ++j
        global SF_t`j' `"`r(`name')'"'
    }
end

capture program drop sf_compare
program define sf_compare
    local matrices : r(matrices)
    local scalars : r(scalars)
    local macros : r(macros)
    assert "`matrices'"=="$SF_matrices"
    assert "`scalars'"=="$SF_scalars"
    assert "`macros'"=="$SF_macros"
    local j=0
    foreach name of local matrices {
        local ++j
        assert mreldif(r(`name'),SF_m`j')<1e-11
    }
    local j=0
    foreach name of local scalars {
        local ++j
        assert r(`name')==SF_s`j' | reldif(r(`name'),SF_s`j')<1e-11
    }
    local j=0
    foreach name of local macros {
        local ++j
        assert `"`r(`name')'"'==`"${SF_t`j'}"'
    }
end

set seed 41126
set obs 384
generate long rowid=_n
generate int g=ceil(_n/32)
generate int h=mod(floor((_n-1)/4),8)+1
generate int f=mod(_n-1,4)+1
generate int f2=mod(floor((_n-1)/2),3)+1
egen int a=group(g h)
generate str8 gs="c"+string(g,"%02.0f")
generate str24 fs=cond(f==1,".",cond(f==2,"Inf",cond(f==3,"nan","quote's [space]")))
generate str8 f2s=cond(f2==1,"a b",cond(f2==2,"a-b","a/b"))
generate str20 astr="group "+string(a,"%03.0f")+" !"
generate double x=rnormal()
generate double z=rnormal()
generate double w=rnormal()
generate double y=.7*x-.2*z+.1*w+f/5+rnormal()
generate byte keep=mod(rowid,11)!=0
replace fs="00 excluded" if !keep
replace f=99 if !keep
replace f2s="00 excluded" if !keep
replace f2=99 if !keep
replace astr="00 excluded" if !keep
replace a=999 if !keep
replace fs="" in 17
replace f=. in 17
replace f2s="" in 19
replace f2=. in 19
replace astr="" in 23
replace a=. in 23
replace z=. in 29
replace y=. in 31
clonevar long_fixed_effect_identifier_aaa=fs
clonevar long_fixed_effect_identifier_aab=f2s
generate double shuffle=runiform()
sort shuffle
generate long incoming=_n
quietly datasignature
local data_signature "`r(datasignature)'"

foreach spec in feone femultiple femixedtypes absorb both nested long twoway twomixed twoabsorb twoboth twolong {
    local cluster gs
    local numeric "fevar(f)"
    local strings "fevar(fs)"
    local extra "jackknife regtable"
    if "`spec'"=="femultiple" {
        local numeric "fevar(f f2)"
        local strings "fevar(fs f2s)"
    }
    if "`spec'"=="femixedtypes" {
        local numeric "fevar(f f2)"
        local strings "fevar(fs f2)"
    }
    if inlist("`spec'","absorb","twoabsorb") {
        local numeric "absorb(a)"
        local strings "absorb(astr)"
    }
    if inlist("`spec'","both","twoboth") {
        local numeric "fevar(f f2) absorb(a)"
        local strings "fevar(fs f2s) absorb(astr)"
    }
    if "`spec'"=="nested" {
        local numeric "fevar(g)"
        local strings "fevar(gs)"
    }
    if inlist("`spec'","long","twolong") {
        local numeric "fevar(f f2)"
        local strings "fevar(long_fixed_effect_identifier_aaa long_fixed_effect_identifier_aab)"
    }
    if substr("`spec'",1,3)=="two" {
        local cluster "gs h"
        local extra "regtable"
    }
    if "`spec'"=="twomixed" local extra "mixed regtable"
    set sortseed 41126
    quietly summclust y x z w, cluster(`cluster') `numeric' `extra' sample(keep) rho(.4) table rankcheck addmeans nograph
    sf_save
    if substr("`spec'",1,3)!="two" matrix SF_reg=regresstab
    assert incoming==_n
    set sortseed 41126
    quietly summclust y x z w, cluster(`cluster') `strings' `extra' sample(keep) rho(.4) table rankcheck addmeans nograph
    sf_compare
    if substr("`spec'",1,3)!="two" assert mreldif(regresstab,SF_reg)<1e-11
    if "`spec'"=="feone" {
        matrix actual=r(cv3)
        local names : colnames actual
        assert "`names'"=="x z w fs_1 fs_2 fs_3 fs_4"
    }
    if inlist("`spec'","long","twolong") {
        matrix actual=r(cv1)
        local names : colnames actual
        local expected "x z w long_fixed_effect_identifie_f1_1 long_fixed_effect_identifie_f1_2 long_fixed_effect_identifie_f1_3 long_fixed_effect_identifie_f1_4 long_fixed_effect_identifie_f2_1 long_fixed_effect_identifie_f2_2"
        assert "`names'"=="`expected'"
        local unique : list uniq names
        assert "`names'"=="`unique'"
    }
    assert incoming==_n
    quietly datasignature
    assert "`r(datasignature)'"=="`data_signature'"
    display "PASS numeric/string equivalence and unchanged input: `spec'"
}

/* String absorption must retain the same nesting restrictions and restore order. */
capture quietly summclust y x z w, cluster(gs) absorb(fs) sample(keep) nograph
assert _rc==198
assert incoming==_n
capture quietly summclust y x z w, cluster(gs h) absorb(gs) sample(keep) nograph
assert _rc==198
assert incoming==_n

/* Mapping uses the complete-case sample, not extra categories outside it. */
quietly summclust y x z w, cluster(gs) fevar(fs) sample(keep) nograph
matrix sample_cv3=r(cv3)
matrix sample_pl=r(partlev)
preserve
quietly keep if keep & !missing(y,x,z,w,gs,fs)
quietly summclust y x z w, cluster(gs) fevar(fs) nograph
assert mreldif(r(cv3),sample_cv3)<1e-11
assert mreldif(r(partlev),sample_pl)<1e-11
restore

/* Exercise the default graph path with string fixed effects. */
quietly summclust y x z w, cluster(gs) fevar(fs) sample(keep)
graph drop _all
assert incoming==_n
quietly datasignature
assert "`r(datasignature)'"=="`data_signature'"

display "SUMMCLUST STRING FIXED EFFECTS CERTIFICATION PASSED"
log close
exit, clear
