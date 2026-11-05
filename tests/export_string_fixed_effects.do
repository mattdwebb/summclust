version 13
clear all
set more off
capture log close _all
capture mkdir output
log using output/stata_string_fe_exports.log, text replace
run Stata/summclust.ado
mata:
void sfe_write(string scalar spec) {
    external real matrix beta
    real scalar f,j,k,c
    real matrix a
    string rowvector vars
    string scalar name
    string matrix rn,cn
    f=fopen("output/stata_string_fe.csv","a")
    vars=tokens(st_local("matrices"))
    for(k=1;k<=cols(vars);k++) {
        name=vars[k]; a=st_matrix("r("+name+")")
        if (name=="b" & st_local("mode")=="one") a=beta'
        rn=st_matrixrowstripe("r("+name+")"); cn=st_matrixcolstripe("r("+name+")")
        if (name=="b" & st_local("mode")=="one") {
            rn=J(1,2,""); cn=st_matrixcolstripe("r(cv1)")
        }
        for(j=1;j<=rows(a);j++) for(c=1;c<=cols(a);c++) {
            fput(f,spec+","+name+","+strofreal(j)+","+strofreal(c)+","+strofreal(a[j,c],"%24.17g")+","+rn[j,2]+","+cn[c,2])
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
file open f using output/stata_string_fe.csv, write replace
file write f "specification,statistic,row,column,value,row_label,column_label" _n
file close f
file open f using output/stata_string_fe_keys.csv, write replace
file write f "specification,dimension,index,identifier1,identifier2" _n
file close f
foreach mode in one two mixed {
    foreach design in fe absorb both {
        local spec `mode'_`design'
        import delimited validation/fixtures/synthetic.csv, clear asdouble
        generate str30 fstr="Effect / "+string(fe,"%02.0f")
        generate str30 hstr="Period "+string(h,"%02.0f")
        generate str30 astr="Cell: "+string(a,"%04.0f")
        replace fstr="" if mod(row,19)==0
        replace hstr="" if mod(row,23)==0
        replace astr="" if mod(row,29)==0
        replace fstr="Excluded level" if keep==0
        local fes ""
        local ab ""
        if "`design'"=="fe" local fes "fstr hstr"
        if "`design'"=="both" local fes "fstr"
        if "`design'"!="fe" local ab "astr"
        local clusters "gs hs"
        local options ""
        if "`mode'"=="one" {
            local clusters "gs"
            local options "jackknife"
        }
        if "`mode'"=="mixed" local options "mixed"
        summclust y x z w, cluster(`clusters') fevar(`fes') absorb(`ab') sample(keep==1) rho(.5) rankcheck table nograph `options'
        local matrices : r(matrices)
        local scalars : r(scalars)
        if "`mode'"=="one" local matrices "`matrices' b"
        mata: sfe_write("`spec'")
        keep if keep==1
        foreach variable in y x z w `clusters' `fes' `ab' {
            drop if missing(`variable')
        }
        local dimensions "g"
        if "`mode'"!="one" local dimensions "g h i"
        foreach dim of local dimensions {
            preserve
            local keyvars "gs"
            if "`dim'"=="h" local keyvars "hs"
            if "`dim'"=="i" local keyvars "gs hs"
            keep `keyvars'
            duplicates drop
            sort `keyvars'
            file open f using output/stata_string_fe_keys.csv, write append
            forvalues j=1/`=_N' {
                local first : word 1 of `keyvars'
                local id1=`first'[`j']
                local id2=""
                if "`dim'"=="i" local id2=hs[`j']
                file write f "`spec',`dim',`j',`id1',`id2'" _n
            }
            file close f
            restore
        }
    }
}
display "STRING FIXED EFFECT EXPORTS PASSED"
log close
exit, clear
