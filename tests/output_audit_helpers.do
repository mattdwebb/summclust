/* Independent references and complete r() snapshots. No package mathematics
   are reused: regress/predict supply the reference fits and hat diagonals. */
capture program drop audit_snapshot
program define audit_snapshot
    args destination
    local scalars : r(scalars)
    local macros : r(macros)
    local matrices : r(matrices)
    tempname handle value
    file open `handle' using "`destination'", write text replace
    foreach name of local scalars {
        file write `handle' "scalar `name' " %24.17e (r(`name')) _n
    }
    foreach name of local macros {
        file write `handle' "macro `name' `r(`name')'" _n
    }
    foreach name of local matrices {
        matrix `value'=r(`name')
        local cols : colfullnames `value'
        local rows : rowfullnames `value'
        file write `handle' "matrix `name' `rows' : `cols'" _n
        forvalues i=1/`=rowsof(`value')' {
            forvalues j=1/`=colsof(`value')' {
                file write `handle' "value `name' `i' `j' " %24.17e (el(`value',`i',`j')) _n
            }
        }
    }
    file close `handle'
end

capture program drop audit_model
program define audit_model
    syntax varlist(min=3), CLuster(varlist min=2 max=2) ID(string) OUT(string) [FEvar(varlist)]
    gettoken y rhs : varlist
    gettoken x controls : rhs
    gettoken g h : cluster
    local factors ""
    foreach f of local fevar {
        local factors "`factors' i.`f'"
    }
    marksample used
    markout `used' `cluster' `fevar', strok
    quietly summclust `varlist', cluster(`cluster') fevar(`fevar') nograph
    audit_snapshot "`out'/`id'_returns.txt"
    tempname b ds df max1 max3 L1 U1 L3 U3
    matrix `b'=r(b)
    matrix `ds'=r(dimsummary)
    scalar `df'=r(df_r)
    scalar `max1'=r(cv1_max_se)
    scalar `max3'=r(cv3_max_se)
    scalar `L1'=real("`r(cv1lci)'")
    scalar `U1'=real("`r(cv1uci)'")
    scalar `L3'=real("`r(cv3lci)'")
    scalar `U3'=real("`r(cv3uci)'")
    assert reldif(real("`r(cv1t)'"),el(`b',1,1)/`max1')<1e-13
    assert reldif(real("`r(cv3t)'"),el(`b',1,1)/`max3')<1e-13
    assert reldif(real("`r(cv1p)'"),2*ttail(`df',abs(el(`b',1,1)/`max1')))<1e-13
    assert reldif(real("`r(cv3p)'"),2*ttail(`df',abs(el(`b',1,1)/`max3')))<1e-13
    foreach t in g h intersection {
        tempname pl_`t' lv_`t' b_`t' v1_`t' v3_`t'
        matrix `pl_`t''=r(partlev_`t')
        matrix `lv_`t''=r(lever_`t')
        matrix `b_`t''=r(betanog_`t')
        matrix `v1_`t''=r(cv1_`t')
        matrix `v3_`t''=r(cv3_`t')
    }
    quietly regress `y' `x' `controls' `factors' if `used'
    tempname bf max_pl_error max_l_error max_b_error max_v_error
    scalar `bf'=_b[`x']
    assert abs(`bf'-el(`b',1,1))<1e-8
    tempvar hat resx sq ig ih ii
    quietly predict double `hat' if `used', hat
    quietly regress `x' `controls' `factors' if `used'
    quietly predict double `resx' if `used', residuals
    quietly generate double `sq'=`resx'^2
    quietly summarize `sq', meanonly
    local ss=r(sum)
    quietly egen long `ig'=group(`g') if `used'
    quietly egen long `ih'=group(`h') if `used'
    quietly egen long `ii'=group(`cluster') if `used'
    scalar `max_pl_error'=0
    scalar `max_l_error'=0
    scalar `max_b_error'=0
    scalar `max_v_error'=0
    local index=0
    foreach t in g h intersection {
        local ++index
        local code : word `index' of `ig' `ih' `ii'
        quietly regress `y' `x' `controls' `factors' if `used', vce(cluster `code')
        assert abs(_se[`x']^2-el(`v1_`t'',1,1))<1e-8
        tempvar pl lv first
        quietly egen double `pl'=total(`sq') if `used', by(`code')
        quietly replace `pl'=`pl'/`ss'
        quietly egen double `lv'=total(`hat') if `used', by(`code')
        quietly egen byte `first'=tag(`code') if `used'
        quietly summarize `pl' if `first'==1
        assert abs(r(sd)/abs(r(mean))-el(`ds',`index',3))<1e-9
        preserve
        quietly keep if `first'==1
        sort `code'
        mata: st_numscalar("`max_pl_error'",max((st_numscalar("`max_pl_error'"),max(abs(st_data(.,"`pl'")-st_matrix("`pl_`t''"))))))
        mata: st_numscalar("`max_l_error'",max((st_numscalar("`max_l_error'"),max(abs(st_data(.,"`lv'")-st_matrix("`lv_`t''"))))))
        mata: st_numscalar("__audit_gstar",1/sum(st_data(.,"`pl'"):^2))
        assert abs(__audit_gstar-el(`ds',`index',6))<1e-7
        restore
        quietly summarize `code', meanonly
        local J=r(max)
        tempname v3
        scalar `v3'=0
        forvalues j=1/`J' {
            quietly regress `y' `x' `controls' `factors' if `used' & `code'!=`j'
            scalar `max_b_error'=max(`max_b_error',abs(_b[`x']-el(`b_`t'',`j',1)))
            scalar `v3'=`v3'+(_b[`x']-`bf')^2*(`J'-1)/`J'
        }
        scalar `max_v_error'=max(`max_v_error',abs(`v3'-el(`v3_`t'',1,1)))
    }
    assert `max_pl_error'<1e-9 & `max_l_error'<1e-8
    assert `max_b_error'<1e-7 & `max_v_error'<1e-8
    assert abs(`L1'-(`bf'-invttail(`df',.025)*`max1'))<1e-8
    assert abs(`U1'-(`bf'+invttail(`df',.025)*`max1'))<1e-8
    assert abs(`L3'-(`bf'-invttail(`df',.025)*`max3'))<1e-8
    assert abs(`U3'-(`bf'+invttail(`df',.025)*`max3'))<1e-8
    display "MODEL PASSED `id': beta=" %16.12f `bf' " CV1max=" %16.12f `max1' " CV3max=" %16.12f `max3'
    display "Maximum errors: PL=" %10.3e `max_pl_error' " leverage=" %10.3e `max_l_error' " deletion=" %10.3e `max_b_error' " CV3=" %10.3e `max_v_error'
end
