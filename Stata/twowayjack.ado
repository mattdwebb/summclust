*! version 2.0.0 23aug2026
version 13

capture program drop twowayjack
program define twowayjack, rclass
    version 13
    syntax varlist(min=3) [if] [in], CLuster(varlist min=2 max=2) ///
        [FEVAR(varlist) ABSORB(varname) GSTAR Rho(real -1) SAMple(string) ///
         TABle ADDmeans JACKknife REGtable NOGraph RANKcheck]

    marksample touse, strok
    tempvar finalsample
    quietly generate byte `finalsample' = `touse'
    if `"`sample'"' != "" quietly replace `finalsample' = 0 if !(`sample')

    local opts "cluster(`cluster') sample(`finalsample'==1)"
    if `"`fevar'"' != "" local opts `"`opts' fevar(`fevar')"'
    if `"`absorb'"' != "" local opts `"`opts' absorb(`absorb')"'
    if `"`gstar'"' != "" local opts `"`opts' gstar"'
    if `rho' != -1 local opts `"`opts' rho(`rho')"'
    if `"`table'"' != "" local opts `"`opts' table"'
    if `"`addmeans'"' != "" local opts `"`opts' addmeans"'
    if `"`jackknife'"' != "" local opts `"`opts' jackknife"'
    if `"`regtable'"' != "" local opts `"`opts' regtable"'
    if `"`nograph'"' != "" local opts `"`opts' nograph"'
    if `"`rankcheck'"' != "" local opts `"`opts' rankcheck"'

    summclust `varlist', `opts'
    return add
end
