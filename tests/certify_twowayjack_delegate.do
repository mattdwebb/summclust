version 13
clear all
set more off
capture log close _all
log using certify_twowayjack_delegate_results.log, replace text
run summclust.ado
run twowayjack.ado
set seed 5102
set obs 360
generate int g=ceil(_n/24)
generate int h=mod(_n-1,12)+1
generate double x=rnormal()
generate double z=rnormal()
generate double y=.6*x-.2*z+rnormal()
generate byte keep=mod(_n,4)!=0
quietly summclust y x z, cluster(g h) sample(keep==1) nograph rankcheck
tempname V1 V3 D
matrix `V1'=r(cv1)
matrix `V3'=r(cv3)
matrix `D'=r(dimsummary)
local se1=real("`r(cv1se)'")
local se3=real("`r(cv3se)'")
local GG=r(G)
local HH=r(H)
quietly twowayjack y x z if keep, cluster(g h) nograph rankcheck
mata: st_numscalar("__dv1",max(abs(st_matrix("r(cv1)"):-st_matrix("`V1'")))); st_numscalar("__dv3",max(abs(st_matrix("r(cv3)"):-st_matrix("`V3'")))); st_numscalar("__dd",max(abs(st_matrix("r(dimsummary)"):-st_matrix("`D'"))))
assert scalar(__dv1)<1e-10 & scalar(__dv3)<1e-10 & scalar(__dd)<1e-10
assert abs(real("`r(cv1se)'")-`se1')<1e-10
assert abs(real("`r(cv3se)'")-`se3')<1e-10
assert r(G)==`GG' & r(H)==`HH'
display as result _newline "THIN TWOWAYJACK DELEGATE CERTIFICATION PASSED"
log close
exit, clear
