version 13
clear all
set more off
capture log close _all
log using verify_net_install_twoway_results.log, replace text
capture mkdir __net_install_tw
local source "`c(pwd)'"
local oldplus "`c(sysdir_plus)'"
sysdir set PLUS "`source'/__net_install_tw"
net install summclust, from("`source'") replace
/* Change directory so the ado files in the source do not shadow PLUS. */
cd __net_install_tw
discard
which summclust
which twowayjack
set seed 51
set obs 120
generate int g=ceil(_n/12)
generate int h=mod(_n-1,6)+1
generate double x=rnormal()
generate double z=rnormal()
generate double y=.5*x-.2*z+rnormal()
quietly summclust y x z, cluster(g h) nograph
assert r(G)==10 & r(H)==6
quietly twowayjack y x z if mod(_n,2), cluster(g h) nograph
assert r(G)==10 & r(H)==3
display as result _newline "LOCAL NET INSTALL TWO-WAY VERIFICATION PASSED"
sysdir set PLUS "`oldplus'"
cd "`source'"
log close
exit, clear
