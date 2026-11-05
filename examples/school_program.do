* Run from the top-level summclust folder after downloading the repository.
version 13
clear all
set more off
set linesize 100
set sortseed 11052026
capture mkdir output
log using output/school_program_stata.log, text replace
adopath ++ "Stata"
import delimited using "examples/school_program.csv", clear

* The data are fictional: student observations in 24 schools over five years.
* The coefficient of interest is on program, the second variable below.
summclust test_score program baseline_score, cluster(school) fevar(year) table
graph export "output/school_program.png", width(1600) replace

display "SCHOOL PROGRAM EXAMPLE COMPLETED"
log close
