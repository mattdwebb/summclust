*! version 5.1.1 23sep2026
/*------------------------------------*/
/*summclust */
/*written by Matt Webb */
/*version 5.1.1 2026-09-23 */
/*------------------------------------*/
version 13

/* Optimized implementation; validation scope is documented in docs/AUDIT_2026-09-23.md. */

cap program drop summclust
program define summclust, rclass

	local __summclust_cmdline `"`0'"'
	syntax varlist(min=3), cluster(varlist min=1 max=2) [fevar(varlist) absorb(varname) gstar  Rho(real -1) SAMple(string)  TABle ADDmeans JACKknife REGtable NOGraph RANKcheck]

	local __summclust_ncluster : word count `cluster'
	if `__summclust_ncluster' == 2 {
		_summclust_twoway `__summclust_cmdline'
		return add
		exit
	}

	local numvars = wordcount("`varlist'")
	local yvar = word("`varlist'",1)

	local CTEMP ""

	local start = 1 + `numvars' + 1

	/*split indepvars into variables of interest and controls*/
    local XTEMP = word("`varlist'",2)

    forvalues i = 3/`numvars' {

      local wordi = word("`varlist'",`i')

      local CTEMP = "`CTEMP' `wordi' "
    }

	local FEVAR "`fevar'"

	local SMPL "`sample'"

	local svars = "`addmeans'"

	/*add if to sample*/
	if "`sample'" != "" {
	    local SMPL = "if `sample'"
	}

	/* Check rho and distinguish an omitted option from every valid value. */
	local rho_specified = (`rho' != -1)
	if `rho_specified' {
		if inrange(`rho',0.00001,1)==0 {
			disp as error "The value for rho, `rho', is invalid. rho must be between 0 and 1."
			exit 198
		}
		local gstarind = 1
	}
	else {
		/* Retain the historical internal default when rho() is omitted. */
		local rho = 0.1234567
		local gstarind = 0
	}

	/*flag to calculate g-star*/
	if "`gstar'" != "" {
		local gstarind = 1
	}

	/*flag for absorb used properly*/
	local nocv3 = 0

	/*flags for singular clusters*/
	local singclust = 0
	local singclustfull = 0
	local colsdrop = 0
	local countfevar = 0

	/*declare tempvars*/
	tempvar ytilde
	tempvar temp_sd

	local fevars " "

	/*xvar labels for regtable*/
	local name_xvar = "`CTEMP'"

	local catvars = " "
	foreach catvar in `FEVAR' {
	    local catvars = " `catvars' " +  "i.`catvar'"
	}


	local absvars = " "
	foreach avar in `absorb' {
	    local absvars = " `absvars' " +  "i.`avar'"
	}

	tempvar temp_sample
	mark `temp_sample'
	markout `temp_sample' `yvar' `XTEMP' `CTEMP' `catvars' `absvars' `cluster',strok

	/*impose the "sample option" restriction*/
	if "`sample'" != "" {

		tempvar temp_samp
		qui gen `temp_samp' = 0
		qui replace `temp_samp' = 1 `SMPL'
		qui replace `temp_sample' = `temp_sample'*`temp_samp'
	}

	quietly count if `temp_sample'==1
	if r(N)==0 error 2000

	/*for fevar*/
	if `"`fevar'"' != "" {
		/* Build only coefficient labels here. The indicator matrix itself */
		/* is created once, preallocated, in Mata after the sample is sorted. */
		local names_fe = ""
		local q = 0
		foreach fvar in `FEVAR' {
			local q = `q' + 1
			local countfevar = `countfevar' + 1
			qui levelsof `fvar' if `temp_sample'==1, local(fevarlevel)
			local nlevels : word count `fevarlevel'
			local fpos = 0
			foreach flevel in `fevarlevel' {
				local fpos = `fpos' + 1
				if (`q'==1 & `"`absorb'"'=="") | `fpos'<`nlevels' {
					local names_fe "`names_fe' `fvar'_`flevel'"
				}
			}
		}
		local XVAR = "`CTEMP'"
	}

	/* Column order of the Mata design and all full covariance matrices. */
	if `"`fevar'"' == "" & `"`absorb'"' == "" {
		local coefnames "`XTEMP' _cons `CTEMP'"
	}
	else {
		local coefnames "`XTEMP' `CTEMP'"
		if `"`fevar'"' != "" local coefnames "`coefnames' `names_fe'"
	}

	/*sort the data by the clustering variable*/
		/*make a sequential-numeric cluster indicator*/

	tempvar temp_indexa temp_order
	qui gen long `temp_order' = _n
	qui egen `temp_indexa' = group(`cluster') if `temp_sample'==1
	qui summ `temp_indexa'
	local G = r(max)
	if `G' < 2 {
		display as error "summclust requires at least two clusters in the estimation sample."
		exit 498
	}

	/* Keep absorb groups contiguous within clusters. This permits fast, */
	/* in-place demeaning below and fixes ordering when absorb is nested.  */
	if `"`absorb'"' != "" {
		qui sort `temp_indexa' `absorb' `temp_order'
	}
	else {
		qui sort `temp_indexa' `temp_order'
	}
	qui putmata CVAR=`temp_indexa' if `temp_sample'==1, replace

	mata: numobs = rows(CVAR)

	mata: info = panelsetup(CVAR,1)
	mata: ng = info[.,2] - info[.,1] :+ 1

	mata: xtilde= J(numobs,0,.)
	local xtilde = " "

	local XVAR = "`CTEMP'"

	/*dump everything into mata here, then demean by absorb*/

	qui gen double `ytilde' = `yvar' if `temp_sample'==1

	foreach xvar in `XTEMP' {

		qui putmata `xvar' = `xvar' if `temp_sample'==1, replace

		mata: xtilde = xtilde, `xvar'

		local xtilde = "`xtilde'" + " ``varname''"

	}


	if `"`absorb'"' != "" {

		/* Verify that cluster is constant within absorb. Use the numeric */
		/* cluster index so string cluster variables are also supported. */
		cap drop `temp_sd'
		qui bysort `absorb': egen `temp_sd' = sd(`temp_indexa') if `temp_sample'==1
		qui replace `temp_sd' = 0 if `temp_sample'==1 & missing(`temp_sd')
		qui summarize `temp_sd' if `temp_sample'==1, meanonly
		if r(max) > 1e-10 {
			disp as error "Cluster variable not constant within absorb variable.  Use fevar instead."
			exit 198
		}

		/* Restore the exact row order used to form CVAR and xtilde. */
		qui sort `temp_indexa' `absorb' `temp_order'
		qui putmata absorb = `absorb' if `temp_sample'==1, replace
	}

	/*maybe make this for when there is neither absorb and fevar*/
	if `"`absorb'"' == "" & `"`fevar'"' == "" {
		/*add constant*/
		tempvar tempcons
		qui gen `tempcons'  = 1
		local xtilde = "`xtilde'" + " `tempcons' "

		mata: ones = J(rows(xtilde),1,1)
		mata: xtilde = xtilde, ones
	}

	/*put everything in mata*/
		if `"`absorb'"' != "" {
			qui sort `temp_indexa' `absorb' `temp_order'
		}
		else {
			qui sort `temp_indexa' `temp_order'
		}

		qui putmata Y = `ytilde' if `temp_sample'==1, replace

		qui putmata XOTHER = (`CTEMP') if `temp_sample'==1, replace

		if `"`fevar'"' != "" {
			qui putmata FECODES = (`FEVAR') if `temp_sample'==1, replace
			local fe_has_absorb = (`"`absorb'"' != "")
			mata: FEDUM = _summclust_dummies(FECODES,`fe_has_absorb')
			mata: X = xtilde, XOTHER, FEDUM
		}
		else {
			mata: X = xtilde, XOTHER
		}

		mata: k = cols(X)
		mata: st_numscalar("k", k)
		mata: st_local("k", strofreal(k))
		mata: numobs = rows(X)
		mata: st_numscalar("numobs", numobs)
		mata: st_local("numobs", strofreal(numobs))

		mata: sk = rank(X)
		mata: st_numscalar("sk", sk)
		mata: st_local("sk", strofreal(sk))

		local coeffnum = `sk'

		mata: iota =J(numobs,1,1)

	/*for absorb option*/
	if `"`absorb'"' != "" {

		/*need to partition everything by absorb variable*/
		mata:  abinfo = panelsetup(absorb,1)

		mata : abng = abinfo[.,2] - abinfo[.,1] :+ 1

		mata: numab = rows(abinfo)

		mata: st_numscalar("H", numab)
		mata: st_local("H", strofreal(numab))

		local coeffnum = `sk' + `H'

		/* Demean in place. The original repeatedly appended matrices, */
		/* which copies the accumulated data at every absorb group.     */
		mata: Y = _summclust_within(Y,abinfo)
		mata: X = _summclust_within(X,abinfo)

		mata: st_local("sk", strofreal(rank(X)))
		local coeffnum = `sk' + `H'
	} /*end of  absorb if*/
	if `numobs' <= `coeffnum' {
		quietly sort `temp_order'
		display as error "Residual degrees of freedom are not positive after fixed effects."
		exit 2001
	}


		/*--------------------------------------------------------------*/
		/* Fast core: one Mata loop over clusters, no per-cluster       */
		/* Stata-to-Mata crossings and no dynamically growing matrices. */
		/*--------------------------------------------------------------*/
		local gfact = (`G'*(`numobs'-1))/((`G'-1)*(`numobs'-`coeffnum'))
		local dojack = ("`jackknife'" != "")
		local dogstar = (`gstarind' > 0)
		local dorankcheck = ("`rankcheck'" != "")

		mata: _summclust_core(X,Y,info,`gfact',`dojack',`dogstar',`dorankcheck',`sk')

		/* Build the problem-cluster list without G separate Mata calls. */
		local probclust = ""
		local probclustfull = ""
		tempname __sc_bad
		mata: st_matrix("`__sc_bad'",badcluster')
		forvalues g=1/`G' {
			if el(`__sc_bad',1,`g') != 0 {
				local probclustfull "`probclustfull' `g'"
			}
		}

		local singfrac = (`singclustfull'/`G')*100

		/* CV1 inference */
		local cv1t = `BETATEMP'/`setemp'
		local cv1p = 2*min(ttail(`G'-1,`cv1t'),1-ttail(`G'-1,`cv1t'))
		local cv1lci = `BETATEMP'-invttail(`G'-1,0.025)*`setemp'
		local cv1uci = `BETATEMP'+invttail(`G'-1,0.025)*`setemp'

		/* CV3 inference */
		local cv3t = `BETATEMP'/`cv3se'
		local cv3p = 2*min(ttail(`G'-1,`cv3t'),1-ttail(`G'-1,`cv3t'))
		local cv3lci = `BETATEMP'-invttail(`G'-1,0.025)*`cv3se'
		local cv3uci = `BETATEMP'+invttail(`G'-1,0.025)*`cv3se'

		if "`jackknife'" != "" {
			local cv3jt = `BETATEMP'/`cv3jse'
			local cv3jp = 2*min(ttail(`G'-1,`cv3jt'),1-ttail(`G'-1,`cv3jt'))
			local cv3jlci = `BETATEMP'-invttail(`G'-1,0.025)*`cv3jse'
			local cv3juci = `BETATEMP'+invttail(`G'-1,0.025)*`cv3jse'
		}

		if `singclustfull' >= 1 & `gdrop' > 0 {
			local cv3tdrop = `BETATEMP'/`cv3sedrop'
			local cv3pdrop = 2*min(ttail(`gdrop'-1,`cv3tdrop'),1-ttail(`gdrop'-1,`cv3tdrop'))
			local cv3lcidrop = `BETATEMP'-invttail(`gdrop'-1,0.025)*`cv3sedrop'
			local cv3ucidrop = `BETATEMP'+invttail(`gdrop'-1,0.025)*`cv3sedrop'

			if "`jackknife'" != "" {
				local cv3jtdrop = `BETATEMP'/`cv3jsedrop'
				local cv3jpdrop = 2*min(ttail(`gdrop'-1,`cv3jtdrop'),1-ttail(`gdrop'-1,`cv3jtdrop'))
				local cv3jlcidrop = `BETATEMP'-invttail(`gdrop'-1,0.025)*`cv3jsedrop'
				local cv3jucidrop = `BETATEMP'+invttail(`gdrop'-1,0.025)*`cv3jsedrop'
			}
		}
		else if `singclustfull' >= 1 {
			local cv3tdrop = .
			local cv3pdrop = .
			local cv3lcidrop = .
			local cv3ucidrop = .
			local cv3jtdrop = .
			local cv3jpdrop = .
			local cv3jlcidrop = .
			local cv3jucidrop = .
		}

		local SUMVAR "ng Leverage pLeverage betag"

		/*matrices to store output*/
			mata: clustsum= J(7,4,.)
			mata: bonus = J(6,4,.)

		/* Dropped-subsample matrices were formed in the fast Mata core. */
		if `singclustfull' >= 1 & `colsdrop' > 0 {
			local SUMVAR "ng Leverage pLeverage betag betadropg"
			mata: clustsum = J(7,5,.)
			mata: bonus = J(6,5,.)
		}


		mata: betag = betags[1,.]'

		local s = 0

		foreach svar in `SUMVAR' {

			  local s = `s' + 1

			  local tempsvar = "temp_`svar'"

			  tempvar `tempsvar'

			  qui getmata ``tempsvar'' = `svar', force

			  qui summ ``tempsvar'' in 1/`G', det

			  /*min*/
				local min = r(min)
				mata: clustsum[1,`s'] = `min'

			  /*q1*/
				local q1 = r(p25)
				mata: clustsum[2,`s'] = `q1'

			  /*median*/
				local median = r(p50)
				mata: clustsum[3,`s'] = `median'

			  /*mean*/
				local mean = r(mean)
				mata: clustsum[4,`s'] = `mean'

				local graph_`s' = `mean'

			   /*q3*/
				local q1 = r(p75)
				mata: clustsum[5,`s'] = `q1'

			  /*max*/
				local max = r(max)
				mata: clustsum[6,`s'] = `max'

			  /*coeff of variation*/

				  mata: meandiff = `svar' :- `mean'
				  mata: meandiff2 = meandiff:*meandiff
				  mata: meandiffsum = colsum(meandiff2)
				  mata: denom= 1/((`G'-1)*(`mean'^2))
				  mata: scalvar = denom*meandiffsum
				  mata: scalvar = sqrt(abs(scalvar))

				  mata: clustsum[7,`s'] = scalvar

			  /*svar options*/
			  if "`svars'"!="" {

				/*harmonic mean*/
					mata: sinv = 1:/`svar'
					mata: meaninv = mean(sinv)
					mata: harmonic = 1/meaninv

					mata: bonus[1,`s'] = harmonic

					mata: harmrat = harmonic/`mean'

					mata: bonus[2,`s'] = harmrat

				/*geometric mean*/

					mata: logvar = log(`svar')
					mata: meanlog = mean(logvar)

					mata: geo = exp(meanlog)

					mata: bonus[3,`s'] = geo

					mata: georat = geo/`mean'

					mata: bonus[4,`s'] = georat

				/*quadratic mean*/

					mata: squares = `svar':*`svar'
					mata: sqmean = mean(squares)
					mata: qdmean = sqrt(sqmean)

					mata: bonus[5,`s'] = qdmean

					mata: quadrat = qdmean/`mean'

					mata: bonus[6,`s'] = quadrat

		  } /*end of svars if*/

		 mata: st_matrix("`svar'", `svar')

		} /*end of sum var loop*/

		/*replace harmonic and geometric  values for betanog*/
			mata: bonus[1,4] = .
			mata: bonus[2,4] = .
			mata: bonus[3,4] = .
			mata: bonus[4,4] = .

			mata: st_matrix("clustsum", clustsum)
			matrix rownames clustsum = min q1 median mean q3 max coefvar

			if `singclustfull' >=1  & `colsdrop' >= 1 {
				matrix colnames clustsum = Ng Leverage "Partial L." "all beta no g" "kept beta no g"
					mata: bonus[1,5] = .
					mata: bonus[2,5] = .
					mata: bonus[3,5] = .
					mata: bonus[4,5] = .
			}
			else {
				matrix colnames clustsum = Ng Leverage "Partial L." "beta no g"
			}
			/*rowlabels for cluster by cluster matrix*/
				qui levelsof `cluster' if `temp_sample'==1, local(CLUSTERNAMES)

				mata: cnames = J(`G',1,".")
				local f = 0
				foreach cnames in `CLUSTERNAMES'{
				    local f = `f' + 1
					mata: cnames[`f',1] = "`cnames'"
					local clustlabel`f' = "`cnames'"
				}

		/****************************/
		/*gstar option*/
		/***************************/
			/*check if fixed effects invariant within clusters*/
			if `gstarind' == 1 {

				if `"`fevar'"' != "" {

					foreach fvar in `FEVAR' {

						if "`fvar'" == "`cluster'" {
							local gstarind = 2
						}
						else {

							cap drop  `temp_sd'
							qui bysort `fvar': egen `temp_sd' = sd(`temp_indexa') if `temp_sample'==1
							qui replace `temp_sd' = 0 if `temp_sample'==1 & missing(`temp_sd')
							qui summarize `temp_sd' if `temp_sample'==1, meanonly

							if r(max) <= 1e-10 {
								local gstarind = 2
							}
						}
					} /*end of fvar*/
				} /*end of fevar if*/

				if `"`absorb'"' != "" {

					if "`absorb'" == "`cluster'" {
						local gstarind = 2
					}
					else {

						/*check for invariance of cluster variable*/
					    cap drop  `temp_sd'
						qui bysort `absorb': egen `temp_sd' = sd(`temp_indexa') if `temp_sample'==1
						qui replace `temp_sd' = 0 if `temp_sample'==1 & missing(`temp_sd')
						qui summarize `temp_sd' if `temp_sample'==1, meanonly

						if r(max) <= 1e-10 {
							local gstarind = 2
						}
					}
				} /*end of absorb if*/
			} /*end of if 1*/

			if `gstarind' == 1 {
				mata: _summclust_gstar_all(`rho')
			}
			else if `gstarind' == 2 {
				mata: _summclust_gstar_zero()
			}

		/*output display*/
		disp ""
		disp ""
		disp "SUMMCLUST - MacKinnon, Nielsen, and Webb"
		disp " "
		disp "Cluster summary statistics for `XTEMP' when clustered by `cluster'."
		disp "There are `numobs' observations within `G' `cluster' clusters."

		/*singular clusters warnings*/
		if `singclustfull' >= 1 & `colsdrop' >= 1{
			mata problems = J(0,1,".")
			local clustnames = ""
			foreach problem in `probclustfull' {
				mata: problems = problems \ cnames[`problem',1]
				local clustnames = "`clustnames' " + "`clustlabel`problem'' "
			}

			local singfrac = round(`singfrac',0.01)

			disp " "
			disp as error "***************************************************************************"
			disp "WARNING -  Elements of beta undefined when certain cluster(s) are omitted"
			disp "***************************************************************************"
			disp "Two standard errors are calculated: "
			disp "The first standard error uses a generalized inverse."
			disp "The second standard error drops the singularities."
			disp "***************************************************************************"
			disp "There are `singclustfull' problem clusters, out of `G' clusters."
			disp "The problematic `cluster' cluster(s) are: `clustnames'"

			if `singfrac' < 20 {
				disp "As only" %7.2f `singfrac' " % of subsamples are singular, dropping them is preferred."
			}
			else {
			disp "As" %7.2f `singfrac' "% of subsamples are singular, the generalized inverse is preferred."
			}
			disp " "
			disp " "

		} /*end of singclust if*/

		else if `singclustfull' >= 1 & `colsdrop' == 0{
			mata problems = J(0,1,".")
			local clustnames = ""
			foreach problem in `probclustfull' {
				mata: problems = problems \ cnames[`problem',1]
				local clustnames = "`clustnames' " + "`clustlabel`problem'' "
			}

			local singfrac = round(`singfrac',0.01)

			disp " "
			disp as error "***************************************************************************"
			disp "WARNING -  Elements of beta undefined when certain cluster(s) are omitted"
			disp "***************************************************************************"
			disp "The standard error uses a generalized inverse."
			disp "***************************************************************************"
			disp "There are `singclustfull' problem clusters, out of `G' clusters."
			disp "All clusters are problematic."

			disp " "
			disp " "

		} /*end of singclust else if*/

		if "`jackknife'"!="" {
			matrix define cvstuff = J(3,6,.)
		matrix rownames cvstuff = CV1 CV3 CV3J
			local RSPECCV "&|&&|"
		}
		else {
			matrix define cvstuff = J(2,6,.)
			matrix rownames cvstuff = CV1 CV3
			local RSPECCV "&|&|"
		}

		/*check if absorb used improperly*/
		if `nocv3'==1 {

			matrix cvstuff[1,1] = `BETATEMP'
			matrix cvstuff[1,2] = `setemp'
			matrix cvstuff[1,3] = `cv1t'
			matrix cvstuff[1,4] = `cv1p'
			matrix cvstuff[1,5] = `cv1lci'
			matrix cvstuff[1,6] = `cv1uci'

			matrix cvstuff[2,1] = `BETATEMP'

			disp" "
			disp as error "*********************************************************************"
			disp "WARNING -  `cluster' is not constant within the absorb variable, `absorb' ".
			disp "Beta no g, and CV3 are not computed correctly in this case. "
			disp "Use the fevar option for `absorb' instead, if the above are desired."
			disp "*********************************************************************"
			disp" "

			/*also replace the beta-no-g, partial leverage, and leverage*/
			forvalues i = 1/ 7 {
				matrix clustsum[`i',2] = .
				matrix clustsum[`i',3] = .
				matrix clustsum[`i',4] = .
			}

		} /*end of nocv3 */
		else {

			matrix cvstuff[1,1] = `BETATEMP'
			matrix cvstuff[1,2] = `setemp'
			matrix cvstuff[1,3] = `cv1t'
			matrix cvstuff[1,4] = `cv1p'
			matrix cvstuff[1,5] = `cv1lci'
			matrix cvstuff[1,6] = `cv1uci'

			matrix cvstuff[2,1] = `BETATEMP'
			matrix cvstuff[2,2] = `cv3se'
			matrix cvstuff[2,3] = `cv3t'
			matrix cvstuff[2,4] = `cv3p'
			matrix cvstuff[2,5] = `cv3lci'
			matrix cvstuff[2,6] = `cv3uci'

			if "`jackknife'"!="" {

				matrix cvstuff[3,1] = `BETATEMP'
				matrix cvstuff[3,2] = `cv3jse'
				matrix cvstuff[3,3] = `cv3jt'
				matrix cvstuff[3,4] = `cv3jp'
				matrix cvstuff[3,5] = `cv3jlci'
				matrix cvstuff[3,6] = `cv3juci'

				if `nocv3'== 1 {
					forvalues j= 2 /6 {
						matrix cvstuff[3,`j'] = .
					}
				}
			}
		} /*end else*/

		/*singular clusters*/
		if `singclustfull' >=1 {

			if "`jackknife'" != "" {
				matrix define cvstuffdrop = J(2,6,.)
				matrix rownames cvstuffdrop = CV3 CV3J
				local RSPECCVDROP "&|&|"

				matrix cvstuffdrop[2,1] = `BETATEMP'
				matrix cvstuffdrop[2,2] = `cv3jsedrop'
				matrix cvstuffdrop[2,3] = `cv3jtdrop'
				matrix cvstuffdrop[2,4] = `cv3jpdrop'
				matrix cvstuffdrop[2,5] = `cv3jlcidrop'
				matrix cvstuffdrop[2,6] = `cv3jucidrop'
			}
			else {

				matrix define cvstuffdrop = J(1,6,.)
				matrix rownames cvstuffdrop = CV3
				local RSPECCVDROP "|&|"

			}

			matrix cvstuffdrop[1,1] = `BETATEMP'
			matrix cvstuffdrop[1,2] = `cv3sedrop'
			matrix cvstuffdrop[1,3] = `cv3tdrop'
			matrix cvstuffdrop[1,4] = `cv3pdrop'
			matrix cvstuffdrop[1,5] = `cv3lcidrop'
			matrix cvstuffdrop[1,6] = `cv3ucidrop'

		} /* end of singclustfull if */

		matrix colnames cvstuff = "Coeff" "Sd. Err." "t-stat" "P value" CI-lower CI-upper

		matlist cvstuff, title(Regression Output) rowtitle(s.e.) ///
			cspec(& %-4s w6 | %9.6f w10 & %9.6f & %6.4f w7 & %6.4f w7 & %9.6f w10 & %9.6f w10 &) ///
			rspec("`RSPECCV'")

		if `nocv3'== 1 {
			disp as error "WARNING - CV3 estimates not available as absorb option improperly specified."
		}
		else {

			if `singclustfull' >=1 & `colsdrop' >= 1 {

				matrix colnames cvstuffdrop = "Coeff" "Sd. Err." "t-stat" "P value" CI-lower CI-upper

				matlist cvstuffdrop, title(Regression Output -- Dropping Singular Omit-One-Cluster Subsamples) rowtitle(s.e.) ///
				cspec(& %-4s w6 | %9.6f w10 & %9.6f & %6.4f w7 & %6.4f w7 & %9.6f w10 & %9.6f w10 &) ///
				rspec("`RSPECCVDROP'")
			}
		}

		mata: cvstuff=st_matrix("cvstuff")

		if `singclustfull' >=1  & `colsdrop' >= 1  {
			matlist clustsum , title("Cluster Variability") rowtitle(Statistic) ///
			cspec(& %-10s | %8.2f o2 & %9.6f  o2 & %9.6f  w10 & %9.6f  o2 &  %9.6f  o2 &) ///
			rspec(&-&&&&&-&)
		}
		else  {
			matlist clustsum , title("Cluster Variability") rowtitle(Statistic) ///
			cspec(& %-10s | %8.2f o4 & %9.6f  o4 & %9.6f  w10 & %9.6f  o4 &) ///
			rspec(&-&&&&&-&)
		}

		/*error for no-beta-g*/
			if `nocv3'==1 {
				disp as error "Some results not available as absorb option is improperly specified."
				disp as text ""
			}

		/*-------------------*/
		/*optional output*/
		/*-------------------*/

		/*gstar*/
		if `gstarind'==1 {

			disp " "
			disp "Effective Number of Clusters"
			disp "-----------------------------"
			disp  "G*(0)  = " %6.3f `gstarzero'
			if `rho_specified' {
				disp   "G*(`rho') = " %6.3f `gstarrho'
			}
			disp  "G*(1)  = " %6.3f `gstarone'
			disp "-----------------------------"
		}
		else if `gstarind'==2 {
			disp " "
			disp "Effective Number of Clusters"
			disp "-----------------------------"
			disp  "G*(0)  = " %6.3f `gstarzero'
			disp "-----------------------------"
			disp as error "G*(rho) and G*(1) are not available."
			disp as error "There are fixed effects at the cluster or subcluster level."
		}

		/*additional summary statistics*/
		if "`svars'"!=""{

			mata: st_matrix("bonus", bonus)
			matrix rownames bonus = "Harmonic Mean" "Harmonic Ratio" "Geometric Mean" "Geometric Ratio" "Quadratic Mean" "Quadratic Ratio"

			if `singclustfull' >= 1 & `colsdrop' >= 1 {
				matrix colnames bonus = Ng Leverage "Partial L." "all beta no g" "kept beta no g"
				matlist bonus, title("Alternative Sample Means and Ratios to Arithmetic Mean") ///
				cspec(& %-10s w15 | %11.3f o2 & %9.6f  w10  & %9.6f  w10 & %9.6f  o2 & %9.6f  o2 &) rspec(&|&&&&&|)
			}
			else {
				matrix colnames bonus = Ng Leverage "Partial L." "beta no g"
				matlist bonus, title("Alternative Sample Means and Ratios to Arithmetic Mean") ///
				cspec(& %-10s w15 | %11.3f o4 & %9.6f  w10  & %9.6f  w10 & %9.6f  o4 &) rspec(&|&&&&&|)
			}





		} /*end of svars*/

		mata: scall = (ng, Leverage, pLeverage, betag)

		/*display std errors for all coefficients*/
		if "`regtable'" != "" {

			/* CV3 and CV3J are covariance matrices. */
			if "`jackknife'" == "" {

				mata fullcv3 = diagonal(cv3)
				mata: fullcv3 = sqrt(fullcv3)

				mata: regresstab = J(rows(fullcv3),6,.)
				mata: regresstab[.,1] =beta
				mata: regresstab[.,2] = fullcv3
				mata: tstats = beta :/ fullcv3
				mata: regresstab[.,3] = tstats

				mata: pvals = J(rows(fullcv3),2,.)

				forvalues i = 1 / `sk' {

					mata: pvals[`i',1] =  ttail(`G'-1,regresstab[`i',3])
					mata: pvals[`i',2] = 1 -  ttail(`G'-1,regresstab[`i',3])
					mata: regresstab[`i',4] = 2*min(pvals[`i',.])
				}

				mata: regresstab[.,5] = regresstab[.,1] - invttail(`G'-1,.025)* regresstab[.,2]
				mata: regresstab[.,6] = regresstab[.,1] + invttail(`G'-1,.025)* regresstab[.,2]

				mata: st_matrix("regresstab",regresstab)

				/*dimensions of this matrix*/
				mata: rtrows = rows(regresstab)
				mata: st_numscalar("rtrows", rtrows)
				mata: st_local("rtrows", strofreal(rtrows))

				matrix colnames regresstab = "Coeff" "Sd. Err." "t-stat" "P value" CI-lower CI-upper

				matrix rownames regresstab = `coefnames'

				if `rtrows' <= 50 {
					local rspecreg: display "-" _dup(`rtrows') "&" "-"
					matlist regresstab, title(Linear Regression -- CV3) rowtitle(`yvar') ///
					cspec(& %-12s w6 | %9.6f w10 & %9.6f & %6.4f w7 & %6.4f w7 & %9.6f w10 & %9.6f w10 &) ///
					rspec(`rspecreg')
				}
				else{
					disp "Linear Regression -- CV3"
					matrix list regresstab
				}

				mata: st_matrix("cv3",cv3)
				return matrix cv3 = cv3

			} /*end of no jackknife if*/

			if "`jackknife'" != "" {

				mata: fullcv3 = sqrt(diagonal(cv3j))
				mata: st_matrix("cv3j",cv3j)
				return matrix cv3j = cv3j

				mata: regresstab = J(rows(fullcv3),6,.)
				mata: regresstab[.,1] =beta
				mata: regresstab[.,2] = fullcv3
				mata: tstats = beta :/ fullcv3
				mata: regresstab[.,3] = tstats

				mata: pvals = J(rows(fullcv3),2,.)

				forvalues i = 1 / `sk' {

					mata: pvals[`i',1] =  ttail(`G'-1,regresstab[`i',3])
					mata: pvals[`i',2] = 1 -  ttail(`G'-1,regresstab[`i',3])
					mata: regresstab[`i',4] = 2*min(pvals[`i',.])
				}

				mata: regresstab[.,5] = regresstab[.,1] - invttail(`G'-1,.025)* regresstab[.,2]
				mata: regresstab[.,6] = regresstab[.,1] + invttail(`G'-1,.025)* regresstab[.,2]

				mata: st_matrix("regresstab",regresstab)

				matrix colnames regresstab = "Coeff" "Sd. Err." "t-stat" "P value" CI-lower CI-upper

				/*dimensions of this matrix*/
				mata: rtrows = rows(regresstab)
				mata: st_numscalar("rtrows", rtrows)
				mata: st_local("rtrows", strofreal(rtrows))

				matrix rownames regresstab = `coefnames'

				if `rtrows' <= 50 {

					local rspecreg: display "-" _dup(`rtrows') "&" "-"

					matlist regresstab, title(Linear Regression -- CV3J) rowtitle(`yvar') ///
					cspec(& %-12s w6 | %9.6f w10 & %9.6f & %6.4f w7 & %6.4f w7 & %9.6f w10 & %9.6f w10 &) ///
					rspec(`rspecreg')
				}
				else{
					disp "Linear Regression -- CV3J"
					matrix list regresstab
				}

			} /*end of yes jackknife if */

			/*singular clusters*/
			if `singclustfull' >= 1 {

					/* CV3 and CV3J are covariance matrices. */
				if "`jackknife'" == "" {

					mata fullcv3drop = diagonal(cv3drop)
					mata: fullcv3drop = sqrt(fullcv3drop)

					mata: regresstabdrop = J(rows(fullcv3drop),6,.)
					mata: regresstabdrop[.,1] =beta
					mata: regresstabdrop[.,2] = fullcv3drop
					mata: tstatsdrop = beta :/ fullcv3drop
					mata: regresstabdrop[.,3] = tstatsdrop

					mata: pvalsdrop = J(rows(fullcv3drop),2,.)

					forvalues i = 1 / `sk' {

						mata: pvalsdrop[`i',1] =  ttail(`gdrop'-1,regresstabdrop[`i',3])
						mata: pvalsdrop[`i',2] = 1 -  ttail(`gdrop'-1,regresstabdrop[`i',3])
						mata: regresstabdrop[`i',4] = 2*min(pvalsdrop[`i',.])
					}

					mata: regresstabdrop[.,5] = regresstabdrop[.,1] - invttail(`gdrop'-1,.025)* regresstabdrop[.,2]
					mata: regresstabdrop[.,6] = regresstabdrop[.,1] + invttail(`gdrop'-1,.025)* regresstabdrop[.,2]

					mata: st_matrix("regresstabdrop",regresstabdrop)

					/*dimensions of this matrix*/
					mata: rtrows = rows(regresstabdrop)
					mata: st_numscalar("rtrows", rtrows)
					mata: st_local("rtrows", strofreal(rtrows))

					matrix colnames regresstabdrop = "Coeff" "Sd. Err." "t-stat" "P value" CI-lower CI-upper

					matrix rownames regresstabdrop = `coefnames'

					if `rtrows' <= 50 {
						local rspecreg: display "-" _dup(`rtrows') "&" "-"

						matlist regresstabdrop, title(Linear Regression -- CV3 -- Omitting Singular Omit-One-Cluster Subsamples) ///
						rowtitle(`yvar') ///
						cspec(& %-12s w6 | %9.6f w10 & %9.6f & %6.4f w7 & %6.4f w7 & %9.6f w10 & %9.6f w10 &) ///
						rspec(`rspecreg')
					}
					else{
						disp "Linear Regression -- CV3 -- Omitting Singular Omit-One-Cluster Subsamples"
						matrix list regresstabdrop
					}

					mata: st_matrix("cv3drop",cv3drop)
					return matrix cv3drop = cv3drop

				} /* end of no jackknife if */

				if "`jackknife'" != "" {

					mata: fullcv3drop = sqrt(diagonal(cv3jdrop))
					mata: st_matrix("cv3jdrop",cv3jdrop)
					return matrix cv3jdrop = cv3jdrop

					mata: regresstabdrop = J(rows(fullcv3drop),6,.)
					mata: regresstabdrop[.,1] =beta
					mata: regresstabdrop[.,2] = fullcv3drop
					mata: tstatsdrop = beta :/ fullcv3drop
					mata: regresstabdrop[.,3] = tstatsdrop

					mata: pvalsdrop = J(rows(fullcv3drop),2,.)

					forvalues i = 1 / `sk' {

						mata: pvalsdrop[`i',1] =  ttail(`gdrop'-1,regresstabdrop[`i',3])
						mata: pvalsdrop[`i',2] = 1 -  ttail(`gdrop'-1,regresstabdrop[`i',3])
						mata: regresstabdrop[`i',4] = 2*min(pvalsdrop[`i',.])
					}

					mata: regresstabdrop[.,5] = regresstabdrop[.,1] - invttail(`gdrop'-1,.025)* regresstabdrop[.,2]
					mata: regresstabdrop[.,6] = regresstabdrop[.,1] + invttail(`gdrop'-1,.025)* regresstabdrop[.,2]

					mata: st_matrix("regresstabdrop",regresstabdrop)
					matrix colnames regresstabdrop = "Coeff" "Sd. Err." "t-stat" "P value" CI-lower CI-upper

					/*dimensions of this matrix*/
					mata: rtrows = rows(regresstabdrop)
					mata: st_numscalar("rtrows", rtrows)
					mata: st_local("rtrows", strofreal(rtrows))

					matrix rownames regresstabdrop = `coefnames'

					if `rtrows' <= 50 {

						local rspecreg: display "-" _dup(`rtrows') "&" "-"

						matlist regresstabdrop, title(Linear Regression -- CV3J -- Omitting Singular Omit-One-Cluster Subsamples) ///
						rowtitle(`yvar') ///
						cspec(& %-12s w6 | %9.6f w10 & %9.6f & %6.4f w7 & %6.4f w7 & %9.6f w10 & %9.6f w10 &) ///
						rspec(`rspecreg')
					}
					else{
						disp "Linear Regression -- CV3J -- Omitting Singular Omit-One-Cluster Subsamples"
						matrix list regresstabdrop
					}

				}/*end of yes jackknife if*/

			} /*end of singclustfull*/

		} /*end of regtable*/

		if "`table'"!="" {

			if `G'< 53 {

				mata: st_matrix("scall", scall)
				matrix rownames scall = `CLUSTERNAMES'
				matrix colnames scall =  Ng Leverage "Partial L." "beta no g"

				local TITLE "Cluster by Cluster Statistics"
				local Gm1 = `G'- 1

				local rspec : display "rspec(&-"_dup(`Gm1') "&" "-)"

				matlist scall, title("`TITLE'") rowtitle("`cluster'") ///
					cspec(& %-10s | %8.0f o4 & %9.6f  o4 & %9.6f  w10 & %9.6f  o4 &) `rspec'

				return matrix scall = scall
			}
			else {

				disp " "
				disp "Cluster by Cluster Statistics"

				disp as error "The number of clusters, `G' exceeds 52, results reported as matrix. "

				disp as text ""

				mata: scall

			}

		} /*end of table*/

		/* Restore the caller's incoming observation order. */
		quietly sort `temp_order'

		/*stuff to return*/
			tempname __sc_cv1 __sc_cv3 __sc_cv3j __sc_legacyflag __sc_rankfail
			mata: st_matrix("`__sc_cv1'",crve1)
			mata: st_matrix("`__sc_cv3'",cv3)
			matrix rownames `__sc_cv1' = `coefnames'
			matrix colnames `__sc_cv1' = `coefnames'
			matrix rownames `__sc_cv3' = `coefnames'
			matrix colnames `__sc_cv3' = `coefnames'
			return matrix cv1 = `__sc_cv1'
			return matrix cv3 = `__sc_cv3'
			if "`jackknife'" != "" {
				mata: st_matrix("`__sc_cv3j'",cv3j)
				matrix rownames `__sc_cv3j' = `coefnames'
				matrix colnames `__sc_cv3j' = `coefnames'
				return matrix cv3j = `__sc_cv3j'
			}
			if "`rankcheck'" != "" {
				mata: st_matrix("`__sc_legacyflag'",badcluster)
				mata: st_matrix("`__sc_rankfail'",rankfail)
				return matrix legacyflag = `__sc_legacyflag'
				return matrix rankfail = `__sc_rankfail'
			}
			return matrix ng = ng
			return matrix lever = Leverage
			return matrix partlev = pLeverage
			return matrix betanog = betag

		/*nograph - default is to draw a graph*/
		if "`nograph'" == "" {
			tempvar ng
			qui getmata `ng' = ng, force
			label variable `ng' "Observations Per Cluster"

			tempvar lever
			qui getmata `lever' = Leverage, force
			label variable `lever' "Leverage"

			tempvar plever
			qui getmata `plever' = pLeverage, force
			label variable `plever' "Partial Leverage"

			tempvar betanog
			qui getmata `betanog' = betag, force
			label variable `betanog' "Omit-One-Cluster Coefficients"

			/* Keep component graphs in memory and use only the G rows. */
			tempname __g_ln __g_pn __g_lb __g_pb

			qui twoway scatter `lever' `ng' in 1/`G', ///
				name(`__g_ln',replace) nodraw scheme(s1mono) ///
				xline(`graph_1') yline(`graph_2')
			qui twoway scatter `plever' `ng' in 1/`G', ///
				name(`__g_pn',replace) nodraw scheme(s1mono) ///
				xline(`graph_1') yline(`graph_3')
			qui twoway scatter `lever' `betanog' in 1/`G', ///
				name(`__g_lb',replace) nodraw scheme(s1mono) ///
				xline(`graph_4') yline(`graph_2')
			qui twoway scatter `plever' `betanog' in 1/`G', ///
				name(`__g_pb',replace) nodraw scheme(s1mono) ///
				xline(`graph_4') yline(`graph_3')

			graph combine `__g_ln' `__g_pn' `__g_lb' `__g_pb', ///
				title("Cluster Specific Statistics For `G' `cluster' Clusters") ///
				scheme(s1mono)

			qui capture graph drop `__g_ln'
			qui capture graph drop `__g_pn'
			qui capture graph drop `__g_lb'
			qui capture graph drop `__g_pb'
		}
		else {
			disp "Graph option suppressed."
		}

		/*additional stuff to return*/
		/*scalars*/

			return local beta =  `BETATEMP'

			/*cv1*/
			return local cv1se = `setemp'
			return local cv1t = `cv1t'
			return local cv1p = `cv1p'
			return local cv1lci = `cv1lci'
			return local cv1uci = `cv1uci'

			/*cv3*/
			return local cv3se = `cv3se'
			return local cv3t = `cv3t'
			return local cv3p = `cv3p'
			return local cv3lci = `cv3lci'
			return local cv3uci = `cv3uci'

			/*cv3drop*/
			if `singclustfull' >= 1 {

				return local cv3sedrop     = `cv3sedrop'
				return local cv3tdrop   = `cv3tdrop'
				return local cv3pdrop   = `cv3pdrop'
				return local cv3lcidrop = `cv3lcidrop'
				return local cv3ucidrop = `cv3ucidrop'

			}

			/*jackknife*/
			if "`jackknife'"!="" {

				return local beta = `BETATEMP'
				return local cv3Jse = `cv3jse'
				return local cv3Jt  = `cv3jt'
				return local cv3Jp = `cv3jp'
				return local cv3Jlci  = `cv3jlci'
				return local cv3Juci = `cv3juci'

				if `singclustfull' >= 1 {

					return local cv3sedrop     = `cv3jsedrop'
					return local cv3tdrop   = `cv3jtdrop'
					return local cv3pdrop   = `cv3jpdrop'
					return local cv3lcidrop = `cv3jlcidrop'
					return local cv3ucidrop = `cv3jucidrop'

				}


			}

			/*gstar*/
			if `gstarind'==1 {
				return local gstarzero = `gstarzero'
				return local gstarone = `gstarone'
				return local gstarrho = `gstarrho'
			}
			else if `gstarind'==2 {
				return local gstarzero = `gstarzero'
			}

end

capture program drop _summclust_twoway
program define _summclust_twoway, rclass
	version 13
	syntax varlist(min=3), cluster(varlist min=2 max=2) [fevar(varlist) absorb(varname) gstar Rho(real -1) SAMple(string) TABle ADDmeans JACKknife REGtable NOGraph RANKcheck]

	if `"`jackknife'"' != "" {
		display as error "jackknife requests the centered CV3J estimator, which remains supported with one-way clustering."
		display as error "Centered two-way CV3J is not studied in MacKinnon, Nielsen, and Webb and is deliberately not implemented; use native two-way CV3 instead."
		exit 198
	}
	local rho_specified = (`rho' != -1)
	if `rho_specified' & !inrange(`rho',0.00001,1) {
		display as error "rho() must be between 0 and 1."
		exit 198
	}
	if `rho' == -1 local rho = 0.1234567

	local yvar  : word 1 of `varlist'
	local xvar  : word 2 of `varlist'
	local nvars : word count `varlist'
	local controls ""
	forvalues j=3/`nvars' {
		local controls `controls' `: word `j' of `varlist''
	}
	local clusterg : word 1 of `cluster'
	local clusterh : word 2 of `cluster'

	tempname TWY TWX0 TWXC TWG TWH TWI TWFE TWAB TWABINFO TWX
	tempvar touse order idg idh idi ycopy
	mark `touse'
	markout `touse' `yvar' `xvar' `controls' `fevar' `absorb' `clusterg' `clusterh', strok
	if `"`sample'"' != "" {
		tempvar samplesel
		quietly generate byte `samplesel' = 0
		quietly replace `samplesel' = 1 if `sample'
		quietly replace `touse' = `touse' * `samplesel'
	}
	quietly count if `touse'
	if r(N)==0 error 2000
	local N = r(N)

	quietly generate long `order' = _n
	quietly egen long `idg' = group(`clusterg') if `touse'
	quietly egen long `idh' = group(`clusterh') if `touse'
	quietly egen long `idi' = group(`clusterg' `clusterh') if `touse'
	quietly summarize `idg' if `touse', meanonly
	local G = r(max)
	quietly summarize `idh' if `touse', meanonly
	local H = r(max)
	quietly summarize `idi' if `touse', meanonly
	local I = r(max)
	if `G' < 2 {
		display as error "Two-way summclust requires at least two nonempty clusters in the first dimension (`clusterg'); found `G'."
		exit 498
	}
	if `H' < 2 {
		display as error "Two-way summclust requires at least two nonempty clusters in the second dimension (`clusterh'); found `H'."
		exit 498
	}
	if `I' < 2 {
		display as error "Two-way summclust requires at least two nonempty G-H intersections; found `I'."
		exit 498
	}

	/* Encode fixed effects so numeric and string identifiers share one path. */
	local fecodes ""
	local fenames ""
	foreach fv of local fevar {
		tempvar fcode
		quietly egen long `fcode' = group(`fv') if `touse'
		local fecodes `fecodes' `fcode'
		quietly summarize `fcode' if `touse', meanonly
		local flast = r(max)
		local ffirst = 1
		if `"`absorb'"' != "" local flast = `flast'-1
		else if `"`fenames'"' != "" local flast = `flast'-1
		if `flast' >= `ffirst' {
			forvalues ell=`ffirst'/`flast' {
				local fenames `fenames' `fv'_`ell'
			}
		}
	}

	/* Demeaning is valid for every delete-cluster component only when */
	/* each absorbed group lies inside one G-H intersection. */
	local abcode ""
	if `"`absorb'"' != "" {
		tempvar abcode sdg sdh
		quietly egen long `abcode' = group(`absorb') if `touse'
		quietly bysort `abcode': egen double `sdg' = sd(`idg') if `touse'
		quietly bysort `abcode': egen double `sdh' = sd(`idh') if `touse'
		quietly replace `sdg' = 0 if `touse' & missing(`sdg')
		quietly replace `sdh' = 0 if `touse' & missing(`sdh')
		quietly summarize `sdg' if `touse', meanonly
		local badg = r(max)
		quietly summarize `sdh' if `touse', meanonly
		if `badg' > 1e-10 | r(max) > 1e-10 {
			quietly sort `order'
			display as error "absorb() groups must be nested within both clustering dimensions for native two-way CV3."
			display as error "Use fevar() for fixed effects that cross either clustering dimension."
			exit 198
		}
	}

	/* Determine dimension-specific availability of G*(1) and G*(rho). */
	local g1ok_g = 1
	local g1ok_h = 1
	local g1ok_i = 1
	if `"`absorb'"' != "" {
		local g1ok_g = 0
		local g1ok_h = 0
		local g1ok_i = 0
	}
	foreach fc of local fecodes {
		foreach d in g h i {
			tempvar sdfe
			quietly bysort `fc': egen double `sdfe' = sd(`id`d'') if `touse'
			quietly replace `sdfe' = 0 if `touse' & missing(`sdfe')
			quietly summarize `sdfe' if `touse', meanonly
			if r(max) <= 1e-10 local g1ok_`d' = 0
		}
	}

	/* Put the estimation sample into Mata exactly once. */
	if `"`absorb'"' != "" quietly sort `abcode' `order'
	else quietly sort `order'
	quietly generate double `ycopy' = `yvar' if `touse'
	quietly putmata `TWY'=`ycopy' `TWX0'=`xvar' `TWXC'=(`controls') `TWG'=`idg' `TWH'=`idh' `TWI'=`idi' if `touse', replace
	if `"`fecodes'"' != "" quietly putmata `TWFE'=(`fecodes') if `touse', replace
	if `"`absorb'"' != "" quietly putmata `TWAB'=`abcode' if `touse', replace

	local coefnames `xvar'
	mata: `TWX' = `TWX0'
	if `"`absorb'"' == "" & `"`fevar'"' == "" {
		mata: `TWX' = `TWX', J(rows(`TWX'),1,1)
		local coefnames `coefnames' _cons
	}
	mata: `TWX' = `TWX', `TWXC'
	local coefnames `coefnames' `controls'
	if `"`fecodes'"' != "" {
		local hasabsorb = (`"`absorb'"' != "")
		mata: `TWX' = `TWX', _summclust_dummies(`TWFE',`hasabsorb')
		local coefnames `coefnames' `fenames'
	}
	if `"`absorb'"' != "" {
		mata: `TWABINFO' = panelsetup(`TWAB',1)
		mata: `TWY' = _summclust_within(`TWY',`TWABINFO')
		mata: `TWX' = _summclust_within(`TWX',`TWABINFO')
		mata: st_local("tw_nabsorb",strofreal(rows(`TWABINFO')))
	}
	else local tw_nabsorb = 0

	mata: st_local("tw_k",strofreal(cols(`TWX')))
	mata: st_local("tw_rank",strofreal(rank(`TWX')))
	local coeffnum = `tw_rank' + `tw_nabsorb'
	if `N' <= `coeffnum' {
		quietly sort `order'
		capture mata: mata drop `TWY'
		capture mata: mata drop `TWX0'
		capture mata: mata drop `TWXC'
		capture mata: mata drop `TWG'
		capture mata: mata drop `TWH'
		capture mata: mata drop `TWI'
		capture mata: mata drop `TWFE'
		capture mata: mata drop `TWAB'
		capture mata: mata drop `TWABINFO'
		capture mata: mata drop `TWX'
		display as error "Residual degrees of freedom are not positive after fixed effects."
		exit 2001
	}

	local dorank = (`"`rankcheck'"' != "")
	tempname b V1g V1h V1i V1 V3g V3h V3i V3 dsum regsum tabg tabh
	tempname ngg ngh ngi lvg lvh lvi plg plh pli bog boh boi rfg rfh rfi
	tempname V1copy V3copy V1gcopy V1hcopy V1icopy V3gcopy V3hcopy V3icopy
	tempname ngicopy lvicopy plicopy boicopy rficopy
	capture noisily mata: _summclust_twoway_core(`TWX',`TWY',`TWG',`TWH',`TWI',`coeffnum',`dorank',`rho',("`b'","`V1g'","`V1h'","`V1i'","`V1'","`V3g'","`V3h'","`V3i'","`V3'","`dsum'","`regsum'","`ngg'","`ngh'","`ngi'","`lvg'","`lvh'","`lvi'","`plg'","`plh'","`pli'","`bog'","`boh'","`boi'","`rfg'","`rfh'","`rfi'"))
	local tw_core_rc = _rc
	quietly sort `order'
	capture mata: mata drop `TWY'
	capture mata: mata drop `TWX0'
	capture mata: mata drop `TWXC'
	capture mata: mata drop `TWG'
	capture mata: mata drop `TWH'
	capture mata: mata drop `TWI'
	capture mata: mata drop `TWFE'
	capture mata: mata drop `TWAB'
	capture mata: mata drop `TWABINFO'
	capture mata: mata drop `TWX'
	if `tw_core_rc' exit `tw_core_rc'

	foreach M in b V1g V1h V1i V1 V3g V3h V3i V3 {
		matrix colnames ``M'' = `coefnames'
	}
	foreach M in V1g V1h V1i V1 V3g V3h V3i V3 {
		matrix rownames ``M'' = `coefnames'
	}
	matrix `V1copy' = `V1'
	matrix `V3copy' = `V3'
	matrix `V1gcopy' = `V1g'
	matrix `V1hcopy' = `V1h'
	matrix `V1icopy' = `V1i'
	matrix `V3gcopy' = `V3g'
	matrix `V3hcopy' = `V3h'
	matrix `V3icopy' = `V3i'
	matrix `ngicopy' = `ngi'
	matrix `lvicopy' = `lvi'
	matrix `plicopy' = `pli'
	matrix `boicopy' = `boi'
	matrix `rficopy' = `rfi'
	matrix rownames `dsum' = G H I
	matrix colnames `dsum' = cv_N cv_leverage cv_partlev cv_betanog clusters Gstar0
	matrix rownames `regsum' = `coefnames'
	matrix colnames `regsum' = beta CV1max_se CV3max_se
	if `"`table'"' != "" {
		matrix `tabg' = `ngg', `lvg', `plg', `bog'
		matrix `tabh' = `ngh', `lvh', `plh', `boh'
		matrix colnames `tabg' = N leverage partial_leverage beta_minus
		matrix colnames `tabh' = N leverage partial_leverage beta_minus
	}

	local df = min(`G',`H')-1
	local cv1_conventional_var = el(`V1',1,1)
	local cv3_conventional_var = el(`V3',1,1)
	local cv1_conventional_se = cond(`cv1_conventional_var'>=0,sqrt(`cv1_conventional_var'),.)
	local cv3_conventional_se = cond(`cv3_conventional_var'>=0,sqrt(`cv3_conventional_var'),.)
	local cv1t = `tw_beta'/`tw_cv1se'
	local cv3t = `tw_beta'/`tw_cv3se'
	local cv1p = 2*ttail(`df',abs(`cv1t'))
	local cv3p = 2*ttail(`df',abs(`cv3t'))
	local crit = invttail(`df',.025)
	local cv1lci = `tw_beta'-`crit'*`tw_cv1se'
	local cv1uci = `tw_beta'+`crit'*`tw_cv1se'
	local cv3lci = `tw_beta'-`crit'*`tw_cv3se'
	local cv3uci = `tw_beta'+`crit'*`tw_cv3se'

	display as text _newline "Native two-way cluster variance estimation"
	display as text "Clustering: G=" as result "`clusterg'" as text " (" as result %8.0f `G' as text "), H=" as result "`clusterh'" as text " (" as result %8.0f `H' as text "), intersections=" as result %8.0f `I'
	display as text "Inference df = min(G,H) - 1 = " as result %8.0f `df'
	if `"`clusterg'"' == `"`clusterh'"' display as text "note: the clustering dimensions are identical; conventional two-way results reduce to one-way clustering."
	if min(`G',`H') == 2 display as error "warning: only two clusters occur in the smaller dimension; inference uses 1 degree of freedom and may be unreliable."
	display as text ""
	display as text "Focal inference using MacKinnon-Nielsen-Webb max standard errors"
	display as text %9s "s.e." " |" %11s "Coefficient" %11s "Std. err." %9s "t" %9s "P>|t|" %11s "[95% conf." %11s "interval]"
	display as text "----------+-----------------------------------------------------------------"
	display as text %9s "CV1 max" " |" as result %11.6f `tw_beta' %11.6f `tw_cv1se' %9.4f `cv1t' %9.4f `cv1p' %11.6f `cv1lci' %11.6f `cv1uci'
	display as text %9s "CV3 max" " |" as result %11.6f `tw_beta' %11.6f `tw_cv3se' %9.4f `cv3t' %9.4f `cv3p' %11.6f `cv3lci' %11.6f `cv3uci'
	display as text "----------+-----------------------------------------------------------------"
	display as text "Max selectors: CV1=" as result "`tw_cv1which'" as text "; CV3=" as result "`tw_cv3which'"
	display as text "Conventional focal CV1: variance=" as result %11.6f `cv1_conventional_var' as text "  s.e.=" as result %11.6f `cv1_conventional_se'
	display as text "Conventional focal CV3: variance=" as result %11.6f `cv3_conventional_var' as text "  s.e.=" as result %11.6f `cv3_conventional_se'
	if `cv1_conventional_var' < 0 display as error "warning: the conventional two-way CV1 focal variance is negative; its conventional standard error is undefined. MNW max inference remains available."
	if `cv3_conventional_var' < 0 display as error "warning: the conventional two-way CV3 focal variance is negative; its conventional standard error is undefined. MNW max inference remains available."
	display as text "Conventional matrices use V_G + V_H - V_intersection and are not used for the max-se inference above."
	matlist `dsum', names(rows) format(%8.3f) title("Dimension-specific diagnostics")
	if `"`table'"' != "" {
		if `G' <= 52 matlist `tabg', names(rows) format(%10.4f) title("Cluster diagnostics: G = `clusterg' (sorted cluster order)")
		else display as text "note: G-level table suppressed because `clusterg' has more than 52 clusters; use r(table_g)."
		if `H' <= 52 matlist `tabh', names(rows) format(%10.4f) title("Cluster diagnostics: H = `clusterh' (sorted cluster order)")
		else display as text "note: H-level table suppressed because `clusterh' has more than 52 clusters; use r(table_h)."
		display as text "note: intersection-level detail is not printed; use the *_intersection matrices in r()."
	}
	if `"`regtable'"' != "" matlist `regsum', names(rows) format(%11.6f) title("Coefficient-wise MNW max standard errors; df=min(G,H)-1 (not a covariance matrix)")
	if `"`addmeans'"' != "" display as text "note: addmeans has no additional effect for the two-way dimension summary."
	if `"`nograph'"' == "" display as text "note: two-way diagnostic graphs are not drawn; use the dimension-specific diagnostic matrices in r()."
	if `"`rankcheck'"' != "" & (`tw_rankfail_g'+`tw_rankfail_h'+`tw_rankfail_i')>0 {
		display as error "warning: rankcheck found rank-deficient omit-one subsamples: G=`tw_rankfail_g', H=`tw_rankfail_h', I=`tw_rankfail_i'."
	}

	return matrix b = `b'
	return matrix cv1 = `V1'
	return matrix cv3 = `V3'
	return matrix cv1_conventional = `V1copy'
	return matrix cv3_conventional = `V3copy'
	return matrix cv1_g = `V1gcopy'
	return matrix cv1_h = `V1hcopy'
	return matrix cv1_intersection = `V1icopy'
	return matrix cv3_g = `V3gcopy'
	return matrix cv3_h = `V3hcopy'
	return matrix cv3_intersection = `V3icopy'
	return matrix V_cv1_g = `V1g'
	return matrix V_cv1_h = `V1h'
	return matrix V_cv1_i = `V1i'
	return matrix V_cv3_g = `V3g'
	return matrix V_cv3_h = `V3h'
	return matrix V_cv3_i = `V3i'
	return matrix dimsummary = `dsum'
	return matrix ng_g = `ngg'
	return matrix ng_h = `ngh'
	return matrix ng_i = `ngi'
	return matrix lever_g = `lvg'
	return matrix lever_h = `lvh'
	return matrix lever_i = `lvi'
	return matrix partlev_g = `plg'
	return matrix partlev_h = `plh'
	return matrix partlev_i = `pli'
	return matrix betanog_g = `bog'
	return matrix betanog_h = `boh'
	return matrix betanog_i = `boi'
	return matrix rankfail_g = `rfg'
	return matrix rankfail_h = `rfh'
	return matrix rankfail_i = `rfi'
	return matrix ng_intersection = `ngicopy'
	return matrix lever_intersection = `lvicopy'
	return matrix partlev_intersection = `plicopy'
	return matrix betanog_intersection = `boicopy'
	return matrix rankfail_intersection = `rficopy'
	if `"`table'"' != "" {
		return matrix table_g = `tabg'
		return matrix table_h = `tabh'
	}
	return scalar N = `N'
	return scalar k = `tw_k'
	return scalar rank = `tw_rank'
	return scalar G = `G'
	return scalar H = `H'
	return scalar I = `I'
	return scalar N_clust_g = `G'
	return scalar N_clust_h = `H'
	return scalar N_clust_intersection = `I'
	return scalar df_r = `df'
	return scalar cv1_conventional_var = `cv1_conventional_var'
	return scalar cv3_conventional_var = `cv3_conventional_var'
	return scalar cv1_conventional_se = `cv1_conventional_se'
	return scalar cv3_conventional_se = `cv3_conventional_se'
	return scalar cv1max_ties = `tw_cv1ties'
	return scalar cv3max_ties = `tw_cv3ties'
	return scalar max_tie_tolerance = 1e-10
	return scalar cv1_max_var = `tw_cv1var'
	return scalar cv3_max_var = `tw_cv3var'
	return scalar cv1_max_se = `tw_cv1se'
	return scalar cv3_max_se = `tw_cv3se'
	return scalar cv1var_max = `tw_cv1var'
	return scalar cv3var_max = `tw_cv3var'
	return scalar cv1se_max = `tw_cv1se'
	return scalar cv3se_max = `tw_cv3se'
	return scalar gstarzero_g = `tw_gstar0_g'
	return scalar gstarzero_h = `tw_gstar0_h'
	return scalar gstarzero_i = `tw_gstar0_i'
	return scalar gstarzero_intersection = `tw_gstar0_i'
	if (`"`gstar'"' != "" | `rho_specified') & `g1ok_g' {
		return scalar gstarone_g = `tw_gstar1_g'
		return scalar gstarrho_g = `tw_gstarrho_g'
	}
	if (`"`gstar'"' != "" | `rho_specified') & `g1ok_h' {
		return scalar gstarone_h = `tw_gstar1_h'
		return scalar gstarrho_h = `tw_gstarrho_h'
	}
	if (`"`gstar'"' != "" | `rho_specified') & `g1ok_i' {
		return scalar gstarone_i = `tw_gstar1_i'
		return scalar gstarrho_i = `tw_gstarrho_i'
		return scalar gstarone_intersection = `tw_gstar1_i'
		return scalar gstarrho_intersection = `tw_gstarrho_i'
	}
	return local beta `tw_beta'
	return local cv1se `tw_cv1se'
	return local cv1t `cv1t'
	return local cv1p `cv1p'
	return local cv1lci `cv1lci'
	return local cv1uci `cv1uci'
	return local cv3se `tw_cv3se'
	return local cv3t `cv3t'
	return local cv3p `cv3p'
	return local cv3lci `cv3lci'
	return local cv3uci `cv3uci'
	return local cv1max_component `tw_cv1which'
	return local cv3max_component `tw_cv3which'
	return local cluster_g `clusterg'
	return local cluster_h `clusterh'
	return local vce "two-way max-se"
end

capture mata: mata drop _summclust_within()
capture mata: mata drop _summclust_dummies()
capture mata: mata drop _summclust_focal_identified()
capture mata: mata drop _summclust_omitsolve()
capture mata: mata drop _summclust_core()
capture mata: mata drop _summclust_gstar_all()
capture mata: mata drop _summclust_gstar_zero()
capture mata: mata drop _summclust_tw_cv()
capture mata: mata drop _summclust_tw_gstar()
capture mata: mata drop _summclust_tw_component()
capture mata: mata drop _summclust_twoway_core()

mata

real matrix _summclust_dummies(real matrix FE, real scalar hasabsorb)
{
	real scalar q, l, nkeep, total, outcol
	real matrix levels, D

	total = 0
	for (q=1; q<=cols(FE); q++) {
		levels = uniqrows(sort(FE[.,q],1))
		total = total + rows(levels) - (hasabsorb | q>1)
	}

	D = J(rows(FE),total,0)
	outcol = 0
	for (q=1; q<=cols(FE); q++) {
		levels = uniqrows(sort(FE[.,q],1))
		nkeep = rows(levels) - (hasabsorb | q>1)
		for (l=1; l<=nkeep; l++) {
			outcol = outcol + 1
			D[.,outcol] = (FE[.,q]:==levels[l])
		}
	}
	return(D)
}

real matrix _summclust_within(real matrix A, real matrix pinfo)
{
	real matrix means, rows_h
	real scalar h

	means = panelsum(A,pinfo) :/ (pinfo[.,2]-pinfo[.,1]:+1)
	for (h=1; h<=rows(pinfo); h++) {
		rows_h = pinfo[h,1]::pinfo[h,2]
		A[rows_h,.] = A[rows_h,.] - J(rows(rows_h),1,1)*means[h,.]
	}
	return(A)
}

real scalar _summclust_focal_identified(real matrix A)
{
    real matrix Ai, target
    Ai = invsym(A)
    if (hasmissing(Ai)) return(0)
    target = J(1,cols(A),0)
    target[1,1] = 1
    return(max(abs(Ai[1,.]*A-target)) < 1e-7)
}

real matrix _summclust_omitsolve(
	real matrix A,
	real matrix b,
	real scalar forceginv,
	real scalar checkcond)
{
	real matrix result

	if (forceginv) return(invsym(A)*b)
	if (checkcond) {
		if (cond(A)>1e8) return(invsym(A)*b)
	}
	result = cholsolve(A,b)
	if (hasmissing(result)) result = invsym(A)*b
	return(result)
}

void _summclust_core(
	real matrix X,
	real matrix Y,
	real matrix info,
	real scalar gfact,
	real scalar dojack,
	real scalar dogstar,
	real scalar dorankcheck,
	real scalar fullrank)
{
	real scalar Gm, Km, g, betatemp, cv1se, cv3se, cv3jse, forceginv
	real scalar gdrop, colsdrop, cv3sedrop, cv3jsedrop, xttxtt
	real matrix C, CtC, Ctx, cgamma, xtt, xw
	real matrix Xg, Yg, Ug, XgXg, XgYg, Zg, left, right, beta_g
	real matrix betadiff, betagbar, betadiffj
	real matrix betadiffdrop, betagdropbar, betadiffjdrop
	external real matrix XpX, XpY, XpXi, beta, U, ZZp
	external real matrix betags, Leverage, badcluster, rankfail, pLeverage
	external real matrix gamzero_all, gamone_all, crve1
	external real matrix cv3, cv3j, betadropgs, cv3drop, cv3jdrop
	external real matrix betag, betadropg

	Gm = rows(info)
	Km = cols(X)
	XpX  = quadcross(X,X)
	XpY  = quadcross(X,Y)
	forceginv = (cond(XpX)>1e7)
	if (forceginv) {
		/* Match 4.22's cluster-summed accumulation when rounding can */
		/* materially affect a generalized inverse. */
		XpX = J(Km,Km,0)
		XpY = J(Km,1,0)
		for (g=1; g<=Gm; g++) {
			Xg = panelsubmatrix(X,g,info)
			Yg = panelsubmatrix(Y,g,info)
			XpX = XpX+quadcross(Xg,Xg)
			XpY = XpY+quadcross(Xg,Yg)
		}
	}
	XpXi = invsym(XpX)
	beta = XpXi*XpY
	U    = Y-X*beta

	ZZp        = J(Km,Km,0)
	betags     = J(Km,Gm,.)
	Leverage   = J(Gm,1,.)
	badcluster = J(Gm,1,0)
	rankfail   = J(Gm,1,.)

	/* Residualize the focal regressor on every other design column. */
	C = X[.,2..Km]
	CtC = quadcross(C,C)
	Ctx = quadcross(C,X[.,1])
	cgamma = cholsolve(CtC,Ctx)
	if (hasmissing(cgamma)) cgamma = invsym(CtC)*Ctx
	xtt = X[.,1]-C*cgamma
	xttxtt = cross(xtt,xtt)
	pLeverage = panelsum(xtt:^2,info):/xttxtt

	if (dogstar) {
		xw = X*XpXi[.,1]
		gamzero_all = panelsum(xw:^2,info)
		gamone_all  = panelsum(xw,info):^2
	}
	else {
		gamzero_all = J(0,1,.)
		gamone_all  = J(0,1,.)
	}

	for (g=1; g<=Gm; g++) {
		Xg = panelsubmatrix(X,g,info)
		Yg = panelsubmatrix(Y,g,info)
		Ug = panelsubmatrix(U,g,info)
		XgXg = quadcross(Xg,Xg)
		XgYg = quadcross(Xg,Yg)
		Zg   = Xg'*Ug
		ZZp  = ZZp+Zg*Zg'

		left = XpX-XgXg
		right = XpY-XgYg
		beta_g = _summclust_omitsolve(left,right,forceginv,(Km<=50))
		betags[.,g] = beta_g
		Leverage[g,1] = trace(XgXg*XpXi)

		/* Compatibility flag: this is not a mathematical rank test. */
		badcluster[g,1] = (sum(abs(beta_g):<=1e-8)>0)
		if (dorankcheck) rankfail[g,1] = (rank(left)<fullrank)
	}

	crve1 = gfact:*(XpXi*ZZp*XpXi)
	betatemp = beta[1,1]
	cv1se = sqrt(crve1[1,1])

	betadiff = beta*J(1,Gm,1)-betags
	cv3 = ((Gm-1)/Gm):*(betadiff*betadiff')
	cv3se = sqrt(cv3[1,1])

	if (dojack) {
		betagbar = mean(betags')
		betadiffj = betags'-J(Gm,1,1)*betagbar
		cv3j = ((Gm-1)/Gm):*quadcross(betadiffj,betadiffj)
		cv3jse = sqrt(cv3j[1,1])
	}
	else {
		cv3j = J(Km,Km,.)
		cv3jse = .
	}

	betadropgs = select(betags,(badcluster:==0)')
	gdrop = cols(betadropgs)
	colsdrop = gdrop
	if (gdrop>0) {
		betadropg = betadropgs[1,.]'
		betadiffdrop = beta*J(1,gdrop,1)-betadropgs
		cv3drop = ((gdrop-1)/gdrop):*(betadiffdrop*betadiffdrop')
		cv3sedrop = sqrt(cv3drop[1,1])
		if (dojack) {
			betagdropbar = mean(betadropgs')
			betadiffjdrop = betadropgs'-J(gdrop,1,1)*betagdropbar
			cv3jdrop = ((gdrop-1)/gdrop):*quadcross(betadiffjdrop,betadiffjdrop)
			cv3jsedrop = sqrt(cv3jdrop[1,1])
		}
		else {
			cv3jdrop = J(Km,Km,.)
			cv3jsedrop = .
		}
	}
	else {
		betadropg  = J(0,1,.)
		cv3drop    = J(Km,Km,.)
		cv3jdrop   = J(Km,Km,.)
		cv3sedrop  = .
		cv3jsedrop = .
	}

	betag = betags[1,.]'
	st_numscalar("BETATEMP",betatemp)
	st_numscalar("setemp",cv1se)
	st_numscalar("cv3se",cv3se)
	st_numscalar("cv3jse",cv3jse)
	st_numscalar("gdrop",gdrop)
	st_numscalar("colsdrop",colsdrop)
	st_numscalar("cv3sedrop",cv3sedrop)
	st_numscalar("cv3jsedrop",cv3jsedrop)
	st_local("BETATEMP",strofreal(betatemp))
	st_local("setemp",strofreal(cv1se))
	st_local("cv3se",strofreal(cv3se))
	st_local("cv3jse",strofreal(cv3jse))
	st_local("singclustfull",strofreal(sum(badcluster)))
	st_local("gdrop",strofreal(gdrop))
	st_local("colsdrop",strofreal(colsdrop))
	st_local("cv3sedrop",strofreal(cv3sedrop))
	st_local("cv3jsedrop",strofreal(cv3jsedrop))
}

void _summclust_gstar_all(real scalar rho)
{
	real scalar Gm, gambarz, gambaro, gambarr
	real scalar gammazero, gammaone, gammarho
	real scalar gstarzero, gstarone, gstarrho
	real matrix gamrho_all
	external real matrix gamzero_all, gamone_all

	Gm = rows(gamzero_all)
	gamrho_all = rho:*gamone_all+(1-rho):*gamzero_all
	gambarz = mean(gamzero_all)
	gambaro = mean(gamone_all)
	gambarr = mean(gamrho_all)
	gammazero = mean(((gamzero_all:-gambarz):/gambarz):^2)
	gammaone  = mean(((gamone_all:-gambaro):/gambaro):^2)
	gammarho  = mean(((gamrho_all:-gambarr):/gambarr):^2)
	gstarzero = Gm/(1+gammazero)
	gstarone  = Gm/(1+gammaone)
	gstarrho  = Gm/(1+gammarho)
	st_local("gstarzero",strofreal(gstarzero))
	st_local("gstarone",strofreal(gstarone))
	st_local("gstarrho",strofreal(gstarrho))
}

void _summclust_gstar_zero()
{
	real scalar Gm, gambarz, gammazero, gstarzero
	external real matrix gamzero_all

	Gm = rows(gamzero_all)
	gambarz = mean(gamzero_all)
	gammazero = mean(((gamzero_all:-gambarz):/gambarz):^2)
	gstarzero = Gm/(1+gammazero)
	st_local("gstarzero",strofreal(gstarzero))
}

struct _summclust_tw_result {
	real scalar J
	real matrix cv1
	real matrix cv3
	real matrix ng
	real matrix leverage
	real matrix partlev
	real matrix betanog
	real matrix rankfail
	real matrix gamzero
	real matrix gamone
}

real scalar _summclust_tw_cv(real matrix z)
{
	real scalar m
	if (rows(z)<2) return(.)
	m = mean(z)
	if (missing(m) | abs(m)<1e-14) return(.)
	return(sqrt(sum((z:-m):^2)/((rows(z)-1)*m^2)))
}

real scalar _summclust_tw_gstar(real matrix gam)
{
	real scalar m, cv2
	m = mean(gam)
	if (missing(m) | abs(m)<1e-14) return(.)
	cv2 = mean(((gam:-m):/m):^2)
	return(rows(gam)/(1+cv2))
}

struct _summclust_tw_result scalar _summclust_tw_component(
	real matrix X,
	real matrix Y,
	real matrix U,
	real matrix cid,
	real matrix XpX,
	real matrix XpY,
	real matrix XpXi,
	real matrix beta,
	real matrix xtt,
	real scalar xttxtt,
	real matrix xw,
	real scalar coeffnum,
	real scalar dorank,
	real scalar fullrank,
	real scalar forceginv,
	string scalar dimension)
{
	struct _summclust_tw_result scalar out
	real scalar j, K, gfact
	real matrix ord, cs, Xs, Ys, Us, xts, xws, info
	real matrix Xj, Yj, Uj, XjXj, XjYj, Zj, ZZp, left, bminus, bdiff

	ord = order(cid,1)
	cs = cid[ord,.]
	Xs = X[ord,.]
	Ys = Y[ord,.]
	Us = U[ord,.]
	xts = xtt[ord,.]
	xws = xw[ord,.]
	info = panelsetup(cs,1)
	out.J = rows(info)
	K = cols(X)
	gfact = out.J*(rows(X)-1)/((out.J-1)*(rows(X)-coeffnum))
	out.ng = info[.,2]-info[.,1]:+1
	out.leverage = J(out.J,1,.)
	out.partlev = panelsum(xts:^2,info):/xttxtt
	out.betanog = J(out.J,1,.)
	out.rankfail = J(out.J,1,.)
	out.gamzero = panelsum(xws:^2,info)
	out.gamone = panelsum(xws,info):^2
	out.cv3 = J(K,K,0)
	ZZp = J(K,K,0)

	for (j=1; j<=out.J; j++) {
		Xj = panelsubmatrix(Xs,j,info)
		Yj = panelsubmatrix(Ys,j,info)
		Uj = panelsubmatrix(Us,j,info)
		XjXj = quadcross(Xj,Xj)
		XjYj = quadcross(Xj,Yj)
		Zj = Xj'*Uj
		ZZp = ZZp + Zj*Zj'
		left = XpX-XjXj
        if (!_summclust_focal_identified(left)) {
            errprintf("Focal coefficient is not numerically identified after deleting %s cluster %g. CV3 is undefined.\n",dimension,j)
            _error(498)
        }
		bminus = _summclust_omitsolve(left,XpY-XjYj,forceginv,(K<=50))
		out.betanog[j,1] = bminus[1,1]
		bdiff = beta-bminus
		out.cv3 = out.cv3 + bdiff*bdiff'
		out.leverage[j,1] = trace(XjXj*XpXi)
		if (dorank) out.rankfail[j,1] = (rank(left)<fullrank)
	}
	out.cv1 = gfact:*(XpXi*ZZp*XpXi)
	out.cv3 = ((out.J-1)/out.J):*out.cv3
	return(out)
}

void _summclust_twoway_core(
	real matrix X,
	real matrix Y,
	real matrix cidg,
	real matrix cidh,
	real matrix cidi,
	real scalar coeffnum,
	real scalar dorank,
	real scalar rho,
	string rowvector names)
{
	struct _summclust_tw_result scalar g, h, i
	real scalar K, q, best, forceginv, cv1var, cv3var, cv1se, cv3se
	real scalar g0, h0, i0, g1, h1, i1, gr, hr, ir
	real matrix XpX, XpY, XpXi, beta, U, C, gamma, xtt, xw
	real matrix V1, V3, dsum, regsum, se1, se3
	string scalar which1, which3

	K = cols(X)
	XpX = quadcross(X,X)
	XpY = quadcross(X,Y)
	forceginv = (cond(XpX)>1e7)
    if (!_summclust_focal_identified(XpX)) {
        errprintf("Focal coefficient is not numerically identified in the full model.\n")
        _error(498)
    }
	XpXi = invsym(XpX)
	beta = XpXi*XpY
	U = Y-X*beta
	C = X[.,2..K]
	gamma = cholsolve(quadcross(C,C),quadcross(C,X[.,1]))
	if (hasmissing(gamma)) gamma = invsym(quadcross(C,C))*quadcross(C,X[.,1])
	xtt = X[.,1]-C*gamma
	xw = X*XpXi[.,1]

	g = _summclust_tw_component(X,Y,U,cidg,XpX,XpY,XpXi,beta,xtt,cross(xtt,xtt),xw,coeffnum,dorank,rank(X),forceginv,"G")
	h = _summclust_tw_component(X,Y,U,cidh,XpX,XpY,XpXi,beta,xtt,cross(xtt,xtt),xw,coeffnum,dorank,rank(X),forceginv,"H")
	i = _summclust_tw_component(X,Y,U,cidi,XpX,XpY,XpXi,beta,xtt,cross(xtt,xtt),xw,coeffnum,dorank,rank(X),forceginv,"intersection")

    if (hasmissing(g.cv1) | hasmissing(h.cv1) | hasmissing(i.cv1) |
        hasmissing(g.cv3) | hasmissing(h.cv3) | hasmissing(i.cv3)) {
        errprintf("A two-way variance component is missing; combined inference is undefined.\n")
        _error(498)
    }
	V1 = g.cv1+h.cv1-i.cv1
	V3 = g.cv3+h.cv3-i.cv3
	which1 = "two-way"
	cv1var = V1[1,1]
	if (!missing(g.cv1[1,1]) & (missing(cv1var) | g.cv1[1,1]>cv1var)) {
		cv1var = g.cv1[1,1]
		which1 = "G"
	}
	if (!missing(h.cv1[1,1]) & (missing(cv1var) | h.cv1[1,1]>cv1var)) {
		cv1var = h.cv1[1,1]
		which1 = "H"
	}
	which3 = "two-way"
	cv3var = V3[1,1]
	if (!missing(g.cv3[1,1]) & (missing(cv3var) | g.cv3[1,1]>cv3var)) {
		cv3var = g.cv3[1,1]
		which3 = "G"
	}
	if (!missing(h.cv3[1,1]) & (missing(cv3var) | h.cv3[1,1]>cv3var)) {
		cv3var = h.cv3[1,1]
		which3 = "H"
	}
	cv1se = sqrt(cv1var)
	cv3se = sqrt(cv3var)

	g0 = _summclust_tw_gstar(g.gamzero)
	h0 = _summclust_tw_gstar(h.gamzero)
	i0 = _summclust_tw_gstar(i.gamzero)
	g1 = _summclust_tw_gstar(g.gamone)
	h1 = _summclust_tw_gstar(h.gamone)
	i1 = _summclust_tw_gstar(i.gamone)
	gr = _summclust_tw_gstar(rho:*g.gamone+(1-rho):*g.gamzero)
	hr = _summclust_tw_gstar(rho:*h.gamone+(1-rho):*h.gamzero)
	ir = _summclust_tw_gstar(rho:*i.gamone+(1-rho):*i.gamzero)

	dsum = ((_summclust_tw_cv(g.ng),_summclust_tw_cv(g.leverage),_summclust_tw_cv(g.partlev),_summclust_tw_cv(g.betanog),g.J,g0)
		\ (_summclust_tw_cv(h.ng),_summclust_tw_cv(h.leverage),_summclust_tw_cv(h.partlev),_summclust_tw_cv(h.betanog),h.J,h0)
		\ (_summclust_tw_cv(i.ng),_summclust_tw_cv(i.leverage),_summclust_tw_cv(i.partlev),_summclust_tw_cv(i.betanog),i.J,i0))
	se1 = J(K,1,.)
	se3 = J(K,1,.)
	for (q=1; q<=K; q++) {
		best = V1[q,q]
		if (!missing(g.cv1[q,q]) & (missing(best) | g.cv1[q,q]>best)) best=g.cv1[q,q]
		if (!missing(h.cv1[q,q]) & (missing(best) | h.cv1[q,q]>best)) best=h.cv1[q,q]
		if (!missing(best) & best>=0) se1[q,1]=sqrt(best)
		best = V3[q,q]
		if (!missing(g.cv3[q,q]) & (missing(best) | g.cv3[q,q]>best)) best=g.cv3[q,q]
		if (!missing(h.cv3[q,q]) & (missing(best) | h.cv3[q,q]>best)) best=h.cv3[q,q]
		if (!missing(best) & best>=0) se3[q,1]=sqrt(best)
	}
	regsum = (beta,se1,se3)

	st_matrix(names[1],beta')
	st_matrix(names[2],g.cv1)
	st_matrix(names[3],h.cv1)
	st_matrix(names[4],i.cv1)
	st_matrix(names[5],V1)
	st_matrix(names[6],g.cv3)
	st_matrix(names[7],h.cv3)
	st_matrix(names[8],i.cv3)
	st_matrix(names[9],V3)
	st_matrix(names[10],dsum)
	st_matrix(names[11],regsum)
	st_matrix(names[12],g.ng)
	st_matrix(names[13],h.ng)
	st_matrix(names[14],i.ng)
	st_matrix(names[15],g.leverage)
	st_matrix(names[16],h.leverage)
	st_matrix(names[17],i.leverage)
	st_matrix(names[18],g.partlev)
	st_matrix(names[19],h.partlev)
	st_matrix(names[20],i.partlev)
	st_matrix(names[21],g.betanog)
	st_matrix(names[22],h.betanog)
	st_matrix(names[23],i.betanog)
	st_matrix(names[24],g.rankfail)
	st_matrix(names[25],h.rankfail)
	st_matrix(names[26],i.rankfail)

	st_local("tw_beta",strofreal(beta[1,1],"%21.17g"))
	st_local("tw_cv1var",strofreal(cv1var,"%21.17g"))
	st_local("tw_cv3var",strofreal(cv3var,"%21.17g"))
	st_local("tw_cv1se",strofreal(cv1se,"%21.17g"))
	st_local("tw_cv3se",strofreal(cv3se,"%21.17g"))
    st_local("tw_cv1ties",strofreal((abs(V1[1,1]-cv1var)<=1e-10*max((abs(cv1var),1e-300)))+2*(abs(g.cv1[1,1]-cv1var)<=1e-10*max((abs(cv1var),1e-300)))+4*(abs(h.cv1[1,1]-cv1var)<=1e-10*max((abs(cv1var),1e-300)))))
    st_local("tw_cv3ties",strofreal((abs(V3[1,1]-cv3var)<=1e-10*max((abs(cv3var),1e-300)))+2*(abs(g.cv3[1,1]-cv3var)<=1e-10*max((abs(cv3var),1e-300)))+4*(abs(h.cv3[1,1]-cv3var)<=1e-10*max((abs(cv3var),1e-300)))))
	st_local("tw_cv1which",which1)
	st_local("tw_cv3which",which3)
	st_local("tw_gstar0_g",strofreal(g0,"%21.17g"))
	st_local("tw_gstar0_h",strofreal(h0,"%21.17g"))
	st_local("tw_gstar0_i",strofreal(i0,"%21.17g"))
	st_local("tw_gstar1_g",strofreal(g1,"%21.17g"))
	st_local("tw_gstar1_h",strofreal(h1,"%21.17g"))
	st_local("tw_gstar1_i",strofreal(i1,"%21.17g"))
	st_local("tw_gstarrho_g",strofreal(gr,"%21.17g"))
	st_local("tw_gstarrho_h",strofreal(hr,"%21.17g"))
	st_local("tw_gstarrho_i",strofreal(ir,"%21.17g"))
	st_local("tw_rankfail_g",strofreal(dorank ? sum(g.rankfail) : 0))
	st_local("tw_rankfail_h",strofreal(dorank ? sum(h.rankfail) : 0))
	st_local("tw_rankfail_i",strofreal(dorank ? sum(i.rankfail) : 0))
}

end

/*--------------------------------------*/
/* Change Log */
/*--------------------------------------*/
*1.004 - mata cv table, CV3J addition, absorb check, return stuff
*1.005 - sample option, make gstar locals conditional
*1.006 - gstar1 caution, large table as matrix
*1.007 - geometric mean
*1.008 - absorb-jack conflict bug corrected
*1.009 - added version number for ssc submission
*1.010 - check for singularities and string cluster variable correction
*2.000 - replace globals with locals, tempvars, absorb-treat conflict resolution
*2.001 - additional global replacements, gstar formatting
*2.002 - fixed constant issue, d.o.f. for cv1
*2.003 - small formatting changes to tables
*3.000 - added regtable and nograph options, fixed bugs
*3.001 - fixed absorb errors
*4.000 - changed handling for singularities
*4.001 - changed absorb calculations to panelsum
*4.002 - debugged absorb and regtable errors
*4.003 - debugged marksample errors
*4.100 - debugged additional absorb error
*4.200 - debugged singular flag error
*4.201 - debugged singular flag absorb error
*4.202 - additional returned values
*4.203 - "singular subsample" terminology used in output
*4.204 - several fixes
*4.208 - cleaned up singular clusters warnings
*4.210 - remove Lg from output when absorb not properly used
*4.22 - regtable bug
*4.22-fast1 - one Mata cluster loop; preallocated beta(-g) and leverage matrices
*             - cholsolve() with invsym() fallback for omit-one systems
*             - vectorized CV3/CV3J, partial leverage, and G* calculations
*             - in-place absorb demeaning, bounded summaries, in-memory graphs
*5.0.0 - optimized Stata-13-compatible implementation
*        - corrected partial leverage and absorb ordering
*        - preallocated fevar dummies and condition-aware solves
*        - rankcheck and full covariance matrices added to r()
*        - corrected CV3J covariance returns and table metadata
*5.1.0 - native two-way G/H/intersection CV1 and CV3 components
*        - scalar max-se inference with min(G,H)-1 degrees of freedom
*        - dimension-specific diagnostics and rank checks
*5.1.1 - audited integration; focal identification guards and private native Mata arrays
*        - corrected native FE labels, double outcome copy, and absorbed CV1 degrees of freedom
*        - explicit maximum-component tie masks and portable certification runner
