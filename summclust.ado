/*------------------------------------*/
/*summclust */
/*written by Matt Webb */
/*version 2.003 2022-05-05 */
/*------------------------------------*/
version 13

cap program drop summclust
program define summclust, rclass

	syntax varlist(numeric) , xvar(varlist)  yvar(varname) cluster(varname)  [fevar(varlist) absorb(varname) gstar  Rho(real 0.1234567) SAMple(string)  TABle SVars JACKknife  ]
	 
	local XTEMP "`varlist'"
	
	local CTEMP "`xvar'"
	
	local FEVAR "`fevar'"
	
	local SMPL "`sample'"
	
	
	/*add if to sample*/
	if "`sample'" != "" {
	    local SMPL = "if `sample'"
	}
	
		
	/*check rho argument*/
	if inrange(`rho',0.00001,1)==0 {
		
		disp as error "The value for rho, `rho', is invalid. rho must be between 0 and 1."
		exit 198
	}
	else {
		/*turns gstar on if rho option is specified*/
		if inrange(`rho',0.1234,0.1235) == 0 {
			local gstarind = 1
		}
		else {
			local gstarind = 0
		}
		
	}
	
	/*flag to calculate g-star*/
	if "`gstar'" != "" {
		local gstarind = 1
	}
	
	/*flag for abosrb used properly*/
	local nocv3 = 0 
	
	/*flag for singular clusters*/
	local singclust = 0
	
	/*declare tempvars*/
	tempvar ytilde
	tempvar temp_sd
	
	local fevars " "
	
	/*run the regression to get the estimation sample*/
	
	local catvars = " "
	foreach catvar in `FEVAR' {
	    local catvars = " `catvars' " +  "i.`catvar'"
	}
	
		/*get est sample and number of vars*/
		qui reg `yvar' `XTEMP' `CTEMP' `catvars' `SMPL', cluster(`cluster')
				 
		tempvar temp_sample 
		qui gen `temp_sample' = e(sample)==1
		
	
	/*for fevar*/
	if `"`fevar'"' != "" {
		
		local tempvars = " "
		local j=0
		foreach fvar in `FEVAR' {
			qui levelsof `fvar' , local(fevarlevel)
			foreach flevel in `fevarlevel' {
			    
				local j = `j'+1
				local tj = "t_`j'"
				
				tempvar `tj'
				qui gen ``tj'' = `fvar'==`flevel'
			}
			
		}
		
		
		local fevars = "`' `t_1'-`t_`j' ' "				
		local XVAR = "`xvar'"
		local CTEMP = "`xvar' `fevars'"
	}
		
		/*sort the data by the clustering variable*/
			/*make a sequential-numeric cluster indicator*/
			
		tempvar temp_indexa	 
		qui egen `temp_indexa' = group(`cluster') if `temp_sample'==1
		qui summ `temp_indexa '
		local G = r(max)
				
		qui sort `temp_indexa'
		qui putmata CVAR=`temp_indexa' if `temp_sample'==1, replace
		
		mata: numobs = rows(CVAR)
		
		mata: xtilde= J(numobs,0,.)
		local xtilde = " "
	
	
	/*for absorb option*/
	if `"`absorb'"' != "" {
		
		local FETEMP "i.`absorb'"
		
		/*partial out all the fixed effects from outcome variable*/
			qui xi: reg `yvar' `FETEMP' if `temp_sample'==1
			qui predict `ytilde' if `temp_sample'==1, res
					
		/*partial out all the fixed effects from variables of interest*/
		foreach xvar in `XTEMP' {
			
			local varname = "xtilde_`xvar'"
			tempvar `varname'
		
			qui reg `xvar' `FETEMP' if `temp_sample'==1
			qui predict ``varname'' if `temp_sample'==1, res
			
			qui putmata `varname' = ``varname''  if `temp_sample'==1, replace
			mata: xtilde = xtilde, `varname'
			
			local xtilde = "`xtilde'" + " ``varname''"
				
		}
		
		/*partial out all the absorbed fixed effects from control variables*/
		local CTEMPER = " " 
		
		foreach cvar in  `XVAR' {
			
			local resname = "res_`cvar'"
			tempvar `resname'
			
			qui reg `cvar' `FETEMP' if `temp_sample'==1
			qui predict ``resname'' if `temp_sample'==1, res
			local CTEMPER = "`CTEMPER' " + " ``resname'' "
						
		}		
		
		/*partial out all the absorbed fixed effects from the other fixed effects*/	
		forvalues i=1/`j' {
			
			local tj = "t_`i'"
			local resname = "res_`i'"
		
			tempvar `resname'
			
			qui reg ``tj'' `FETEMP' if `temp_sample'==1
			qui predict ``resname'' if `temp_sample'==1, res
			qui replace ``tj'' = ``resname''
		}
		local CTEMP = "`CTEMPER' `fevars'"
		
	} /*end of if*/
	else {
		/*if no absorbed variable*/
		qui gen `ytilde' = `yvar' 
		
		foreach xvar in `XTEMP' {
		    local varname = "xtilde_`xvar'"
			tempvar `varname'
			qui gen ``varname''  = `xvar' if `temp_sample'==1
			
			qui putmata `varname' = ``varname'' if `temp_sample'==1, replace
		
			mata: xtilde = xtilde, `varname'
		
			local xtilde = "`xtilde'" + " ``varname''"
		}
				
		/*add constant*/
			tempvar tempcons 
			qui gen `tempcons'  = 1
			local xtilde = "`xtilde'" + " `tempcons' "
		
	} /*end of else*/
	
	*disp "xtilde is `xtilde'"
	

		/*count number of x variables*/
			*qui putmata xtilde = xtilde* if `temp_sample'==1, replace
			mata: xnum = cols(xtilde)
				
		/*check if the fixed effects are constant within cluster*/
		if `"`absorb'"' != "" {
			tempvar temp_indexb 
						
			qui egen `temp_indexb'  = group(`absorb')
		
			qui bysort  `temp_indexb': egen `temp_sd' = sd(`temp_indexa')
			
			qui summ `temp_sd'
			
			local sd_mean = r(mean)
			
			if inrange(`sd_mean',-0.0001,0.0001)==0{
								
				local nocv3 = 1
			}
		}
			
	/*put everything in mata*/
		sort `temp_indexa'
				
		qui putmata X = (`xtilde' `CTEMP' ) if `temp_sample'==1, replace
		
		mata: k = cols(X)
		mata: st_numscalar("k", k)
		mata: st_local("k", strofreal(k))
		mata: numobs = rows(X)
		mata: st_numscalar("numobs", numobs)
		mata: st_local("numobs", strofreal(numobs))	
		
		mata: sk = rank(X)
		mata: st_numscalar("sk", sk)
		mata: st_local("sk", strofreal(sk))	
					
				
		mata: iota =J(rows(X),1,1)
		qui putmata Y = `ytilde' if `temp_sample'==1, replace
		
		mata: info = panelsetup(CVAR,1)
			
		/*create all the submatrices*/
			forvalues g=1/`G' {
			
				mata: X`g' = panelsubmatrix(X,`g',info)
				mata: Y`g' = panelsubmatrix(Y,`g',info)
				mata: iota`g' = panelsubmatrix(iota,`g',info)
				
				mata: X`g'pX`g' = quadcross(X`g',X`g')
				mata: X`g'pY`g' = quadcross(X`g',Y`g')
			}
		
		/*summ over the sub-blocks*/
			mata: XpX = X1pX1 
			mata: XpY = X1pY1
			
			forvalues g=2/`G' {
			    
				mata: XpX = XpX + X`g'pX`g'
				mata: XpY = XpY + X`g'pY`g'
			
			}
			mata: XpXi = invsym(XpX)
			mata: beta = XpXi*XpY
			
		/*generate scores*/
			forvalues g=1/`G' {
			 
				mata: W`g' = Y`g' - X`g' *beta
				*mata: Z`g' = quadcross(X`g',W`g')
				mata: Z`g' = X`g''* W`g'
				mata: Z`g'Z`g'p = Z`g'*Z`g''
			}
			
			mata: ZZp = Z1Z1p
			
			forvalues g=2/`G'{
			
				mata: ZZp = ZZp + Z`g'Z`g'p
			
			}
			
		mata: crve1 = XpXi * ZZp * XpXi
		local gfact = (`G'*(`numobs'-1))/((`G'-1)*(`numobs'-`sk'))
		mata: crve1 =  `gfact':* crve1
		
		
		mata: betatemp = beta[1,1]
		
		mata: st_numscalar("BETATEMP", betatemp)
		mata: st_local("BETATEMP", strofreal(betatemp))	
	 
		mata: cv1se  =  sqrt(crve1[1,1])
		
		mata: st_numscalar("setemp", cv1se)
		mata: st_local("setemp", strofreal(cv1se))	
			 
		local cv1t = `BETATEMP'/`setemp'
		 
		local cv1p = 2*min(ttail(`G'-1, `cv1t'), 1-ttail(`G'-1, `cv1t'))
		 
		local cv1lci = `BETATEMP' - invttail(`G'-1,0.025) * `setemp'
		local cv1uci = `BETATEMP' + invttail(`G'-1,0.025) * `setemp'
		 
		 			
		/*beta-g, Hg, Lg*/
		mata: betags = J(k,0,.)
		mata: Leverage = J(`G',1,.)
		mata: cv3block = J(k,`G',.)
		
		/*string for problem clusters*/
		local probclust = ""
		
		forvalues g=1/`G'{
		
			mata: left = XpX - X`g'pX`g' 
			mata: right = XpY -  X`g'pY`g'
			mata: leftinv = invsym(left)
			mata: beta`g' = leftinv *right
			
			mata: betags = betags, beta`g'
						
			mata: H2`g' = X`g''*X`g'*XpXi
			
			mata: L`g' = trace(H2`g')
			mata: Leverage[`g',1] =  L`g' 
			
			/*cv3*/
			mata: betatemp = (beta :- beta`g')
			mata: cv3block[.,`g'] = betatemp
			
			/*check betag*/
			mata: st_numscalar("betag",beta`g'[1,1] )
			mata: st_local("betag", strofreal(beta`g'[1,1]))
			if inrange(`betag',-0.00000000001,0.00000000001) {
			    local singclust = 1
				local probclust = "`probclust'" + " `g'"
			}  		
		}
				
		/*jackknife*/
		
		if "`jackknife'" != "" {
			
			mata: betagbar = mean(betags')
			mata: betadiff = betags' :- betagbar
			mata: betadiffsum = betadiff[1,.]'* betadiff[1,.]
			forvalues g  = 2 /`G'{
			    mata: betadiffsum = betadiffsum +  betadiff[`g',.]'* betadiff[`g',.]
				
			}
			
			mata: cv3j = sqrt( ((`G'-1)/`G') :* betadiffsum)
			
			mata: cv3jse = cv3j[1,1]
			mata: st_numscalar("cv3jse", cv3jse)
			mata: st_local("cv3jse", strofreal(cv3jse))
			
			local cv3jt = `BETATEMP'/`cv3jse'
			local cv3jp = 2*min(ttail(`G'-1, `cv3jt'), 1-ttail(`G'-1, `cv3jt'))
			
			local cv3jlci = `BETATEMP' - invttail(`G'-1,0.025) * `cv3jse'
			local cv3juci = `BETATEMP' + invttail(`G'-1,0.025) * `cv3jse'
			
		}
	 
		/*partial leverage*/
		qui reg ``varname'' ``tempcons'' `CTEMP' if `temp_sample'==1, noc
				
		local xtt = "xtt"
		tempvar xtt
		qui predict `xtt' if `temp_sample'==1, resid
		
		qui putmata xtt = `xtt' if `temp_sample'==1, replace
		
		mata: xttxtt = cross(xtt,xtt)
		
		mata: pLeverage = J(`G',1,.)
		mata: cv3sum = J(k,k,0)
		
		forvalues g=1/`G'{
			mata: xtt`g' = panelsubmatrix(xtt,`g',info)
			mata: pL`g' = cross(xtt`g',xtt`g') / xttxtt
			mata: pLeverage[`g',1] =  pL`g'
			
			mata: cv3sum = cv3sum + cv3block[.,`g']*cv3block[.,`g']'
		
		}
				
		mata: cv3 = ((`G'-1)/`G'):*cv3sum
		mata: cv3se = sqrt(cv3[1,1])
		mata: st_numscalar("cv3se", cv3se)
		mata: st_local("cv3se", strofreal(cv3se))
		
		local cv3t = `BETATEMP'/`cv3se'
		
		local cv3p = 2*min(ttail(`G'-1, `cv3t'), 1-ttail(`G'-1, `cv3t'))
		
		local cv3lci = `BETATEMP' - invttail(`G'-1,0.025) * `cv3se'
		local cv3uci = `BETATEMP' + invttail(`G'-1,0.025) * `cv3se'
				
		mata: ng = info[.,2] - info[.,1] :+ 1
		mata: betag = betags[1,.]'

		/*matrices to store output*/
			mata: clustsum= J(7,4,.)
			mata: bonus = J(6,4,.)
		
		local SUMVAR "ng Leverage pLeverage betag"
				
		local s = 0
		
		foreach svar in `SUMVAR' {
			
			  local s = `s' + 1
			  
			  local tempsvar = "temp_`svar'"
			  
			  tempvar `tempsvar'
		
			  qui getmata ``tempsvar'' = `svar', force
			  
			  qui summ ``tempsvar'', det
			  			  
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
		    matrix colnames clustsum = Ng Leverage "Partial L." "beta no g" 
			
			/*rowlabels for cluster by cluster matrix*/
				qui levelsof `cluster' `SMPL', local(CLUSTERNAMES)
				
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
							qui bysort  `fvar': egen `temp_sd' = sd(`cluster')
							
							qui summ `temp_sd'
			
							local sd_mean = r(mean)
			
							if inrange(`sd_mean',-0.0001,0.0001)==1{
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
					    
						qui bysort  `absorb': egen `temp_sd' = sd(`cluster')
						
						qui summ `temp_sd'
		
						local sd_mean = r(mean)
		
						if inrange(`sd_mean',-0.0001,0.0001)==1{
							local gstarind = 2
						}
					}
				} /*end of absorb if*/
			} /*end of if 1*/
		
			if `gstarind' == 1 {
		    		
				forvalues g= 1/ `G'{
					
					mata: w`g' = XpXi[.,1]
		
					mata: gamzero`g' =  w`g'' *  X`g'pX`g' * w`g'
					
					mata: gamone`g' = (iota`g'' * X`g' * w`g')' * (iota`g'' * X`g' * w`g')
					
					mata: gamrho`g' = `rho' * gamone`g' + (1-`rho')*gamzero`g'
					
				}
				
				mata: gamzt = 0
				mata: gamrt = 0
				mata: gamot = 0
				
				forvalues g = 1 / `G'{
					
					mata: gamzt = gamzt + gamzero`g'
					mata: gamot = gamot + gamone`g'
					mata: gamrt = gamrt + gamrho`g'
					
				}
				
				mata: gambarz = gamzt/ `G'
				mata: gambaro = gamot / `G'
				mata: gambarr = gamrt / `G'
							
				mata: gammazero = 0
				mata: gammaone = 0
				mata: gammarho = 0
				
				forvalues g = 1 / `G'{
					
					mata: tempone = ((gamone`g' - gambaro)/gambaro)^2
					mata: tempzero = ((gamzero`g' - gambarz)/gambarz)^2
					mata: temprho = ((gamrho`g'- gambarr)/gambarr)^2
					
					mata: gammazero = gammazero + tempzero
					mata: gammaone = gammaone + tempone
					mata: gammarho = gammarho + temprho
					
				}
				
				mata: gammazero = gammazero / `G'
				mata: gammaone = gammaone / `G'
				mata: gammarho = gammarho / `G'
				
				mata: gstarone = `G'/(1 + gammaone)
				mata: gstarzero = `G'/(1 + gammazero)
				mata: gstarrho = `G'/(1 + gammarho)
			
				mata: st_numscalar("gstarzero", gstarzero)
				mata: st_local("gstarzero", strofreal(gstarzero))
						
				mata: st_numscalar("gstarone", gstarone)
				mata: st_local("gstarone", strofreal(gstarone))
				
				mata: st_numscalar("gstarrho", gstarrho)
				mata: st_local("gstarrho", strofreal(gstarrho))
				
		} 
		else if `gstarind' == 2 {
		    
				forvalues g= 1/ `G'{
					
					mata: w`g' = XpXi[.,1]
		
					mata: gamzero`g' =  w`g'' *  X`g'pX`g' * w`g'
					
				}
				
				mata: gamzt = 0
								
				forvalues g = 1 / `G'{
					
					mata: gamzt = gamzt + gamzero`g'
										
				}
				
				mata: gambarz = gamzt/ `G'
											
				mata: gammazero = 0
				
				forvalues g = 1 / `G'{
					
					mata: tempzero = ((gamzero`g' - gambarz)/gambarz)^2
					mata: gammazero = gammazero + tempzero
					
				}
				
				mata: gammazero = gammazero / `G'
				mata: gstarzero = `G'/(1 + gammazero)
			
				mata: st_numscalar("gstarzero", gstarzero)
				mata: st_local("gstarzero", strofreal(gstarzero))	
		}
				  
		/*output display*/
		disp ""
		disp ""
		disp "SUMMCLUST - MacKinnon, Nielsen, and Webb"
		disp " "
		disp "Cluster summary statistics for `XTEMP' when clustered by `cluster'."
		disp "There are `numobs' observations within `G' `cluster' clusters."
				
		/*singular clusters warning*/
		if `singclust' == 1 {
			mata problems = J(0,1,".")
			local clustnames = ""
			foreach problem in `probclust' {
				mata: problems = problems \ cnames[`problem',1]
				local clustnames = "`clustnames' " + "`clustlabel`problem'' "
			}
			
			disp " "
			disp as error "*********************************************************************"
			disp "WARNING -  Beta undefined when certain cluster(s) are omitted"
			disp "*********************************************************************"
			disp "The coefficient on `XTEMP' is undefined when certain cluster(s) are omitted."
			disp "It is replaced by 0.00. The CV3 standard error is not meaningful."
			disp "The problematic `cluster' cluster(s) are: `clustnames'"
			disp " "
		}
		
		
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
			
			/*also replace the beta-no-g and leverage*/
			forvalues i = 1/ 7 {
				matrix clustsum[`i',2] = .
				matrix clustsum[`i',4] = .
			}
			
			
		}
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
		}
		
		matrix colnames cvstuff = "Coeff" "Sd. Err." "t-stat" "P value" CI-lower CI-upper
			
		matlist cvstuff, title(Regression Output) rowtitle(s.e.) ///
			cspec(& %-4s w6 | %9.6f w10 & %9.6f & %6.4f w7 & %6.4f w7 & %9.6f w10 & %9.6f w10 &) ///
			rspec("`RSPECCV'")
				
		if `nocv3'== 1 {
			disp as error "WARNING - CV3 estimates not available as absorb option improperly specified."
		}
		
		mata: cvstuff=st_matrix("cvstuff")
		
		matlist clustsum , title("Cluster Variability") rowtitle(Statistic) ///
			cspec(& %-10s | %8.2f o4 & %9.6f  o4 & %9.6f  w10 & %9.6f  o4 &) ///
			rspec(&-&&&&&-&)
		
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
			if inrange(`rho',0.1234,0.1235)==0 {
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
		
		/*additional summary statistcs*/
		if "`svars'"!=""{
			
			mata: st_matrix("bonus", bonus)
			
				matrix colnames bonus = Ng Leverage "Partial L." "beta no g" 
				matrix rownames bonus = "Harmonic Mean" "Harmonic Ratio" "Geometric Mean" "Geometric Ratio" "Quadratic Mean" "Quadratic Ratio"
				
			matlist bonus, title("Alternative Sample Means and Ratios to Arithmetic Mean") ///
				cspec(& %-10s w15 | %11.3f o4 & %9.6f  w10  & %9.6f  w10 & %9.6f  o4 &) rspec(&|&&&&&|)
			
		} /*end of svars*/		
			
		mata: scall = (ng, Leverage, pLeverage, betag)
		
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
			
		/*stuff to return*/
			return matrix ng = ng
			return matrix lever = Leverage
			return matrix partlev = pLeverage
			return matrix betanog = betag
			
		/*scalars*/
			/*cv1*/
			return local beta =  `BETATEMP'
			return local cv1se = `setemp'
			return local cv1t = `cv1t'
			return local cv1p = `cv1p'
			return local cv1lci = `cv1lci'
			return local cv1uci = `cv1uci' 
			
			/*cv3*/
			return local beta =  `BETATEMP'
			return local cv3se = `cv3se'
			return local cv3t = `cv3t'
			return local cv3p = `cv3p'
			return local cv3lci = `cv3lci'
			return local cv3uci = `cv3uci' 
			
			/*jackknife*/
			if "`jackknife'"!="" {	
			
				return local beta = `BETATEMP'
				return local cv3Jse = `cv3jse'
				return local cv3Jt  = `cv3jt'
				return local cv3Jp = `cv3jp'
				return local cv3Jlci  = `cv3jlci'
				return local cv3Juci = `cv3juci' 
				
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