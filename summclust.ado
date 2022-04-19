/*------------------------------------*/
/*summclust */
/*written by Matt Webb */
/*version 1.0009 2022-04-18 */
/*------------------------------------*/
version 11

cap program drop summclust
program define summclust, rclass

	syntax varlist(numeric) , xvar(varlist)  yvar(varname) cluster(varname)  [fevar(varlist) absorb(varname) gstar  Rho(real 0.1234567) SAMple(string)  TABle SVars JACKknife  ]
	 
	global XTEMP "`varlist'"
	
	global CTEMP "`xvar'"
	
	global FEVAR "`fevar'"
	
	global SMPL "`sample'"
	
	
	/*add if to sample*/
	if "`sample'" != "" {
	    global SMPL = "if `sample'"
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
	
	/*for fevar*/
	if `"`fevar'"' != "" {
		
		foreach fvar in $FEVAR {
			
			/*create temporary FE vars*/
			qui tab `fvar', gen (t__fe_`fvar'_)
		
		}
		
		qui ds t__fe_*
		local fevars = r(varlist)
	 
		global CTEMP "`xvar' `fevars'"	
	}
	
		/*run the regression to get the estimation sample*/
	
		qui reg `yvar' $XTEMP $CTEMP $SMPL, cluster(`cluster')
		 
		qui gen temp_sample = e(sample)==1 
		
	
	/*for absorb option*/
	if `"`absorb'"' != "" {
		
		global FETEMP "i.`absorb'"
		
		/*partial out all the fixed effects from outcome variable*/
			qui xi: reg `yvar' $FETEMP if temp_sample==1
			qui predict ytilde if temp_sample==1, res
		
		/*partial out all the fixed effects from variables of interest*/
		foreach xvar in $XTEMP {
			qui reg `xvar' $FETEMP if temp_sample==1
			qui predict xtilde_`xvar' if temp_sample==1, res
		}
					
		/*partial out all the fixed effects from control variables*/
		foreach cvar in $CTEMP  {
			
			*disp "cvar is `cvar'"
			qui reg `cvar' $FETEMP if temp_sample==1
			qui predict t99_`cvar' if temp_sample==1, res
			
		}
		
		/*check if the fixed effects are constant within cluster*/
			qui egen temp_indexb  = group(`absorb')
		
			*qui bysort `cluster': egen temp_sd = sd(temp_indexb)
			qui bysort  temp_indexb: egen temp_sd = sd(`cluster')
			
			
			qui summ temp_sd
			
			local sd_mean = r(mean)
			
			if inrange(`sd_mean',-0.0001,0.0001)==0{
								
				local nocv3 = 1
			}
		
		
	} /*end of if*/
	else {
		/*if no absorbed variable*/
		qui gen ytilde = `yvar'
		
		foreach xvar in $XTEMP {
			qui gen xtilde_`xvar'  = `xvar' if temp_sample==1
		}
		
		foreach cvar in $CTEMP  {
			qui gen  t99_`cvar' = `cvar' if temp_sample==1
		}
		
		/*add constant*/
			qui gen t99_cons  = 1
		
	} /*end of else*/

		/*sort the data by the clustering variable*/
			/*make a sequential-numeric cluster indicator*/
			 
			qui egen temp_indexa = group(`cluster') if temp_sample==1
			qui summ temp_indexa 
			global G = r(max)
					
			qui sort temp_indexa
			qui putmata CVAR=temp_indexa if temp_sample==1, replace
			
		/*count number of x variables*/
			qui putmata xtilde = xtilde* if temp_sample==1, replace
			mata: xnum = cols(xtilde)
			
	/*put everything in mata*/
		qui putmata X = (xtilde* t99* ) if temp_sample==1, replace
		mata: k = cols(X)
		mata: st_numscalar("k", k)
		mata: st_local("k", strofreal(k))
		mata: numobs = rows(X)
		mata: st_numscalar("numobs", numobs)
		mata: st_local("numobs", strofreal(numobs))	
						
		/*check regression after partialing out*/
		 qui reg ytilde  xtilde* t99* , noc  cluster(`cluster')
		 		 
		 global BETATEMP = _b[xtilde]
		
		 local setemp = _se[xtilde]
		 
		 local cv1t = $BETATEMP/`setemp'
		 
		 local cv1p = 2*min(ttail($G-1, `cv1t'), 1-ttail($G-1, `cv1t'))
		 
		local cv1lci = $BETATEMP - invttail(${G}-1,0.025) * `setemp'
		local cv1uci = $BETATEMP + invttail(${G}-1,0.025) * `setemp'
		 
		 mata: beta = st_matrix("e(b)")
		 mata: beta = beta'	 
		 
		/*get residuals for scores*/
		qui predict uhat if temp_sample==1, residual 
		qui putmata U=uhat if temp_sample==1, replace
		
				
		mata: iota =J(rows(X),1,1)
		qui putmata Y = ytilde if temp_sample==1, replace
		
		mata: info = panelsetup(CVAR,1)
			
		/*create all the submatrices*/
			forvalues g=1/$G {
			
				mata: X`g' = panelsubmatrix(X,`g',info)
				mata: Y`g' = panelsubmatrix(Y,`g',info)
				mata: U`g' = panelsubmatrix(U,`g',info)
				mata: iota`g' = panelsubmatrix(iota,`g',info)
				
				mata: X`g'pX`g' = quadcross(X`g',X`g')
				mata: X`g'pY`g' = quadcross(X`g',Y`g')
				mata: S`g' = quadcross(X`g',U`g')
				mata: S`g'S`g'p = S`g'*S`g''
			}
		
		/*summ over the sub-blocks*/
			mata: XpX = X1pX1 
			mata: XpY = X1pY1
			mata: SSp = S1S1p
			
			forvalues g=2/$G {
			
				mata: XpX = XpX + X`g'pX`g'
				mata: XpY = XpY + X`g'pY`g'
				mata: SSp = SSp + S`g'S`g'p
			
			}
			mata: XpXi = invsym(XpX)
				
		/*beta-g, Hg, Lg*/
		mata: betags = J(k,0,.)
		mata: Leverage = J($G,1,.)
		mata: cv3block = J(k,$G,.)
		
		forvalues g=1/$G {
		
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
					
		}
		
		/*jackknife*/
		
		if "`jackknife'" != "" {
			
			mata: betagbar = mean(betags')
			mata: betadiff = betags' :- betagbar
			mata: betadiffsum = betadiff[1,.]'* betadiff[1,.]
			forvalues g  = 2 /$G {
			    mata: betadiffsum = betadiffsum +  betadiff[`g',.]'* betadiff[`g',.]
				
			}
			
			mata: cv3j = sqrt( (($G-1)/$G) :* betadiffsum)
			
			mata: cv3jse = cv3j[1,1]
			mata: st_numscalar("cv3jse", cv3jse)
			mata: st_local("cv3jse", strofreal(cv3jse))
			
			local cv3jt = $BETATEMP/`cv3jse'
			local cv3jp = 2*min(ttail($G-1, `cv3jt'), 1-ttail($G-1, `cv3jt'))
			
			local cv3jlci = $BETATEMP - invttail(${G}-1,0.025) * `cv3jse'
			local cv3juci = $BETATEMP + invttail(${G}-1,0.025) * `cv3jse'
			
		}
	 

		
		/*partial leverage*/
		qui reg xtilde t99* if temp_sample==1, noc
		qui predict xtt if temp_sample==1, resid
		qui putmata xtt if temp_sample==1, replace
		mata: xttxtt = cross(xtt,xtt)
		
		mata: pLeverage = J($G,1,.)
		mata: cv3sum = J(k,k,0)
		
		
		forvalues g=1/$G {
			mata: xtt`g' = panelsubmatrix(xtt,`g',info)
			mata: pL`g' = cross(xtt`g',xtt`g') / xttxtt
			mata: pLeverage[`g',1] =  pL`g'
			
			mata: cv3sum = cv3sum + cv3block[.,`g']*cv3block[.,`g']'
		
		}
				
		mata: cv3 = (($G-1)/$G):*cv3sum
		mata: cv3se = sqrt(cv3[1,1])
		mata: st_numscalar("cv3se", cv3se)
		mata: st_local("cv3se", strofreal(cv3se))
		
		local cv3t = $BETATEMP/`cv3se'
		
		local cv3p = 2*min(ttail($G-1, `cv3t'), 1-ttail($G-1, `cv3t'))
		
		local cv3lci = $BETATEMP - invttail(${G}-1,0.025) * `cv3se'
		local cv3uci = $BETATEMP + invttail(${G}-1,0.025) * `cv3se'
				
		mata: ng = info[.,2] - info[.,1] :+ 1
		mata: betag = betags[1,.]'

		/*matrices to store output*/
			mata: clustsum= J(7,4,.)
			mata: bonus = J(6,4,.)
		
		global SUMVAR "ng Leverage pLeverage betag"
				
		local s = 0
		
		foreach svar in $SUMVAR {
			
			  local s = `s' + 1
		
			  qui getmata temp_`svar' = `svar', force
			  
			  qui summ temp_`svar', det
			  
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
				  mata: denom= 1/(($G-1)*(`mean'^2))
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
			
		  qui drop temp_`svar'
					
		} /*end of sum var*/
			
			mata: st_matrix("clustsum", clustsum)
			matrix rownames clustsum = min q1 median mean q3 max coefvar
		    matrix colnames clustsum = Ng Leverage "Partial L." "beta no g" 
			
			/*rowlabels for cluster by cluster matrix*/
				qui levelsof `cluster' $SMPL, local(CLUSTERNAMES)
				
				mata: cnames = J($G,1,".")
				local f = 0
				foreach cnames in `CLUSTERNAMES'{
				    local f = `f' + 1
					mata: cnames[`f',1] = "`cnames'"
				}
				
		
		
		
		/****************************/
		/*gstar option*/
		/***************************/
			/*check if fixed effects invariant within clusters*/
			if `gstarind' == 1 {
			 
				if `"`fevar'"' != "" {
					
					foreach fvar in $FEVAR {
						
					
						if "`fvar'" == "`cluster'" {
							local gstarind = 2
						}
						else {
						    qui cap drop temp_sd
							qui bysort  `fvar': egen temp_sd = sd(`cluster')
							
							qui summ temp_sd
			
							local sd_mean = r(mean)
			
							if inrange(`sd_mean',-0.0001,0.0001)==1{
								local gstarind = 2
							}
							
						}
						
						
					}
				} /*end of fevar if*/
				
				if `"`absorb'"' != "" {
						
					if "`absorb'" == "`cluster'" {
						local gstarind = 2
					}
					else {
						qui cap drop temp_sd
						qui bysort  `absorb': egen temp_sd = sd(`cluster')
						
						qui summ temp_sd
		
						local sd_mean = r(mean)
		
						if inrange(`sd_mean',-0.0001,0.0001)==1{
							local gstarind = 2
						}
						
					}
						
				} /*end of absorb if*/
			

			} /*end of if 1*/
		
			if `gstarind' == 1 {
		    		
				forvalues g= 1/ $G {
					
					mata: w`g' = XpXi[.,1]
		
					mata: gamzero`g' =  w`g'' *  X`g'pX`g' * w`g'
					
					mata: gamone`g' = (iota`g'' * X`g' * w`g')' * (iota`g'' * X`g' * w`g')
					
					mata: gamrho`g' = `rho' * gamone`g' + (1-`rho')*gamzero`g'
					
				}
				
				mata: gamzt = 0
				mata: gamrt = 0
				mata: gamot = 0
				
				forvalues g = 1 / $G {
					
					mata: gamzt = gamzt + gamzero`g'
					mata: gamot = gamot + gamone`g'
					mata: gamrt = gamrt + gamrho`g'
					
				}
				
				mata: gambarz = gamzt/ $G
				mata: gambaro = gamot / $G
				mata: gambarr = gamrt / $G
							
				mata: gammazero = 0
				mata: gammaone = 0
				mata: gammarho = 0
				
				forvalues g = 1 / $G {
					
					mata: tempone = ((gamone`g' - gambaro)/gambaro)^2
					mata: tempzero = ((gamzero`g' - gambarz)/gambarz)^2
					mata: temprho = ((gamrho`g'- gambarr)/gambarr)^2
					
					mata: gammazero = gammazero + tempzero
					mata: gammaone = gammaone + tempone
					mata: gammarho = gammarho + temprho
					
				}
				
				mata: gammazero = gammazero / $G
				mata: gammaone = gammaone / $G
				mata: gammarho = gammarho / $G
				
				mata: gstarone = $G/(1 + gammaone)
				mata: gstarzero = $G/(1 + gammazero)
				mata: gstarrho = $G/(1 + gammarho)
			
				
				mata: st_numscalar("gstarzero", gstarzero)
				mata: st_local("gstarzero", strofreal(gstarzero))
						
				mata: st_numscalar("gstarone", gstarone)
				mata: st_local("gstarone", strofreal(gstarone))
				
				mata: st_numscalar("gstarrho", gstarrho)
				mata: st_local("gstarrho", strofreal(gstarrho))
				
		
		} 
		else if `gstarind' == 2 {
		    
				forvalues g= 1/ $G {
					
					mata: w`g' = XpXi[.,1]
		
					mata: gamzero`g' =  w`g'' *  X`g'pX`g' * w`g'
					
				}
				
				mata: gamzt = 0
								
				forvalues g = 1 / $G {
					
					mata: gamzt = gamzt + gamzero`g'
										
				}
				
				mata: gambarz = gamzt/ $G
											
				mata: gammazero = 0
				
				
				forvalues g = 1 / $G {
					
					mata: tempzero = ((gamzero`g' - gambarz)/gambarz)^2
					mata: gammazero = gammazero + tempzero
					
				}
				
				mata: gammazero = gammazero / $G
				mata: gstarzero = $G/(1 + gammazero)
			
				mata: st_numscalar("gstarzero", gstarzero)
				mata: st_local("gstarzero", strofreal(gstarzero))
				
		}
				  
		/*output display*/
		disp ""
		disp ""
		disp "SUMMCLUST - MacKinnon, Nielsen, and Webb"
		disp " "
		disp "Cluster summary statistics for $XTEMP when clustered by `cluster'."
		disp "There are `numobs' observations within ${G} `cluster' clusters."
		
		
		if "`jackknife'"!="" {	
			matrix define cvstuff = J(3,6,.)
			matrix rownames cvstuff = CV1 CV3 CV3J
			global RSPECCV "&|&&|"
		}
		else {
			matrix define cvstuff = J(2,6,.)
			matrix rownames cvstuff = CV1 CV3
			global RSPECCV "&|&|"
		}

		/*check if absorb used improperly*/
		if `nocv3'==1 {
			
			matrix cvstuff[1,1] = $BETATEMP
			matrix cvstuff[1,2] = `setemp'
			matrix cvstuff[1,3] = `cv1t'
			matrix cvstuff[1,4] = `cv1p'
			matrix cvstuff[1,5] = `cv1lci'
			matrix cvstuff[1,6] = `cv1uci' 
			
					
			matrix cvstuff[2,1] = $BETATEMP
			
		
			
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
			
			matrix cvstuff[1,1] = $BETATEMP
			matrix cvstuff[1,2] = `setemp'
			matrix cvstuff[1,3] = `cv1t'
			matrix cvstuff[1,4] = `cv1p'
			matrix cvstuff[1,5] = `cv1lci'
			matrix cvstuff[1,6] = `cv1uci' 
			
					
			matrix cvstuff[2,1] = $BETATEMP
			matrix cvstuff[2,2] = `cv3se'
			matrix cvstuff[2,3] = `cv3t'
			matrix cvstuff[2,4] = `cv3p'
			matrix cvstuff[2,5] = `cv3lci'
			matrix cvstuff[2,6] = `cv3uci' 
			
			if "`jackknife'"!="" {	
				
				matrix cvstuff[3,1] = $BETATEMP
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
			rspec("$RSPECCV")
				
		if `nocv3'== 1 {
			disp as error "WARNING - CV3 estimates not available as absorb option misproperly specified."
		}
		
		mata: cvstuff=st_matrix("cvstuff")
		
		matlist clustsum , title("Cluster Variability") rowtitle(Statistic) ///
			cspec(& %-10s | %8.1f o4 & %9.6f  o4 & %9.6f  w10 & %9.6f  o4 &) ///
			rspec(&-&&&&&-&)
		
		/*error for no-beta-g*/
			if `nocv3'==1 {
				disp as error "Some results not available as absorb option is misproperly specified."
				disp as text ""
			}
		
		/*drop created variables*/
		cap drop temp_indexa ytilde xtilde* t99* uhat
		cap drop xtt
		cap drop t__fe_*
		cap drop temp_sd
		cap drop temp_indexb
		cap drop temp_sample
		
		/*-------------------*/
		/*optional output*/
		/*-------------------*/
		
		/*gstar*/
		if `gstarind'==1 {	
		    
			disp " "
			disp "Effective Number of Clusters"
			disp "-----------------------------"
			disp  "G*(0)  is " " `gstarzero'"
			if inrange(`rho',0.1234,0.1235)==0 {
				disp   "G*(`rho') is `gstarrho' "
			}
			 disp  "G*(1)  is `gstarone' "
			
		}
		else if `gstarind'==2 {
			disp " "
			disp "Effective Number of Clusters"
			disp "-----------------------------"
			disp  "G*(0)  is " " `gstarzero'"
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
			
			if $G < 53 {
			    
				mata: st_matrix("scall", scall)
				matrix rownames scall = `CLUSTERNAMES'
				matrix colnames scall =  Ng Leverage "Partial L." "beta no g"  
			
				global TITLE "Cluster by Cluster Statistics"
				global Gm1 = $G - 1 
			
				
				
				local rspec : display "rspec(&-"_dup(${Gm1}) "&" "-)"
				
				
				matlist scall, title("${TITLE}") rowtitle("`cluster'") ///
					cspec(& %-10s | %8.1f o4 & %9.6f  o4 & %9.6f  w10 & %9.6f  o4 &) `rspec'
				
			
				return matrix scall = scall
			}
			else {
			    
			
				disp " "
				disp "Cluster by Cluster Statistics"
				
				disp as error "The number of clusters, $G exceeds 52, results reported as matrix. "
				
				disp as text ""
				
				*mata: bigtab = cnames, "scall"
				
				*mata: bigtab
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
				return local beta =  $BETATEMP
				return local cv1se = `setemp'
				return local cv1t = `cv1t'
				return local cv1p = `cv1p'
				return local cv1lci = `cv1lci'
				return local cv1uci = `cv1uci' 
				
				/*cv3*/
				return local beta =  $BETATEMP
				return local cv3se = `cv3se'
				return local cv3t = `cv3t'
				return local cv3p = `cv3p'
				return local cv3lci = `cv3lci'
				return local cv3uci = `cv3uci' 
				
				/*jackknife*/
				if "`jackknife'"!="" {	
				
					return local beta = $BETATEMP
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
*1.0004 - mata cv table, CV3J addition, absorb check, return stuff
*1.0005 - sample option, make gstar locals conditional
*1.0006 - gstar1 caution, large table as matrix
*1.0007 - geometric mean
*1.0008 - absorb-jack conflict bug corrected
*1.0009 - added version number for ssc submission