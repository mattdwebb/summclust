/*------------------------------------*/
/*summclust */
/*written by Matt Webb */
/*version 4.210 2023-06-06 */
/*------------------------------------*/
version 13

cap program drop summclust
program define summclust, rclass

	syntax varlist(min=3), cluster(varname) [fevar(varlist) absorb(varname) gstar  Rho(real 0.1234567) SAMple(string)  TABle ADDmeans JACKknife REGtable NOGraph]
	
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
		
	/*for fevar*/
	if `"`fevar'"' != "" {
		
		*regtable names labels
		local names_fe = ""
		
		local tempvars = " "
		local j=0
		
		local q = 0 
		
		foreach fvar in `FEVAR' {
			
			local q = `q' + 1
			
			local countfevar = `countfevar' + 1
			qui levelsof `fvar' if `temp_sample'==1, local(fevarlevel)
			
			local jstart = `j'+1
			
			foreach flevel in `fevarlevel' {
			    
				local j = `j'+1
				local tj = "t_`j'"
				
				tempvar `tj'
				qui gen ``tj'' = `fvar'==`flevel'
				
				/*create labels for these variables for regtable*/
				local name_tj = "`fvar'_`flevel'"
								
				local names_fe = "`names_fe' `name_tj'"	
			
			} /*end of flevel*/
					
			/*drop the last fe for the second and greater fe*/
			local jend = `j'-1
						
			if `q'==1 & `"`absorb'"' == "" {
				local fevars = "`' `t_1'-`t_`j' ' "
				
				local names_fe = subinstr("`names_fe'", "`name_tj'"," ", .)
				
			}
			else if `q'==1 & `"`absorb'"' != "" {
				local fevars = "`' `t_1'-`t_`jend' ' "
				
				local names_fe = subinstr("`names_fe'", "`name_tj'"," ", .)
				
			}
			
			else if `q' > 1 {
				local fevars = "`fevars' `t_`jstart'' - `t_`jend''  "
						 
				local names_fe = subinstr("`names_fe'", "`name_tj'"," ", .)	
			}
			
		} /*end of fvar*/ 
			
		*local fevars = "`' `t_1'-`t_`j' ' "				
		local XVAR = "`CTEMP'"
		local CTEMP = "`CTEMP' `fevars'"
		
	} /*end of fevar*/
		
	/*sort the data by the clustering variable*/
		/*make a sequential-numeric cluster indicator*/
				
	tempvar temp_indexa	 
	qui egen `temp_indexa' = group(`cluster') if `temp_sample'==1
	qui summ `temp_indexa'
	local G = r(max)
			
	qui sort `temp_indexa'
	qui putmata CVAR=`temp_indexa' if `temp_sample'==1, replace
	
	mata: numobs = rows(CVAR)
	
	mata: info = panelsetup(CVAR,1)
	mata: ng = info[.,2] - info[.,1] :+ 1
	
	mata: xtilde= J(numobs,0,.)
	local xtilde = " "
	
	local XVAR = "`CTEMP'"
	
	/*dump everything into mata here, then demean by absorb*/
	
	qui gen `ytilde' = `yvar' if `temp_sample'==1
	
	foreach xvar in `XTEMP' {
				
		qui putmata `xvar' = `xvar' if `temp_sample'==1, replace
	
		mata: xtilde = xtilde, `xvar'
	
		local xtilde = "`xtilde'" + " ``varname''"
		
	}
	
	
	if `"`absorb'"' != "" {
		
		/* want to check that cluster is constant within absorb */
		
		qui putmata absorb =`absorb' if `temp_sample'==1, replace
		
		qui putmata cabsorb = `temp_indexa' if `temp_sample'==1, replace
		
		mata: ainfo = panelsetup(absorb,1)
		
		mata: ang = ainfo[.,2] - ainfo[.,1] :+ 1
		
		mata: absave = panelsum(cabsorb,ainfo) :/ ang
		
		mata: absorb2 = cabsorb:*cabsorb
		
		mata: absorb2ave = panelsum(absorb2, ainfo) :/ ang

		mata absconstant = floatround(absorb2ave)!=floatround( absave:*absave)
		
		mata: st_numscalar("nocv3", absconstant)
		mata: st_local("nocv3", strofreal(absconstant))	
		
		
		/*other way*/
							
		if "`absorb'" == "`cluster'" {
			local nogood = 0
		}
		else {
			
			/*check for invariance of cluster variable*/
			cap drop  `temp_sd'
			qui bysort  `absorb': egen `temp_sd' = sd(`cluster')
			
			qui summ `temp_sd'

			local sd_mean = r(mean)

			if inrange(`sd_mean',-0.0001,0.0001)==1{
				local nogood = 0
			}
			else{
				local nogood = 1
				disp as error "Cluster variable not constant within absorb variable.  Use fevar instead."
				exit 198
				
			}
		}				
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
		sort `temp_indexa'
		
		qui putmata Y = `ytilde' if `temp_sample'==1, replace
		
		qui putmata XOTHER = (`CTEMP' ) if `temp_sample'==1, replace
		
		mata: X = xtilde, XOTHER
				
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
		
		/*demean y */

		mata: ymean = panelsum(Y,abinfo) :/ abng
		
		mata: Ydm = J(0,1,0)
				
		forvalues g=1/`H' {
			
			mata: Ydm`g' = panelsubmatrix(Y,`g',abinfo)
			
			mata: Ydm`g' = Ydm`g' :- ymean[`g',1]
		
			mata: Ydm = Ydm \ Ydm`g'	
		}
		
		mata Y = Ydm
		
		/*demean x*/
		
		mata: xmean = panelsum(X, abinfo) :/ abng
		
		mata: Xdm = J(0,cols(X),0)
				
		forvalues g=1/`H' {
			
			mata: Xdm`g' = panelsubmatrix(X,`g',abinfo)
			
			mata: Xdm`g' = Xdm`g' :- xmean[`g',.]
		
			mata: Xdm = Xdm \ Xdm`g'	
		}
		
		mata X = Xdm	
	
	} /*end of  absorb if*/
	
		/*count number of x variables*/
			mata: xnum = cols(xtilde)
						
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
				mata: Z`g' = X`g''* W`g'
				mata: Z`g'Z`g'p = Z`g'*Z`g''
			}
			
			mata: ZZp = Z1Z1p
			
			forvalues g=2/`G'{
				mata: ZZp = ZZp + Z`g'Z`g'p
			}
		
		
		mata: crve1 = XpXi * ZZp * XpXi
		local gfact = (`G'*(`numobs'-1))/((`G'-1)*(`numobs'-`coeffnum'))
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
		mata: betadropgs = J(k,0,.)

		mata: Leverage = J(`G',1,.)
		mata: cv3block = J(k,`G',.)
		
		/*string for problem clusters*/
		local probclust = ""
		local probclustfull = ""
		
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
			
			/*check for singularities for all elements*/
			mata anyzerosmall =  beta`g':<=  0.00000001
			mata anyzerobig   =  beta`g':>= -0.00000001
			mata anyzero =  anyzerosmall :* anyzerobig 
	
			mata countzero = colsum(anyzero)
				
			mata: st_numscalar("countzero", countzero)
			mata: st_local("countzero", strofreal(countzero))
				
			if `countzero' > 0 {
				
				mata: betadrop`g' = J(rows(beta`g'),1,0)
				
				local singclustfull = `singclustfull'  + 1

				local probclustfull = "`probclustfull'" + " `g'"
			}
			
			else {
				mata: betadrop`g' = beta`g'
			}
			
			mata: betadropgs = betadropgs, betadrop`g'
					
		} /*end of g loop*/
					
		local singfrac = (`singclustfull'/ `G')*100
		if `singclustfull' >=1 {
			
			mata: betadropgs = select(betadropgs, betadropgs[1.,]:!=0)
				
			mata gdrop = cols(betadropgs)
			mata: st_numscalar("gdrop", gdrop)
			mata: st_local("gdrop", strofreal(gdrop))
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
			
			/*dropped singular version*/
			if `singclustfull' >= 1 {
			
				mata: betagdropbar = mean(betadropgs')
				mata: betadiffdrop = betadropgs' :- betagdropbar
				mata: betadiffsumdrop = betadiffdrop[1,.]'* betadiffdrop[1,.]
				forvalues g  = 2 /`gdrop'{
					mata: betadiffsumdrop = betadiffsumdrop +  betadiffdrop[`g',.]'* betadiffdrop[`g',.]
				}
				
				mata: cv3jdrop = sqrt( ((`gdrop'-1)/`gdrop') :* betadiffsumdrop)
				
				mata: cv3jsedrop = cv3jdrop[1,1]
				mata: st_numscalar("cv3jsedrop", cv3jsedrop)
				mata: st_local("cv3jsedrop", strofreal(cv3jsedrop))
				
				local cv3jtdrop = `BETATEMP'/`cv3jsedrop'
				local cv3jpdrop = 2*min(ttail(`gdrop'-1, `cv3jtdrop'), 1-ttail(`gdrop'-1, `cv3jtdrop'))
				
				local cv3jlcidrop = `BETATEMP' - invttail(`gdrop'-1,0.025) * `cv3jsedrop'
				local cv3jucidrop = `BETATEMP' + invttail(`gdrop'-1,0.025) * `cv3jsedrop'
				
			} /*end of singclustfull*/
			
		} /*end of jackknife*/
			 
		/*partial leverage*/
		 
		qui reg ``varname'' ``tempcons'' `CTEMP' if `temp_sample'==1, noc
				
		local xtt = "xtt"
		tempvar xtt
		qui predict `xtt' if `temp_sample'==1, resid
		
		qui putmata xtt = `xtt' if `temp_sample'==1, replace
		
		mata: xttxtt = cross(xtt,xtt)
		
		mata: pLeverage = J(`G',1,.)
		mata: cv3sum = J(k,k,0)
		mata: cv3sumdrop = J(k,k,0)
		
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
		
		local SUMVAR "ng Leverage pLeverage betag"
		
		/*matrices to store output*/
			mata: clustsum= J(7,4,.)
			mata: bonus = J(6,4,.)
		
		/*dropped singular version*/
		if `singclustfull' >=1 {
			
			forvalues g = 1/`gdrop' {
				mata: betatempdrop = (beta :- betadropgs[.,`g'])
				mata: cv3sumdrop = cv3sumdrop +   betatempdrop*betatempdrop'
			}
			
			mata: cv3drop = ((`gdrop'-1)/`gdrop'):*cv3sumdrop
			mata: cv3sedrop = sqrt(cv3drop[1,1])
			mata: st_numscalar("cv3sedrop", cv3sedrop)
			mata: st_local("cv3sedrop", strofreal(cv3sedrop))
			
			local cv3tdrop = `BETATEMP'/`cv3sedrop'
			
			local cv3pdrop = 2*min(ttail(`gdrop'-1, `cv3tdrop'), 1-ttail(`gdrop'-1, `cv3tdrop'))
			
			local cv3lcidrop = `BETATEMP' - invttail(`gdrop'-1,0.025) * `cv3sedrop'
			local cv3ucidrop = `BETATEMP' + invttail(`gdrop'-1,0.025) * `cv3sedrop'	
			
			/*statistics adding dropped betas*/
			mata: colsdrop = cols(betadropgs)
			mata: st_numscalar("colsdrop", colsdrop)
			mata: st_local("colsdrop", strofreal(colsdrop))
			
			if `colsdrop' > 0 {
				
				local SUMVAR "ng Leverage pLeverage betag betadropg"
			
				mata: clustsum= J(7,5,.)
				mata: bonus = J(6,5,.)
			
				mata: betadropg = betadropgs[1,.]'
				
			}
			
		
					
		} /*end of singclustfull */
					
		mata: betag = betags[1,.]'
		
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
					    
						/*check for invariance of cluster variable*/
					    cap drop  `temp_sd'
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
						
			/*the cv3J matrix is the square root already, CV3 is not*/
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
				
				/*include constant row when absorb & fevar options not used */
				if `"`absorb'"' == "" & `"`fevar'"' == "" {
					
					matrix rownames regresstab = "`XTEMP' " "_cons" `name_xvar' `names_fe'
				}
		
				else {
					matrix rownames regresstab = "`XTEMP' " `name_xvar' `names_fe'
				}
				
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
				
				mata fullcv3 = diagonal(cv3j)
												
				mata: cv3jsq = cv3j:*cv3j
				mata: st_matrix("cv3j",cv3jsq)
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
				
				/*include constant row when absorb & fevar options not used */
				if `"`absorb'"' == "" & `"`fevar'"' == "" {
					
					matrix rownames regresstab = "`XTEMP'" "_cons" `name_xvar' `names_fe'
				}
				else {
					matrix rownames regresstab = "`XTEMP' "  `name_xvar' `names_fe'
				}

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
				
			} /*end of yes jackknife if */
			
			/*singular clusters*/
			if `singclustfull' >= 1 {
							
				/*the cv3J matrix is the square root already, CV3 is not*/
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
					
					/*include constant row when absorb & fevar options not used */
					if `"`absorb'"' == "" & `"`fevar'"' == "" {
						
						matrix rownames regresstabdrop = "`XTEMP' " "_cons"  `name_xvar' `names_fe'
					}
					else {
						matrix rownames regresstabdrop = "`XTEMP' " `name_xvar' `names_fe'
					}
					
					if `rtrows' <= 50 {
						local rspecreg: display "-" _dup(`rtrows') "&" "-"
											
						matlist
						regresstabdrop, title(Linear Regression -- CV3 -- Omitting Singular Omit-One-Cluster Subsamples) rowtitle(`yvar') ///
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
					
					mata fullcv3drop = diagonal(cv3jdrop)
													
					mata: cv3jsqdrop = cv3jdrop:*cv3jdrop
					mata: st_matrix("cv3jdrop",cv3jsqdrop)
					return matrix cv3jdrop = cv3jdrop
					
					mata: regresstabdrop = J(rows(fullcv3),6,.)
					mata: regresstabdrop[.,1] =beta
					mata: regresstabdrop[.,2] = fullcv3drop
					mata: tstatsdrop = beta :/ fullcv3drop
					mata: regresstabdrop[.,3] = tstatsdrop
					
					mata: pvalsdrop = J(rows(fullcv3drop),2,.)
					
					forvalues i = 1 / `sk' {
						
						mata: pvalsdrop[`i',1] =  ttail(`G'-1,regresstabdrop[`i',3])
						mata: pvalsdrop[`i',2] = 1 -  ttail(`G'-1,regresstabdrop[`i',3])
						mata: regresstabdrop[`i',4] = 2*min(pvals[`i',.])
					}
					
					mata: regresstabdrop[.,5] = regresstabdrop[.,1] - invttail(`G'-1,.025)* regresstabdrop[.,2]
					mata: regresstabdrop[.,6] = regresstabdrop[.,1] + invttail(`G'-1,.025)* regresstabdrop[.,2]
									
					mata: st_matrix("regresstabdrop",regresstabdrop)
					matrix colnames regresstabdrop = "Coeff" "Sd. Err." "t-stat" "P value" CI-lower CI-upper
					
					/*dimensions of this matrix*/
					mata: rtrows = rows(regresstabdrop)
					mata: st_numscalar("rtrows", rtrows)
					mata: st_local("rtrows", strofreal(rtrows))	
					
					/*include constant row when absorb & fevar options not used */
					if `"`absorb'"' == "" & `"`fevar'"' == "" {
						
						matrix rownames regresstabdrop = "`XTEMP' " "_cons" `name_xvar' `names_fe'
					}
					else {
						matrix rownames regresstabdrop = "`XTEMP' "  `name_xvar' `names_fe'
					}
									
					if `rtrows' <= 50 {
						local rspecreg: display "-" _dup(`rtrows') "&" "-"

						matlist
						regresstabdrop, title(Linear Regression -- CV3 -- Omitting Singular Omit-One-Cluster Subsamples) rowtitle(`yvar') ///
						cspec(& %-12s w6 | %9.6f w10 & %9.6f & %6.4f w7 & %6.4f w7 & %9.6f w10 & %9.6f w10 &) ///
						rspec(`rspecreg')
					}
					else{
						disp "Linear Regression -- CV3 -- Omitting Singular Omit-One-Cluster Subsamples"
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
		
	
		/*stuff to return*/
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

			qui cap graph drop name summclust_temp_ln.gph
			qui cap graph drop name summclust_temp_pn.gph  
			qui cap graph drop name summclust_temp_lb.gph
			qui cap graph drop name summclust_temp_pb.gph 

			qui cap erase summclust_temp_ln.gph
			qui cap erase summclust_temp_pn.gph
			qui cap erase summclust_temp_lb.gph
			qui cap erase summclust_temp_pb.gph

			set graphics off
			qui twoway scatter `lever' `ng', saving(summclust_temp_ln) scheme(s1mono) xline( `graph_1' ) yline( `graph_2' )
			qui twoway scatter `plever' `ng', saving(summclust_temp_pn) scheme(s1mono) xline( `graph_1' ) yline(`graph_3' )
			qui twoway scatter `lever' `betanog', saving(summclust_temp_lb) scheme(s1mono) xline( `graph_4' ) yline( `graph_2' )
			qui twoway scatter `plever' `betanog', saving(summclust_temp_pb) scheme(s1mono) xline(`graph_4' ) yline( `graph_3' )

			set graphics on

			graph combine summclust_temp_ln.gph summclust_temp_pn.gph summclust_temp_lb.gph summclust_temp_pb.gph, title("Cluster Specific Statistics For `G' `cluster' Clusters") scheme(s1mono) 
		}
		else if "`nograph'" != "" {
		    
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
