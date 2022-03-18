# summclust
Stata module for cluster level measures of leverage, influence,  and a cluster jackknife variance estimator.

For a very detailed description see:
MacKinnon, J.G., Nielsen, M.Ø., Webb, M.D., 2022. Leverage, influence, and the jackknife
in clustered regression models: Reliable inference using summclust. QED Working Paper 1483. Queen’s University.    

nlswork example - using regress

  . webuse nlswork, clear
  . reg ln_wage i.grade i.age i.birth_yr union race msp, cluster(ind)

nlswork - using summclust

  . summclust msp, yvar(ln_wage) xvar(union race) fevar(grade age birth_yr) cluster(ind)

adding industry fixed effects using absorb

  . summclust msp, yvar(ln_wage) xvar(union race) fevar(grade age birth_yr) absorb(ind) cluster(ind)

sample restrictions - using sample

  . summclust msp, yvar(ln_wage) xvar(union race i.grade) fevar(age birth_yr) sample(age>=25)
            cluster(ind)

Effective Number of Clusters using gstar or rho.

  . summclust msp, yvar(ln_wage) xvar(union race) fevar(grade age birth_yr) cluster(ind) gstar

  . summclust msp, yvar(ln_wage) xvar(union race) fevar(grade age birth_yr) cluster(ind) rho(0.5)

All Output

  . summclust msp, yvar(ln_wage) xvar(union race) fevar(grade age birth_yr) absorb(dubind)
            cluster(ind) table svars jack rho(0.5)

