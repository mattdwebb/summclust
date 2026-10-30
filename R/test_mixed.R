argv <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", argv[grepl("^--file=", argv)])
here <- dirname(normalizePath(script[1L], winslash = "/"))
source(file.path(here, "summclust.R"))
checks <- 0L
check <- function(ok, label) {
  if (!isTRUE(ok)) stop("FAILED: ", label, call. = FALSE)
  checks <<- checks + 1L
}
near <- function(a, b, tol = 2e-10) max(abs(a-b)) < tol
fails <- function(expr) inherits(tryCatch(force(expr), error = identity), "error")
i <- seq_len(420L)
d <- data.frame(g = (i-1L) %/% 35L, h = (i-1L) %% 7L,
  x = sin(i*.713)+i/191, z = cos(i*.391)+i/193, w = sin(i*.213),
  fe = (i-1L) %% 4L)
d$y <- 1.7 + .8*d$x - .4*d$z + .2*d$w + .3*sin(i*.941)
d$a <- paste(d$g, d$h, sep = ":")
base <- summclust(d, "y", "x", c("z", "w"), c("g", "h"))
off <- summclust(d, "y", "x", c("z", "w"), c("g", "h"), mixed=FALSE)
check(identical(base,off), "omitted and FALSE options identical")
m <- summclust(d, "y", "x", c("z", "w"), c("g", "h"), mixed=TRUE)
cg <- summclust(d,"y","x",c("z","w"),"g")
ch <- summclust(d,"y","x",c("z","w"),"h")
ci <- summclust(d,"y","x",c("z","w"),"a")
check(near(m$cv31,cg$cv3+ch$cv3-ci$cv1), "independent components construct CV31")
check(identical(m$b,base$b), "coefficient invariance")
check(identical(m$cv3_g,base$cv3_g) && identical(m$cv3_h,base$cv3_h), "primary CV3 invariance")
check(identical(m$cv1,base$cv1) && identical(m$cv1_intersection,base$cv1_intersection), "CV1 unchanged")
check(all(is.na(m$cv3)) && all(is.na(m$cv3_intersection)), "uncomputed CV3 unavailable")
check(all(is.na(m$betanog_intersection)) && all(is.na(m$components$intersection$omitted)), "omitted estimates unavailable")
check(m$intersection_delete_one_solves==0L && base$intersection_delete_one_solves==base$I, "actual I solve counts")
check(identical(m$ng_intersection,base$ng_intersection) && near(m$lever_intersection,base$lever_intersection) && near(m$partlev_intersection,base$partlev_intersection), "cheap diagnostics retained")
check(near(m$gstarzero_intersection,base$gstarzero_intersection), "I Gstar retained")
x <- cbind(d$x,1,d$z,d$w)
fit <- stats::lm.fit(x,d$y)
score <- vapply(sort(unique(d$a)),function(a) drop(crossprod(x[d$a==a,,drop=FALSE],fit$residuals[d$a==a])),numeric(ncol(x)))
ai <- solve(crossprod(x))
direct <- m$I/(m$I-1)*(nrow(d)-1)/(nrow(d)-ncol(x))*ai%*%tcrossprod(score)%*%ai
check(near(m$cv1_intersection,direct), "CV1 I score correction independent")
vmax <- max(m$cv31[1,1],cg$cv3[1,1],ch$cv3[1,1])
check(near(m$cv31_max_var,vmax) && near(m$cv31_max_se,sqrt(vmax)),"max comparison")
check(m$df_r==min(m$G,m$H)-1,"df")
tt <- m$b[1,1]/sqrt(vmax)
check(near(m$cv31p,2*pt(-abs(tt),m$df_r)),"p value")
check(near(m$cv31lci,m$b[1,1]-qt(.975,m$df_r)*sqrt(vmax)) && near(m$cv31uci,m$b[1,1]+qt(.975,m$df_r)*sqrt(vmax)),"confidence interval")
check(m$selected_estimator=="CV31" && m$inference_label=="CV31 (mixed, max)","labels")
check(any(grepl("CV31 (mixed, max)",capture.output(print(m)),fixed=TRUE)),"printed label")
check(fails(summclust(d,"y","x","z","g",mixed=TRUE)),"one-way mixed rejects")
check(fails(summclust(d,"y","x","z",c("g","h"),mixed=TRUE,jackknife=TRUE)),"CV3J rejects")
for (bad in list(1,"TRUE",NA,c(TRUE,FALSE),NULL)) check(fails(summclust(d,"y","x","z",c("g","h"),mixed=bad)),"mixed type rejects")
same <- summclust(d,"y","x",c("z","w"),c("g","g"),mixed=TRUE)
check(near(same$cv31,2*cg$cv3-cg$cv1),"identical partitions use 2CV3-CV1")
dn <- d; dn$h <- dn$g %/% 3L
nested <- summclust(dn,"y","x",c("z","w"),c("g","h"),mixed=TRUE)
cn <- summclust(dn,"y","x",c("z","w"),"h")
check(near(nested$cv31,cg$cv3+cn$cv3-cg$cv1),"nested identity")
sw <- summclust(d,"y","x",c("z","w"),c("h","g"),mixed=TRUE)
check(near(sw$cv31,m$cv31),"dimension swap")
shuffle <- rev(seq_len(nrow(d)))
sh <- summclust(d[shuffle,],"y","x",c("z","w"),c("g","h"),mixed=TRUE)
check(near(sh$cv31,m$cv31) && near(sh$partlev_intersection[names(m$partlev_intersection)],m$partlev_intersection),"row permutation")
for (mode in c("fe","absorb","missing","focal","duplicate","unbalanced")) {
  dd <- d; extra <- list(); focal <- "x"; controls <- c("z","w")
  if (mode=="fe") extra$fevar <- "fe"
  if (mode=="absorb") extra$absorb <- "a"
  if (mode=="missing") { dd$y[c(2,17)] <- NA; dd$g <- paste0("g",dd$g); dd$g[8] <- "" }
  if (mode=="focal") { focal <- "z"; controls <- c("x","w") }
  if (mode=="duplicate") { dd$zdup <- 2*dd$z; controls <- c(controls,"zdup") }
  if (mode=="unbalanced") dd <- dd[seq_len(nrow(dd))%%11!=0,]
  args <- c(list(data=dd,outcome="y",focal=focal,controls=controls,cluster=c("g","h")),extra)
  ordinary <- do.call(summclust,args)
  mm <- do.call(summclust,c(args,list(mixed=TRUE,rankcheck=TRUE)))
  check(near(mm$cv31,ordinary$cv3_g+ordinary$cv3_h-ordinary$cv1_intersection),paste(mode,"formula"))
  check(mm$intersection_delete_one_solves==0L && all(is.finite(mm$rankfail_intersection)),paste(mode,"rankcheck without solves"))
}
negative <- read.csv(file.path(dirname(here),"validation","fixtures","negative.csv"))
neg <- summclust(negative,"y","x","z",c("g","h"),mixed=TRUE)
check(neg$cv31_conventional_var<0 && is.na(neg$cv31_conventional_se),"negative conventional variance retained")
check(is.na(neg$inference_conventional$cv31[1,"p"]) && is.finite(neg$cv31p),"negative conventional vs valid max inference")
cat("R mixed certification passed: ",checks," checks.\n",sep="")
