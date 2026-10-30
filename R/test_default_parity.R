# Compare every pre-existing return with the untouched baseline, without mixed.
args <- commandArgs(trailingOnly=FALSE)
script <- sub("^--file=","",args[grepl("^--file=",args)])
root <- dirname(dirname(normalizePath(script[1L],winslash="/")))
old <- new.env(); new <- new.env()
sys.source(file.path(root,"tests","baseline","R","summclust.R"),old)
sys.source(file.path(root,"R","summclust.R"),new)
cases <- read.csv(file.path(root,"validation","fixtures","cases.csv"),stringsAsFactors=FALSE,na.strings="")
splitnames <- function(x) if(is.na(x)||!nzchar(x)) NULL else strsplit(x,"|",fixed=TRUE)[[1L]]
identical_existing <- function(a,b,path="result") {
  if(is.list(a)) {
    for(nm in names(a)) identical_existing(a[[nm]],b[[nm]],paste(path,nm,sep="$"))
  } else if(!identical(a,b)) stop("Changed existing field: ",path,call.=FALSE)
  invisible(TRUE)
}
for(j in seq_len(nrow(cases))) {
  cc <- cases[j,]
  d <- read.csv(file.path(root,"validation","fixtures",cc$dataset),stringsAsFactors=FALSE)
  call <- list(data=d,outcome=cc$outcome,focal=cc$focal,controls=splitnames(cc$controls),
    cluster=splitnames(cc$cluster),fevar=splitnames(cc$fevar),absorb=splitnames(cc$absorb),
    sample=if(is.na(cc$sample))NULL else d[[cc$sample]]==1L,rho=cc$rho)
  for(nm in c("jackknife","rankcheck","table","addmeans","regtable")) call[[nm]] <- as.logical(cc[[nm]])
  a <- tryCatch(do.call(old$summclust,call),error=identity)
  b <- tryCatch(do.call(new$summclust,call),error=identity)
  if(inherits(a,"error")) {
    stopifnot(inherits(b,"error"),conditionMessage(a)==conditionMessage(b),cc$expected_error!=0L)
  } else {
    stopifnot(!inherits(b,"error"),cc$expected_error==0L)
    identical_existing(a,b)
  }
  cat(cc$specification,": exact baseline parity\n",sep="")
}
cat("R DEFAULT PARITY PASSED: ",nrow(cases)," cases\n",sep="")
