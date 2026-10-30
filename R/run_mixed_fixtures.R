args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)])
root <- dirname(dirname(normalizePath(script[1L], winslash = "/")))
dir.create(file.path(root,"output"),showWarnings=FALSE)
source(file.path(root, "R", "summclust.R"))
cases <- read.csv(file.path(root,"validation","fixtures","cases.csv"),stringsAsFactors=FALSE,na.strings="")
cases <- cases[grepl("^two_",cases$specification) & cases$expected_error==0,]
splitnames <- function(x) if (is.na(x)||!nzchar(x)) NULL else strsplit(x,"|",fixed=TRUE)[[1L]]
values <- list(); keys <- list()
for (j in seq_len(nrow(cases))) {
  cc <- cases[j,]
  d <- read.csv(file.path(root,"validation","fixtures",cc$dataset),stringsAsFactors=FALSE)
  call <- list(data=d,outcome=cc$outcome,focal=cc$focal,controls=splitnames(cc$controls),
    cluster=splitnames(cc$cluster),fevar=splitnames(cc$fevar),absorb=splitnames(cc$absorb),
    sample=if(is.na(cc$sample)) NULL else d[[cc$sample]]==1L,rho=cc$rho,mixed=TRUE)
  for (nm in c("jackknife","rankcheck","table","addmeans","regtable")) call[[nm]] <- as.logical(cc[[nm]])
  fit <- do.call(summclust,call)
  stopifnot(fit$intersection_delete_one_solves==0L)
  flat <- flatten_summclust(fit)
  flat$specification <- cc$specification
  values[[length(values)+1L]] <- flat[,c("specification","statistic","row","column","value","row_label","column_label")]
  for (m in seq_along(fit$cluster_ids)) {
    first <- if(m==3L) fit$cluster_pairs[[1L]] else fit$cluster_ids[[m]]
    second <- if(m==3L) fit$cluster_pairs[[2L]] else rep("",length(first))
    keys[[length(keys)+1L]] <- data.frame(specification=cc$specification,
      dimension=c("g","h","i")[m],index=seq_along(first),
      identifier1=as.character(first),identifier2=as.character(second))
  }
  cat(cc$specification,": CV31 var=",fit$cv31_conventional_var,"; I solves=0\n",sep="")
}
options(digits=17)
write.csv(do.call(rbind,values),file.path(root,"output","r_mixed.csv"),row.names=FALSE,na="")
write.csv(do.call(rbind,keys),file.path(root,"output","r_mixed_keys.csv"),row.names=FALSE,na="")
cat("R MIXED FIXTURES PASSED: ",nrow(cases)," cases\n",sep="")
