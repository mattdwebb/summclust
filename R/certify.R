# Reproducible compact certification with captured process output and status.
args <- commandArgs(trailingOnly=FALSE)
script <- sub("^--file=","",args[grepl("^--file=",args)])
here <- dirname(normalizePath(script[1L],winslash="/"))
root <- dirname(here)
dir.create(file.path(root,"output"),showWarnings=FALSE)
rs <- file.path(R.home("bin"),"Rscript.exe")
if(!file.exists(rs)) rs <- file.path(R.home("bin"),"Rscript")
log <- character()
for(file in c("test_summclust.R","test_mixed.R","test_default_parity.R")) {
  result <- system2(rs,c("--vanilla",shQuote(file.path(here,file))),stdout=TRUE,stderr=TRUE)
  code <- attr(result,"status"); if(is.null(code)) code <- 0L
  log <- c(log,paste("SCRIPT",file,"EXIT",code),result,"")
  if(code!=0L) {
    writeLines(log,file.path(root,"output","r_certification.log"))
    stop("Certification failed: ",file)
  }
}
writeLines(c(log,"R CERTIFICATION PASSED"),file.path(root,"output","r_certification.log"))
cat("R CERTIFICATION PASSED; results saved to output/r_certification.log\n")
