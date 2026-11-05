args <- commandArgs(trailingOnly=FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)])
root <- dirname(dirname(normalizePath(script[1L], winslash="/")))
source(file.path(root, "R", "summclust.R"))
dir.create(file.path(root, "output"), showWarnings=FALSE)
values <- keys <- list()
for (mode in c("one", "two", "mixed")) for (design in c("fe", "absorb", "both")) {
  case <- paste(mode, design, sep="_")
  d <- read.csv(file.path(root, "validation", "fixtures", "synthetic.csv"), stringsAsFactors=FALSE)
  d$fstr <- sprintf("Effect / %02d", d$fe)
  d$hstr <- sprintf("Period %02d", d$h)
  d$astr <- sprintf("Cell: %04d", d$a)
  d$fstr[d$row %% 19 == 0] <- ""
  d$hstr[d$row %% 23 == 0] <- ""
  d$astr[d$row %% 29 == 0] <- ""
  d$fstr[d$keep == 0] <- "Excluded level"
  fit <- summclust(d, "y", "x", c("z", "w"), if (mode == "one") "gs" else c("gs", "hs"),
    fevar=if (design == "fe") c("fstr", "hstr") else if (design == "both") "fstr" else NULL,
    absorb=if (design != "fe") "astr" else NULL, sample=d$keep == 1,
    rho=.5, rankcheck=TRUE, table=TRUE, jackknife=mode == "one", mixed=mode == "mixed")
  flat <- flatten_summclust(fit)
  flat$specification <- case
  values[[length(values)+1L]] <- flat[, c("specification", "statistic", "row", "column", "value", "row_label", "column_label")]
  for (j in seq_along(fit$cluster_ids)) {
    first <- if (j == 3L) fit$cluster_pairs[[1L]] else fit$cluster_ids[[j]]
    second <- if (j == 3L) fit$cluster_pairs[[2L]] else rep("", length(first))
    keys[[length(keys)+1L]] <- data.frame(specification=case, dimension=c("g", "h", "i")[j],
      index=seq_along(first), identifier1=as.character(first), identifier2=as.character(second))
  }
  cat("PASS R string FE export", case, "\n")
}
options(digits=17)
write.csv(do.call(rbind, values), file.path(root, "output", "r_string_fe.csv"), row.names=FALSE, na="")
write.csv(do.call(rbind, keys), file.path(root, "output", "r_string_fe_keys.csv"), row.names=FALSE, na="")
