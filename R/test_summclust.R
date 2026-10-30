argv <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", argv[grepl("^--file=", argv)])
here <- if (length(script)) dirname(normalizePath(script[1L], winslash = "/")) else "R"
source(file.path(here, "summclust.R"))
checks <- 0L
check <- function(ok, label) {
  if (!isTRUE(ok)) stop("FAILED: ", label, call. = FALSE)
  checks <<- checks + 1L
}
close_enough <- function(a, b, tol = 2e-10) max(abs(a - b)) < tol
fails <- function(expr) inherits(tryCatch(force(expr), error = identity), "error")

i <- 1:240
d <- data.frame(x = sin(i*.713) + i/191, z = cos(i*.391) + i/193,
                w = sin(i*.213), g = rep(c(2, 5, 11, 20, 22, 37, 41, 88), each = 30),
                h = rep(c("north", "south", "west", "east", "centre"), 48),
                fe = rep(c(3, 7, 12, 20), 60))
d$y <- 1.7 + .8*d$x - .4*d$z + .2*d$w + .3*sin(i*.941)
d$a <- paste(d$g, d$h)
one <- summclust(d, "y", "x", c("z", "w"), "g", jackknife = TRUE, rho = .5, rankcheck = TRUE, table = TRUE, addmeans = TRUE)
x <- cbind(d$x, 1, d$z, d$w)
qrfit <- stats::lm.fit(x, d$y)
check(close_enough(drop(one$b), qrfit$coefficients), "OLS compared with independent QR fit")
check(close_enough(sum(one$lever), ncol(x)), "leverage sums to model rank")
check(close_enough(sum(one$partlev), 1), "partial leverage sums to one")
pr <- stats::lm.fit(x[, -1L], x[, 1L])$residuals
independent_partial <- vapply(sort(unique(d$g)), function(g) sum(pr[d$g == g]^2)/sum(pr^2), numeric(1))
check(close_enough(one$partlev, independent_partial), "focal residualization compared with QR")
check(close_enough(one$gstarzero_full, 1/sum(independent_partial^2)), "Gstar0 identity")
deleted <- vapply(sort(unique(d$g)), function(g) stats::lm.fit(x[d$g != g, ], d$y[d$g != g])$coefficients, numeric(ncol(x)))
deviation <- qrfit$coefficients - deleted
cv3 <- (one$G - 1)/one$G * tcrossprod(deviation)
check(close_enough(one$cv3, cv3), "CV3 compared with separate QR deletion refits")
check(close_enough(one$betanog, deleted[1L, ]), "every deleted focal coefficient")
centered <- sweep(deleted, 1L, rowMeans(deleted))
check(close_enough(one$cv3j, (one$G - 1)/one$G * tcrossprod(centered)), "centered CV3J")
score <- vapply(sort(unique(d$g)), function(g) drop(crossprod(x[d$g == g, ], qrfit$residuals[d$g == g])), numeric(ncol(x)))
inverse <- solve(crossprod(x))
cv1 <- one$G/(one$G - 1) * (nrow(d) - 1)/(nrow(d) - ncol(x)) * inverse %*% tcrossprod(score) %*% inverse
check(close_enough(one$cv1, cv1), "CV1 independent sandwich")
check(close_enough(one$inference$cv1[, "p"], 2*stats::pt(-abs(qrfit$coefficients/sqrt(diag(cv1))), one$G - 1)), "Student-t inference")

shuffle <- c(seq(2, nrow(d), by = 2), seq(1, nrow(d), by = 2))
permuted <- summclust(d[shuffle, ], "y", "x", c("z", "w"), "g", jackknife = TRUE)
check(close_enough(one$cv3, permuted$cv3), "row permutation invariance")
check(close_enough(one$partlev, permuted$partlev), "label-aligned diagnostic invariance")

absorbed <- summclust(d, "y", "x", c("z", "w"), "g", absorb = "g", jackknife = TRUE)
explicit <- summclust(d, "y", "x", c("z", "w"), "g", fevar = "g", jackknife = TRUE)
check(close_enough(absorbed$b[, 1:3], explicit$b[, 1:3]), "absorption vs explicit fixed-effect slopes")
check(close_enough(absorbed$cv1, explicit$cv1[1:3, 1:3]), "absorbed rank finite-sample correction")
check(close_enough(absorbed$cv3, explicit$cv3[1:3, 1:3]), "absorption preserves deletion slopes")

tw <- summclust(d, "y", "x", c("z", "w"), c("g", "h"), rho = .5, rankcheck = TRUE)
check(close_enough(tw$cv1, tw$cv1_g + tw$cv1_h - tw$cv1_intersection), "two-way CV1 inclusion-exclusion")
check(close_enough(tw$cv3, tw$cv3_g + tw$cv3_h - tw$cv3_intersection), "two-way CV3 inclusion-exclusion")
check(tw$df_r == min(tw$G, tw$H) - 1, "two-way inference degrees of freedom")
check(close_enough(tw$cv3_max_var, max(tw$cv3[1, 1], tw$cv3_g[1, 1], tw$cv3_h[1, 1])), "coefficient-specific max rule")
same <- summclust(d, "y", "x", c("z", "w"), c("g", "g"))
check(close_enough(same$cv3, one$cv3), "identical dimensions reduce to one-way")
check(same$cv3max_ties == 7L, "max tie mask")
check(fails(summclust(d, "y", "x", c("z", "w"), c("g", "h"), jackknife = TRUE)), "two-way CV3J rejects")
check(fails(summclust(d, "y", "x", c("z", "w"), c("g", "h"), absorb = "g")), "crossed absorption rejects")

dd <- d; dd$missing_y <- dd$y; dd$missing_y[c(2, 7, 18)] <- NA
dd$gs <- paste0("cluster_", dd$g); dd$gs[29] <- ""
mask <- rep(c(TRUE, TRUE, FALSE), 80)
mf <- summclust(dd, "missing_y", "x", c("z", "w"), "gs", sample = mask)
expected <- mask & !is.na(dd$missing_y) & dd$gs != ""
check(identical(mf$sample_mask, expected), "missing values, empty strings and sample mask")

dd$zdup <- 2*dd$z
singular <- summclust(dd, "y", "x", c("z", "zdup"), "g", rankcheck = TRUE)
check(all(singular$b[, 3] == 0), "scaled duplicate omits smaller pivot column")
check(all(singular$legacyflag == 1L), "legacy zero-coefficient flag is preserved")
check(all(singular$rankfail == 0L), "legacy proxy distinguished from actual rank loss")
check(fails(summclust(dd, "y", "z", c("zdup", "w"), c("g", "h"))), "unidentified focal slope rejects in two-way")
a <- matrix(c(1, 2, 2, 4), 2L)
ai <- .sc_invsym(a)
check(close_enough(a %*% ai %*% a, a) && close_enough(ai %*% a %*% ai, ai), "generalized-inverse identities")
check(close_enough(.sc_invsym(diag(c(1, 1e-17))), diag(c(1, 1e17)), 1), "scale-relative inverse retains small independent column")

check(identical(.sc_precise_sum(c(1e16, 1, -1e16)), 1), "compensated sum preserves cancellation residual")
check(identical(.sc_precise_sum(c(2^900, 2^-900, -2^900)), 2^-900), "extreme exponent range preserves tiny cancellation residual")
check(identical(.sc_mata_dot(c(1e16, 1, -1e16, 1), rep(1, 4)), 2), "two-accumulator dot operation order")

# Two localized nuisance columns disappear under different deletions. The
# third deletion leaves one full-rank selected fit, exposing the diagnostic
# singleton subset without changing the estimator or using stored outputs.
ii <- seq_len(90L)
small <- data.frame(g = rep(1:3, each = 30L), x = sin(ii*.73),
                    z1 = as.numeric(ii <= 30L), z2 = ifelse(ii > 30L & ii <= 60L, cos(ii*.41), 0))
small$y <- 2 + .7*small$x + .5*small$z1 + .3*small$z2 + .1*sin(ii*.17)
singleton <- summclust(small, "y", "x", c("z1", "z2"), "g", rankcheck = TRUE)
selected <- singleton$betanog[singleton$legacyflag == 0L]
check(length(selected) == 1L && is.finite(selected) && selected != 0,
      "diagnostic fixture has exactly one nonzero selected deletion")
direct_cv <- sqrt(sum((selected - mean(selected))^2)/((singleton$G - 1)*mean(selected)^2))
check(identical(unname(singleton$clustsum_full[7L, 5L]), direct_cv) && direct_cv == 0,
      "unrounded singleton dropped-subset dispersion retains original G denominator")
selected <- singleton$betanog
direct_cv <- sqrt(sum((selected - mean(selected))^2)/((singleton$G - 1)*mean(selected)^2))
check(close_enough(singleton$clustsum_full[7L, 4L], direct_cv),
      "full-precision diagnostic agrees with independent original-G formula")
cat("Base-R independent tests passed: ", checks, " checks.\n", sep = "")
