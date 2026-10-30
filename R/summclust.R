# Independent base-R implementation of the documented summclust mathematics.
# Baseline: Matthew D. Webb, summclust 5.1.1, presentation revision 30sep2026.
# See README.md for numerical conventions, attribution and scope.

.sc_rank <- function(a) {
  d <- svd(a, nu = 0L, nv = 0L)$d
  if (!length(d) || d[1L] == 0) return(0L)
  sum(d > min(dim(a)) * d[1L] * .Machine$double.eps)
}

# Accurately add rounded binary64 products, as quadcross does. The fast branch
# splits a bounded binary lattice into exactly summable high/low pieces. The
# general branch maintains an expansion of nonoverlapping rounding residuals.
.sc_precise_sum <- function(v) {
  if (!length(v)) return(0)
  largest <- max(abs(v))
  if (largest == 0) return(0)
  e <- ceiling(log2(largest)); multiplier <- 2^(52 - e)
  if (is.finite(multiplier) && multiplier > 0 && length(v) < 2^26) {
    scaled <- v * multiplier
    if (all(is.finite(scaled)) && !any(v != 0 & scaled == 0) && all(scaled == trunc(scaled))) {
      unit <- 2^(e - 26)
      high <- trunc(v/unit) * unit
      return(sum(high) + sum(v - high))
    }
  }
  partial <- numeric(0)
  for (value in v) {
    out <- numeric(0)
    for (small in partial) {
      if (abs(value) < abs(small)) { old <- value; value <- small; small <- old }
      high <- value + small
      low <- small - (high - value)
      if (low != 0) out <- c(out, low)
      value <- high
    }
    partial <- c(out, value)
  }
  n <- length(partial); high <- partial[n]; n <- n - 1L; low <- 0
  while (n > 0L) {
    x <- high; y <- partial[n]; n <- n - 1L
    high <- x + y; low <- y - (high - x)
    if (low != 0) break
  }
  if (n > 0L && sign(low) == sign(partial[n])) {
    doubled <- low * 2; rounded <- high + doubled
    if (rounded - high == doubled) high <- rounded
  }
  high
}

.sc_lattice <- function(x) {
  largest <- max(abs(x))
  if (largest == 0) return(c(exponent = -Inf, quantum = 0))
  e <- ceiling(log2(largest))
  if (all(x == trunc(x))) return(c(exponent = e, quantum = 0))
  low <- max(-1074, e - 53); high <- 0
  scale <- 2^(-low)
  if (!is.finite(scale) || any(x * scale != trunc(x * scale)))
    return(c(exponent = e, quantum = -1074))
  while (high - low > 1) {
    mid <- floor((low + high)/2); z <- x * 2^(-mid)
    if (all(z == trunc(z))) low <- mid else high <- mid
  }
  c(exponent = e, quantum = low)
}

.sc_quadcross <- function(x, y = NULL) {
  x <- as.matrix(x); same <- is.null(y)
  if (same) y <- x else y <- as.matrix(y)
  out <- crossprod(x, y)
  xi <- vapply(seq_len(ncol(x)), function(j) .sc_lattice(x[, j]), numeric(2))
  yi <- if (same) xi else vapply(seq_len(ncol(y)), function(j) .sc_lattice(y[, j]), numeric(2))
  nbits <- ceiling(log2(max(1, nrow(x))))
  for (j in seq_len(ncol(x))) for (k in if (same) j:ncol(y) else seq_len(ncol(y))) {
    bound <- nbits + xi[1L, j] + yi[1L, k] - xi[2L, j] - yi[2L, k]
    if (bound > 53) {
      out[j, k] <- .sc_precise_sum(x[, j] * y[, k])
      if (same) out[k, j] <- out[j, k]
    }
  }
  out
}

.sc_mata_dot <- function(x, y) {
  odd <- even <- 0
  for (i in seq_along(x)) {
    if (i %% 2L) odd <- odd + x[i] * y[i] else even <- even + x[i] * y[i]
  }
  odd + even
}

# Column-subset generalized inverse, never a Moore-Penrose inverse. The pivot
# convention is separately validated against Mata on singular fixtures.
.sc_invsym <- function(a) {
  k <- nrow(a)
  ans <- matrix(0, k, k)
  if (!k) return(ans)
  a <- (a + t(a))/2
  scale <- sqrt(pmax(diag(a), 0))
  keep <- which(scale > 0)
  if (!length(keep)) return(ans)
  z <- a[keep, keep, drop = FALSE]
  remaining <- seq_along(keep)
  selected <- integer()
  while (length(remaining)) {
    j <- remaining[which.max(diag(z)[remaining])]
    pivot <- z[j, j]
    remaining <- setdiff(remaining, j)
    if (pivot <= 1e-9 * a[keep[j], keep[j]]) next
    selected <- c(selected, j)
    if (length(remaining)) {
      v <- z[remaining, j]
      z[remaining, remaining] <- z[remaining, remaining, drop = FALSE] - outer(v/pivot, v)
      z[lower.tri(z)] <- t(z)[lower.tri(z)]
    }
  }
  ix <- keep[selected]
  if (length(ix)) {
    aa <- a[ix, ix, drop = FALSE]/outer(scale[ix], scale[ix])
    inv <- tryCatch(chol2inv(chol(aa)), error = function(e) solve(aa))
    ans[ix, ix] <- inv/outer(scale[ix], scale[ix])
  }
  ans
}

.sc_cholsolve <- function(a, b) {
  tryCatch({
    k <- nrow(a); z <- a; l <- matrix(0, k, k)
    for (j in seq_len(k)) {
      if (z[j, j] <= 0) stop("singular")
      l[j, j] <- sqrt(z[j, j])
      if (j < k) {
        ii <- (j + 1L):k
        v <- z[ii, j] * (1/l[j, j])
        l[ii, j] <- v
        z[ii, ii] <- z[ii, ii, drop = FALSE] - tcrossprod(v)
      }
    }
    if (any(diag(l) <= 1e-13 * mean(abs(diag(l))))) stop("singular")
    ans <- as.matrix(b)
    for (column in seq_len(ncol(ans))) {
      for (i in seq_len(k)) {
        jj <- seq_len(i - 1L)
        ans[i, column] <- (ans[i, column] - .sc_mata_dot(l[i, jj], ans[jj, column]))/l[i, i]
      }
      for (i in k:1L) {
        jj <- if (i < k) (i + 1L):k else integer()
        ans[i, column] <- (ans[i, column] - .sc_mata_dot(l[jj, i], ans[jj, column]))/l[i, i]
      }
    }
    ans
  }, error = function(e) .sc_invsym(a) %*% b)
}

.sc_levels <- function(x) sort(unique(x), na.last = NA, method = "radix")
.sc_nested <- function(group, cluster) {
  all(vapply(split(cluster, match(group, .sc_levels(group))), function(z) length(unique(z)) == 1L, logical(1)))
}
.sc_cv <- function(x, tw = FALSE) {
  m <- mean(x)
  if (length(x) < 2L || !is.finite(m) || (tw && abs(m) < 1e-14) || m == 0) return(NA_real_)
  sqrt(sum((x - m)^2)/((length(x) - 1) * m^2))
}
.sc_gstar <- function(x, tw = FALSE) {
  m <- mean(x)
  if (!is.finite(m) || m == 0 || (tw && abs(m) < 1e-14)) return(NA_real_)
  length(x)/(1 + mean(((x - m)/m)^2))
}
.sc_inference <- function(beta, se, df) {
  tt <- beta/se
  tt[!is.finite(tt)] <- NA_real_
  if (!is.finite(df) || df <= 0) {
    return(cbind(beta = beta, se = se, t = tt, p = NA_real_, lci = NA_real_, uci = NA_real_))
  }
  critical <- stats::qt(.975, df)
  cbind(beta = beta, se = se, t = tt, p = 2 * stats::pt(-abs(tt), df),
        lci = beta - critical * se, uci = beta + critical * se)
}
.sc_safe_sqrt <- function(x) sqrt(replace(x, x < 0, NA_real_))
.sc_float32 <- function(x) readBin(writeBin(as.double(x), raw(), size = 4L), "double", n = length(x), size = 4L)
# Mata strofreal() defaults to %9.0g. This is an interface-compatibility layer;
# estimates, covariance matrices, and structured inference remain full precision.
.sc_legacy_number <- function(x) {
  vapply(x, function(v) {
    if (!is.finite(v) || v == 0) return(v)
    av <- abs(v)
    if (av < 1e-5 || av >= 1e7) return(as.numeric(sprintf("%.2e", v)))
    digits <- if (av < 1) 7L else max(0L, 6L - floor(log10(av)))
    ans <- as.numeric(sprintf(paste0("%.", digits, "f"), v))
    if (abs(ans) >= 1e7) as.numeric(sprintf("%.2e", v)) else ans
  }, numeric(1), USE.NAMES = FALSE)
}
.sc_legacy_inference <- function(beta, se, df) {
  ans <- .sc_inference(.sc_legacy_number(beta), .sc_legacy_number(se), df)
  if (df > 0) {
    p <- stats::pt(ans[, "t"], df, lower.tail = FALSE)
    ans[, "p"] <- 2 * pmin(p, 1 - p)
  }
  ans
}
.sc_cond <- function(a) {
  s <- svd(a, nu = 0L, nv = 0L)$d
  if (!length(s) || min(s) == 0) Inf else max(s)/min(s)
}
.sc_identified <- function(a, ai = .sc_invsym(a)) {
  max(abs(drop(ai[1L, , drop = FALSE] %*% a) - c(1, rep(0, ncol(a) - 1)))) < 1e-7
}

.sc_component <- function(x, y, id, ids, a, ay, ai, beta, u, partial,
                          weights, coeffnum, fullrank, rankcheck, forceginv,
                          tw = FALSE, dimension = "G", delete_one = TRUE) {
  k <- ncol(x); g <- length(ids)
  omitted <- matrix(NA_real_, k, g)
  ng <- lever <- partlev <- gamzero <- gamone <- numeric(g)
  legacy <- if (delete_one) integer(g) else rep(NA_real_, g)
  rankfail <- rep(NA_real_, g); delete_one_solves <- 0L
  meat <- matrix(0, k, k)
  indices <- split(seq_len(nrow(x)), id)
  partial_total <- sum(partial^2)
  for (j in seq_len(g)) {
    take <- indices[[as.character(j)]]
    xx <- x[take, , drop = FALSE]
    aa <- .sc_quadcross(xx)
    left <- a - aa
    # CV31 needs intersection scores and cheap diagnostics, but no deleted
    # coefficients. Count actual coefficient solves rather than infer savings.
    if (delete_one) {
      rhs <- ay - .sc_quadcross(xx, y[take])
      if (tw && !.sc_identified(left)) {
        stop(sprintf("Focal coefficient is not numerically identified after deleting %s cluster %s; CV3 is undefined.", dimension, as.character(ids[j])), call. = FALSE)
      }
      delete_one_solves <- delete_one_solves + 1L
      if (forceginv || (k <= 50L && .sc_cond(left) > 1e8)) {
        bg <- .sc_invsym(left) %*% rhs
      } else bg <- .sc_cholsolve(left, rhs)
      omitted[, j] <- bg
      legacy[j] <- as.integer(any(abs(bg) <= 1e-8))
    }
    score <- crossprod(xx, u[take])
    meat <- meat + tcrossprod(score)
    ng[j] <- length(take)
    lever[j] <- sum(aa * t(ai))
    partlev[j] <- sum(partial[take]^2)/partial_total
    gamzero[j] <- sum(weights[take]^2)
    gamone[j] <- sum(weights[take])^2
    if (rankcheck) rankfail[j] <- as.integer(.sc_rank(left) < fullrank)
  }
  delta <- beta - omitted
  v1 <- (g * (nrow(x) - 1)/((g - 1) * (nrow(x) - coeffnum))) * ai %*% meat %*% ai
  v3 <- if (delete_one) (g - 1)/g * tcrossprod(delta) else matrix(NA_real_, k, k)
  list(ids = ids, ng = ng, lever = lever, partlev = partlev,
       betanog = omitted[1L, ], omitted = omitted, legacyflag = legacy,
       rankfail = rankfail, cv1 = v1, cv3 = v3, gamzero = gamzero,
       gamone = gamone, G = g, delete_one_solves = delete_one_solves)
}

#' Cluster diagnostics and one-/two-way cluster-robust covariance estimates.
#' data is a data.frame or named list; outcome/focal/controls/cluster/fevar and
#' absorb name columns. sample is a logical mask of original data rows.
#' For the full contract, read R/README.md.
summclust <- function(data, outcome, focal, controls, cluster, fevar = NULL,
                      absorb = NULL, sample = NULL, gstar = FALSE, rho = NULL,
                      jackknife = FALSE, rankcheck = FALSE, table = FALSE,
                      addmeans = FALSE, regtable = FALSE, mixed = FALSE) {
  if (!is.logical(mixed) || length(mixed) != 1L || is.na(mixed))
    stop("mixed must be a single nonmissing logical value.", call. = FALSE)
  if (length(outcome) != 1L || length(focal) != 1L || !length(controls))
    stop("Supply one outcome, one focal regressor, and at least one control, as required by Stata summclust.", call. = FALSE)
  if (!length(cluster) %in% 1:2 || length(absorb) > 1L)
    stop("cluster must name one or two columns; absorb at most one.", call. = FALSE)
  tw <- length(cluster) == 2L
  if (mixed && !tw) stop("mixed is a two-way option; supply two clustering dimensions.", call. = FALSE)
  if (tw && jackknife) stop("Two-way CV3J is not supported; use jackknife = FALSE.", call. = FALSE)
  if (!is.null(rho) && (length(rho) != 1L || !is.finite(rho) || rho < .00001 || rho > 1))
    stop("rho must lie between 0.00001 and 1.", call. = FALSE)
  dogstar <- gstar || !is.null(rho)
  rho_used <- if (is.null(rho)) .1234567 else rho
  columns <- unique(c(outcome, focal, controls, cluster, fevar, absorb))
  if (any(!columns %in% names(data))) stop("Missing data columns: ", paste(setdiff(columns, names(data)), collapse = ", "), call. = FALSE)
  n0 <- length(data[[outcome]])
  if (any(vapply(data[columns], length, integer(1)) != n0)) stop("All selected columns must have the same length.", call. = FALSE)
  for (nm in c(outcome, focal, controls)) if (!is.numeric(data[[nm]])) stop("Regression columns must be numeric: ", nm, call. = FALSE)
  if (!tw && any(!vapply(data[fevar], is.numeric, logical(1)))) stop("One-way fevar columns must be numeric, matching the current Stata path.", call. = FALSE)
  mask <- rep(TRUE, n0)
  for (nm in columns) {
    z <- data[[nm]]
    mask <- mask & !is.na(z)
    if (is.numeric(z)) mask <- mask & is.finite(z)
    if (is.character(z) || is.factor(z)) mask <- mask & as.character(z) != ""
  }
  if (!is.null(sample)) {
    if (!is.logical(sample) || length(sample) != n0) stop("sample must be a logical vector with one entry per input row.", call. = FALSE)
    mask <- mask & !is.na(sample) & sample
  }
  if (!any(mask)) stop("No observations in the estimation sample.", call. = FALSE)
  d <- as.data.frame(lapply(data[columns], function(z) if (is.factor(z)) as.character(z[mask]) else z[mask]), stringsAsFactors = FALSE)
  if (length(absorb)) for (cl in cluster) {
    if (!.sc_nested(d[[absorb]], d[[cl]])) stop("absorb groups must be nested within every clustering dimension; use fevar for crossed effects.", call. = FALSE)
  }
  # Replicate source row ordering, while keeping the returned original-row mask.
  ord <- if (tw) {
    if (length(absorb)) order(d[[absorb]], seq_len(nrow(d)), method = "radix") else seq_len(nrow(d))
  } else {
    if (length(absorb)) order(d[[cluster]], d[[absorb]], seq_len(nrow(d)), method = "radix") else order(d[[cluster]], seq_len(nrow(d)), method = "radix")
  }
  d <- d[ord, , drop = FALSE]
  fe_levels <- lapply(fevar, function(nm) .sc_levels(d[[nm]]))
  fe_counts <- lengths(fe_levels) - as.integer(length(absorb) > 0L | seq_along(fevar) > 1L)
  has_constant <- !length(fevar) && !length(absorb)
  coefficient_names <- c(focal, if (has_constant) "_cons", controls)
  x <- matrix(0, nrow(d), 1L + has_constant + length(controls) + sum(fe_counts))
  x[, 1L] <- d[[focal]]
  pos <- 1L
  if (has_constant) { pos <- pos + 1L; x[, pos] <- 1 }
  x[, pos + seq_along(controls)] <- as.matrix(d[controls])
  pos <- pos + length(controls)
  if (length(fevar)) for (j in seq_along(fevar)) {
    nm <- fevar[j]; levels <- fe_levels[[j]]
    nkeep <- fe_counts[j]
    if (nkeep > 0L) for (q in seq_len(nkeep)) {
      pos <- pos + 1L
      x[, pos] <- as.numeric(d[[nm]] == levels[q])
      coefficient_names <- c(coefficient_names, paste0(nm, "_", if (tw) q else format(levels[q], trim = TRUE, scientific = FALSE)))
    }
  }
  y <- d[[outcome]]; nabsorb <- 0L
  if (length(absorb)) {
    aa <- match(d[[absorb]], .sc_levels(d[[absorb]])); nabsorb <- max(aa)
    for (j in seq_len(nabsorb)) {
      ix <- which(aa == j)
      x[ix, ] <- sweep(x[ix, , drop = FALSE], 2L, colMeans(x[ix, , drop = FALSE]))
      y[ix] <- y[ix] - mean(y[ix])
    }
  }
  colnames(x) <- coefficient_names
  k <- ncol(x); n <- nrow(x); rank <- .sc_rank(x); coeffnum <- rank + nabsorb
  if (n <= coeffnum) stop("Residual degrees of freedom are not positive after fixed effects.", call. = FALSE)
  ids <- lapply(cluster, function(nm) .sc_levels(d[[nm]]))
  cid <- Map(function(nm, lev) match(d[[nm]], lev), cluster, ids)
  if (any(lengths(ids) < 2L)) stop("At least two nonempty clusters are required in each primary dimension.", call. = FALSE)
  if (tw) {
    pairs <- unique(d[cluster]); pairs <- pairs[do.call(order, c(pairs, list(method = "radix"))), , drop = FALSE]
    key <- paste(cid[[1L]], cid[[2L]], sep = ":")
    pairkey <- paste(match(pairs[[1L]], ids[[1L]]), match(pairs[[2L]], ids[[2L]]), sep = ":")
    cid[[3L]] <- match(key, pairkey)
    ids[[3L]] <- vapply(seq_len(nrow(pairs)), function(i) paste(as.character(unlist(pairs[i, ], use.names = FALSE)), collapse = "|"), character(1))
  }
  a <- .sc_quadcross(x); ay <- .sc_quadcross(x, y); forceginv <- .sc_cond(a) > 1e7
  if (!tw && forceginv) {
    a[,] <- 0; ay[,] <- 0
    for (j in seq_along(ids[[1L]])) {
      ix <- which(cid[[1L]] == j)
      a <- a + .sc_quadcross(x[ix, , drop = FALSE]); ay <- ay + .sc_quadcross(x[ix, , drop = FALSE], y[ix])
    }
  }
  ai <- .sc_invsym(a)
  if (tw && !.sc_identified(a, ai)) stop("Focal coefficient is not numerically identified in the full model.", call. = FALSE)
  beta <- drop(ai %*% ay); names(beta) <- coefficient_names
  u <- drop(y - x %*% beta)
  cc <- x[, -1L, drop = FALSE]
  partial <- drop(x[, 1L] - cc %*% .sc_cholsolve(.sc_quadcross(cc), .sc_quadcross(cc, x[, 1L])))
  weights <- drop(x %*% ai[, 1L])
  components <- lapply(seq_along(cid), function(j) .sc_component(x, y, cid[[j]], ids[[j]], a, ay, ai, beta, u, partial, weights, coeffnum, rank, rankcheck, forceginv, tw, c("G", "H", "intersection")[j], delete_one = !(mixed && j == 3L)))
  for (j in seq_along(components)) {
    dimnames(components[[j]]$cv1) <- dimnames(components[[j]]$cv3) <- list(coefficient_names, coefficient_names)
    dimnames(components[[j]]$omitted) <- list(coefficient_names, as.character(ids[[j]]))
    for (nm in c("ng", "lever", "partlev", "betanog", "rankfail", "legacyflag")) names(components[[j]][[nm]]) <- as.character(ids[[j]])
  }
  names(components) <- c("g", "h", "intersection")[seq_along(components)]
  result <- list(b = matrix(beta, nrow = 1L, dimnames = list(NULL, coefficient_names)),
                 beta = beta[1L], N = n, k = k, rank = rank,
                 coefficient_names = coefficient_names, sample_mask = mask,
                 sample_order = which(mask)[ord], cluster_ids = ids,
                 components = components, absorbed_levels = nabsorb,
                 residual_df = n - coeffnum, rho_used = rho_used,
                 mixed = mixed, selected_estimator = if (mixed) "CV31" else if (jackknife) "CV3J" else "CV3")
  g1ok <- vapply(seq_along(cid), function(j) !length(absorb) && !any(vapply(fevar, function(nm) .sc_nested(d[[nm]], cid[[j]]), logical(1))), logical(1))
  if (!tw) {
    co <- components[[1L]]; g <- co$G; result$G <- g; result$df_r <- g - 1
    for (nm in c("cv1", "cv3", "ng", "lever", "partlev", "betanog")) result[[nm]] <- co[[nm]]
    if (rankcheck) { result$legacyflag <- co$legacyflag; result$rankfail <- co$rankfail }
    if (table && g <= 52L) result$scall <- cbind(Ng = co$ng, Leverage = co$lever, partlev = co$partlev, betanog = co$betanog)
    result$inference <- list(cv1 = .sc_inference(beta, .sc_safe_sqrt(diag(co$cv1)), g - 1), cv3 = .sc_inference(beta, .sc_safe_sqrt(diag(co$cv3)), g - 1))
    if (jackknife) {
      centered <- sweep(co$omitted, 1L, rowMeans(co$omitted))
      result$cv3j <- (g - 1)/g * tcrossprod(centered)
      result$inference$cv3j <- .sc_inference(beta, .sc_safe_sqrt(diag(result$cv3j)), g - 1)
    }
    good <- which(co$legacyflag == 0L)
    result$gdrop <- length(good)
    if (any(co$legacyflag != 0L)) {
      dropmatrix <- function(center) {
        if (!length(good)) return(matrix(NA_real_, k, k, dimnames = list(coefficient_names, coefficient_names)))
        z <- co$omitted[, good, drop = FALSE]
        dev <- sweep(z, 1L, if (center) rowMeans(z) else beta)
        (length(good) - 1)/length(good) * tcrossprod(dev)
      }
      result$inference$cv3drop <- .sc_inference(beta, .sc_safe_sqrt(diag(dropmatrix(FALSE))), length(good) - 1)
      if (regtable) result$cv3drop <- dropmatrix(FALSE)
      if (jackknife) {
        result$inference$cv3jdrop <- .sc_inference(beta, .sc_safe_sqrt(diag(dropmatrix(TRUE))), length(good) - 1)
        if (regtable) result$cv3jdrop <- dropmatrix(TRUE)
      }
    }
    if (dogstar) {
      result$gstarzero <- .sc_gstar(co$gamzero)
      if (g1ok[1L]) {
        result$gstarone <- .sc_gstar(co$gamone)
        result$gstarrho <- .sc_gstar(rho_used * co$gamone + (1 - rho_used) * co$gamzero)
      }
    }
    stats <- list(Ng = co$ng, Leverage = co$lever, partlev = co$partlev, betanog = co$betanog)
    if (length(good) < g && length(good)) stats$betadropg <- co$betanog[good]
    # The full-precision companion retains the source diagnostic's original
    # G-1 denominator, including the selected betadropg subset. It removes
    # only float-storage rounding. A nonzero singleton subset has zero
    # dispersion under this convention; a zero/nonfinite mean stays missing.
    result$clustsum <- vapply(stats, function(z) {
      m <- mean(z)
      cv <- if (!is.finite(m) || m == 0 || g <= 1L) NA_real_ else
        sqrt(sum((z - m)^2)/((g - 1) * m^2))
      c(min(z), stats::quantile(z, .25, type = 2), stats::median(z), m,
        stats::quantile(z, .75, type = 2), max(z), cv)
    }, numeric(7))
    rownames(result$clustsum) <- c("min", "q1", "median", "mean", "q3", "max", "coefvar")
    result$clustsum_full <- result$clustsum
    # Stata's legacy getmata-created temporary columns have float storage.
    result$clustsum <- vapply(stats, function(z) {
      f <- .sc_float32(z); m <- mean(f)
      cv <- if (m == 0) NA_real_ else sqrt(sum((z - m)^2)/((g - 1) * m^2))
      c(min(f), stats::quantile(f, .25, type = 2), stats::median(f), m, stats::quantile(f, .75, type = 2), max(f), cv)
    }, numeric(7))
    rownames(result$clustsum) <- rownames(result$clustsum_full)
    if (addmeans) {
      result$bonus <- vapply(stats, function(z) {
        harmonic <- 1/mean(1/replace(z, z == 0, NA_real_), na.rm = TRUE)
        geometric <- exp(mean(log(replace(z, z <= 0, NA_real_)), na.rm = TRUE))
        quadratic <- sqrt(mean(z^2)); m <- mean(.sc_float32(z))
        c(harmonic, harmonic/m, geometric, geometric/m, quadratic, quadratic/m)
      }, numeric(6))
      result$bonus[1:4, 4:ncol(result$bonus)] <- NA_real_
      rownames(result$bonus) <- c("harmonic", "harmonic_ratio", "geometric", "geometric_ratio", "quadratic", "quadratic_ratio")
    }
    result$beta_full <- unname(beta[1L]); result$beta <- .sc_legacy_number(beta[1L])
    result$cv1se_full <- result$inference$cv1[1L, "se"]
    result$cv3se_full <- result$inference$cv3[1L, "se"]
    if (jackknife) result$cv3Jse_full <- result$inference$cv3j[1L, "se"]
    for (nm in intersect(c("gstarzero", "gstarone", "gstarrho"), names(result))) {
      result[[paste0(nm, "_full")]] <- result[[nm]]
      result[[nm]] <- .sc_legacy_number(result[[nm]])
    }
    for (nm in names(result$inference)) {
      prefix <- if (nm == "cv3j") "cv3J" else nm
      z <- .sc_legacy_inference(beta[1L], result$inference[[nm]][1L, "se"], if (grepl("drop$", nm)) length(good) - 1 else g - 1)[1L, ]
      for (s in c("se", "t", "p", "lci", "uci")) result[[paste0(prefix, s)]] <- unname(z[s])
    }
    # The source's historical macro names are overwritten by CV3J drop when
    # jackknife is requested. Separate inference entries retain both meanings.
    if ("cv3drop" %in% names(result$inference)) {
      z <- .sc_legacy_inference(beta[1L], result$inference[[if (jackknife) "cv3jdrop" else "cv3drop"]][1L, "se"], length(good) - 1)[1L, ]
      for (s in c("se", "t", "p", "lci", "uci")) result[[paste0("cv3", s, "drop")]] <- unname(z[s])
    }
  } else {
    g <- components$g$G; h <- components$h$G; ii <- components$intersection$G
    result$G <- result$N_clust_g <- g; result$H <- result$N_clust_h <- h
    result$I <- result$N_clust_intersection <- ii; result$df_r <- min(g, h) - 1L
    result$max_tie_tolerance <- 1e-10; result$cluster_g <- cluster[1L]; result$cluster_h <- cluster[2L]
    result$cluster_pairs <- pairs
    result$inference <- list()
    for (nm in c("cv1", if (mixed) "cv31" else "cv3")) {
      primary <- if (nm == "cv31") "cv3" else nm
      intersection <- if (nm == "cv31") "cv1" else nm
      vv <- components$g[[primary]] + components$h[[primary]] - components$intersection[[intersection]]
      result[[nm]] <- result[[paste0(nm, "_conventional")]] <- vv
      options <- c(vv[1L, 1L], components$g[[primary]][1L, 1L], components$h[[primary]][1L, 1L])
      chosen <- which.max(options); mx <- options[chosen]
      result[[paste0(nm, "_conventional_var")]] <- vv[1L, 1L]
      result[[paste0(nm, "_conventional_se")]] <- .sc_safe_sqrt(vv[1L, 1L])
      result[[paste0(nm, "_max_var")]] <- result[[paste0(nm, "var_max")]] <- mx
      result[[paste0(nm, "_max_se")]] <- result[[paste0(nm, "se_max")]] <- .sc_safe_sqrt(mx)
      result[[paste0(nm, "max_component")]] <- c("two-way", "G", "H")[chosen]
      result[[paste0(nm, "max_ties")]] <- sum(c(1L, 2L, 4L) * (abs(options - mx) <= 1e-10 * max(abs(mx), 1e-300)))
      se <- .sc_safe_sqrt(pmax(diag(vv), diag(components$g[[primary]]), diag(components$h[[primary]])))
      result$inference[[nm]] <- .sc_inference(beta, se, result$df_r)
      z <- result$inference[[nm]][1L, ]
      for (s in c("se", "t", "p", "lci", "uci")) result[[paste0(nm, s)]] <- unname(z[s])
      if (nm == "cv31") for (s in c("t", "p", "lci", "uci")) result[[paste0("cv31_max_", s)]] <- unname(z[s])
    }
    if (mixed) {
      result$cv3 <- result$cv3_conventional <- matrix(NA_real_, k, k,
        dimnames = list(coefficient_names, coefficient_names))
      result$inference_conventional <- list(cv31 = .sc_inference(beta,
        .sc_safe_sqrt(diag(result$cv31)), result$df_r))
      result$cv31_label <- "CV31 (mixed)"
      result$cv31_max_label <- "CV31 (mixed, max)"
    }
    result$intersection_delete_one_solves <- components$intersection$delete_one_solves
    result$estimator_label <- if (mixed) "CV31 (mixed)" else "CV3"
    result$inference_label <- if (mixed) "CV31 (mixed, max)" else "CV3 (max)"
    result$dimsummary <- matrix(NA_real_, 3L, 6L, dimnames = list(c("G", "H", "I"), c("cv_N", "cv_leverage", "cv_partlev", "cv_betanog", "clusters", "Gstar0")))
    for (j in seq_along(components)) {
      co <- components[[j]]; suf <- names(components)[j]
      result[[paste0("delete_one_solves_", suf)]] <- co$delete_one_solves
      for (nm in c("cv1", "cv3", "ng", "lever", "partlev", "betanog", "rankfail")) result[[paste0(nm, "_", suf)]] <- co[[nm]]
      result[[paste0("gstarzero_", suf)]] <- .sc_gstar(co$gamzero, TRUE)
      if (dogstar && g1ok[j]) {
        result[[paste0("gstarone_", suf)]] <- .sc_gstar(co$gamone, TRUE)
        result[[paste0("gstarrho_", suf)]] <- .sc_gstar(rho_used * co$gamone + (1 - rho_used) * co$gamzero, TRUE)
      }
      result$dimsummary[j, ] <- c(vapply(co[c("ng", "lever", "partlev", "betanog")], .sc_cv, numeric(1), tw = TRUE), co$G, .sc_gstar(co$gamzero, TRUE))
      for (nm in c("cv1", "cv3")) result[[paste0("V_", nm, "_", if (j == 3L) "i" else suf)]] <- co[[nm]]
      if (j == 3L) for (nm in c("ng", "lever", "partlev", "betanog", "rankfail", "gstarzero", "gstarone", "gstarrho")) {
        long <- paste0(nm, "_intersection")
        if (!is.null(result[[long]])) result[[paste0(nm, "_i")]] <- result[[long]]
      }
      if (table && j <= 2L) result[[paste0("table_", suf)]] <- cbind(N = co$ng, leverage = co$lever, partlev = co$partlev, betanog = co$betanog)
    }
    result$vce <- if (mixed) "CV31 (mixed, max)" else "two-way max-se"
  }
  if (regtable) result$regtable <- result$inference
  class(result) <- "summclust"
  result
}

print.summclust <- function(x, ...) {
  cat("summclust: ", x$N, " observations; ", x$G, " clusters", sep = "")
  if (!is.null(x$H)) cat(" by ", x$H, "; ", x$I, " intersections", sep = "")
  cat("; inference df = ", x$df_r, "\n", sep = "")
  rows <- lapply(x$inference[intersect(c("cv1", "cv3", "cv3j", "cv31"), names(x$inference))], function(z) z[1L, ])
  if ("cv31" %in% names(rows)) names(rows)[names(rows) == "cv31"] <- "CV31 (mixed, max)"
  print(do.call(rbind, rows))
  invisible(x)
}

# Plain long-format output for cross-language comparisons. Matrices carry
# coefficient names, and diagnostic vectors carry actual cluster identifiers.
flatten_summclust <- function(result, file = NULL) {
  rows <- list(); pos <- 0L
  excluded <- c("sample_mask", "sample_order", "coefficient_names", "cluster_ids", "cluster_pairs", "components", "inference", "regtable")
  for (nm in setdiff(names(result), excluded)) {
    z <- result[[nm]]
    if (!is.numeric(z) && !is.logical(z)) next
    if (is.matrix(z)) {
      rn <- rownames(z); cn <- colnames(z)
      if (is.null(rn)) rn <- as.character(seq_len(nrow(z)))
      if (is.null(cn)) cn <- as.character(seq_len(ncol(z)))
      for (i in seq_len(nrow(z))) for (j in seq_len(ncol(z))) {
        pos <- pos + 1L; rows[[pos]] <- data.frame(statistic = nm, row = i, column = j, value = z[i, j], row_label = rn[i], column_label = cn[j], stringsAsFactors = FALSE)
      }
    } else {
      rn <- names(z); if (length(z) == 1L) rn <- ""
      if (is.null(rn)) rn <- as.character(seq_along(z))
      for (i in seq_along(z)) { pos <- pos + 1L; rows[[pos]] <- data.frame(statistic = nm, row = i, column = 1L, value = z[i], row_label = rn[i], column_label = "", stringsAsFactors = FALSE) }
    }
  }
  out <- do.call(rbind, rows); rownames(out) <- NULL
  if (!is.null(file)) {
    old <- options(digits = 17); on.exit(options(old), add = TRUE)
    utils::write.csv(out, file, row.names = FALSE, na = "")
  }
  out
}
