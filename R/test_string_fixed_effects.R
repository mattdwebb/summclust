argv <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", argv[grepl("^--file=", argv)])
here <- if (length(script)) dirname(normalizePath(script[1L], winslash = "/")) else "R"
source(file.path(here, "summclust.R"))
checks <- 0L
check <- function(ok, label) {
  if (!isTRUE(ok)) stop("FAILED: ", label, call. = FALSE)
  checks <<- checks + 1L
}

# Identical group memberships with labels that cannot serve as coefficient names.
# Every normal level remains observed after exclusions; the extra levels occur
# only in rows outside the requested estimation sample.
i <- seq_len(480L)
numeric_data <- data.frame(
  g = rep(seq_len(12L), each = 40L), h = rep(seq_len(8L), 60L),
  f = 1L + (i - 1L) %/% 8L %% 4L, f2 = 1L + (i - 1L) %/% 3L %% 3L,
  x = sin(i*.731) + i/563, z = cos(i*.419) + i/571,
  w = sin(i*.233))
numeric_data$a <- (numeric_data$g - 1L)*8L + numeric_data$h
numeric_data$y <- .8*numeric_data$x - .3*numeric_data$z + .2*numeric_data$w +
  .15*numeric_data$f + .09*numeric_data$f2 + cos(i*.913)
numeric_data$gs <- sprintf("group %02d", numeric_data$g)
labels_f <- c("A space", "B|pipe", "C\"quote", "D/slash")
labels_f2 <- c("Inf", "NA", "nan")
string_data <- numeric_data
string_data$f <- labels_f[numeric_data$f]
string_data$f2 <- labels_f2[numeric_data$f2]
string_data$a <- sprintf("area %03d", numeric_data$a)
sample <- i %% 7L != 0L
sample[1L] <- FALSE
numeric_data$f[1L] <- numeric_data$f2[1L] <- numeric_data$a[1L] <- 999L
string_data$f[1L] <- string_data$f2[1L] <- string_data$a[1L] <- "ZZ excluded"
numeric_data$y[5L] <- string_data$y[5L] <- NA_real_
numeric_data$f[7L] <- NA_integer_; string_data$f[7L] <- NA_character_
numeric_data$f2[8L] <- NA_integer_; string_data$f2[8L] <- ""
numeric_data$a[9L] <- NA_integer_; string_data$a[9L] <- ""
factor_data <- string_data
for (nm in c("f", "f2", "a")) {
  factor_data[[nm]] <- factor(factor_data[[nm]],
    levels = c("unused first", rev(sort(unique(factor_data[[nm]]), na.last = NA)), "unused last"))
}
original_numeric <- numeric_data
original_string <- string_data
original_factor <- factor_data

specs <- list(
  one_fe = list(cluster = "g", fevar = "f", jackknife = TRUE),
  one_multi_fe = list(cluster = "g", fevar = c("f", "f2")),
  one_string_cluster = list(cluster = "gs", fevar = c("f", "f2")),
  one_absorb = list(cluster = "g", absorb = "a", jackknife = TRUE),
  one_both = list(cluster = "g", fevar = c("f", "f2"), absorb = "a"),
  two_fe = list(cluster = c("g", "h"), fevar = c("f", "f2")),
  two_absorb = list(cluster = c("g", "h"), absorb = "a"),
  two_both = list(cluster = c("g", "h"), fevar = "f2", absorb = "a"),
  mixed_fe = list(cluster = c("g", "h"), fevar = c("f", "f2"), mixed = TRUE),
  mixed_absorb = list(cluster = c("g", "h"), absorb = "a", mixed = TRUE),
  mixed_both = list(cluster = c("g", "h"), fevar = "f2", absorb = "a", mixed = TRUE))
fit <- function(data, spec) do.call(summclust, c(list(data = data,
  outcome = "y", focal = "x", controls = c("z", "w"), sample = sample,
  rho = .5, rankcheck = TRUE, table = TRUE, regtable = TRUE), spec))
for (nm in names(specs)) {
  reference <- fit(numeric_data, specs[[nm]])
  actual <- fit(string_data, specs[[nm]])
  factored <- fit(factor_data, specs[[nm]])
  # Comparing the entire result checks coefficients, full covariances,
  # inference, each deletion, rank/G* diagnostics, row masks, and all labels.
  check(identical(actual, reference), paste(nm, "string/numeric entire-result parity"))
  check(identical(factored, reference), paste(nm, "factor/numeric entire-result parity"))
}
check(identical(numeric_data, original_numeric) && identical(string_data, original_string) &&
      identical(factor_data, original_factor), "input columns, row order, and factor levels unchanged")

# Literal strings that resemble missing/infinite numeric values are categories.
literal <- fit(string_data, specs$one_multi_fe)
check(all(c("f2_1", "f2_2") %in% literal$coefficient_names),
      "literal Inf/NA/nan categories remain observed")
check(literal$N == sum(sample & !is.na(string_data$y) & !is.na(string_data$f) & string_data$f2 != ""),
      "only actual missing values and empty strings are excluded")

# Numeric labels retain their existing one-way names, even when nonconsecutive.
nonconsecutive <- numeric_data
nonconsecutive$f <- c(-2.5, 3.5, 20, 100)[match(nonconsecutive$f, seq_len(4L))]
numeric_fit <- fit(nonconsecutive, specs$one_fe)
check(identical(tail(numeric_fit$coefficient_names, 4L), c("f_-2.5", "f_3.5", "f_20", "f_100")),
      "one-way numeric coefficient names remain based on numeric values")

# Long string-effect names sharing a prefix still have distinct, portable labels.
long_names <- paste0(strrep("q", 31L), c("a", "b"))
long_data <- string_data
names(long_data)[match(c("f", "f2"), names(long_data))] <- long_names
for (cluster in list("g", c("g", "h"))) {
  long_fit <- fit(long_data, list(cluster = cluster, fevar = long_names))
  expected_names <- c(paste0(strrep("q", 27L), "_f1_", 1:4),
                      paste0(strrep("q", 27L), "_f2_", 1:2))
  check(identical(tail(long_fit$coefficient_names, 6L), expected_names) &&
        !anyDuplicated(long_fit$coefficient_names) && all(nchar(long_fit$coefficient_names) <= 32L),
        paste(length(cluster), "way long string-effect labels are distinct and portable"))
  reference <- fit(string_data, list(cluster = cluster, fevar = c("f", "f2")))
  check(identical(unname(long_fit$b), unname(reference$b)) &&
        identical(unname(long_fit$cv1), unname(reference$cv1)) &&
        identical(unname(long_fit$cv3), unname(reference$cv3)),
        paste(length(cluster), "way long names do not change estimates"))
}

# Absorption nesting remains a requirement for every clustering dimension.
bad <- string_data
bad$a <- rep(c("crossed A", "crossed B"), length.out = nrow(bad))
for (cluster in list("g", c("g", "h"))) {
  error <- tryCatch(summclust(bad, "y", "x", c("z", "w"), cluster,
                            absorb = "a", sample = sample), error = identity)
  check(inherits(error, "error") && grepl("nested within every clustering dimension", conditionMessage(error), fixed = TRUE),
        paste(length(cluster), "way crossed string absorption rejected"))
}
cat("R string fixed-effect tests passed: ", checks, " checks.\n", sep = "")
