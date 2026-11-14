## Split-half stability of a 2-component varimax solution.
##
## Repeatedly splits the respondents into two random halves, redoes the analysis
## on each half (same correlation type as PCA.analysis(): polychoric, or Spearman
## when the codes span more than 8 values; pairwise deletion of NAs), and compares
## the loadings of the two halves with each other using Tucker's congruence
## coefficient. Components are matched between halves (axes may swap) and signs
## are ignored (they are arbitrary).
##
## Returns a list with the B x 2 matrix of congruences and a summary table.

library(psych)

split_half_congruence <- function(w, cols, B = 500, seed = 1) {
  
  dat <- as.data.frame(lapply(w[, cols, drop = FALSE], as.numeric))
  nfac <- 2   # matching below is written for two components
  
  cormat <- function(d) {
    rng <- diff(range(as.matrix(d), na.rm = TRUE)) + 1
    R <- if (rng > 8) {
      cor(d, method = "spearman", use = "pairwise.complete.obs")
    } else {
      suppressWarnings(polychoric(d, global = FALSE)$rho)
    }
    cor.smooth(R)
  }
  
  load_vm <- function(d) {
    unclass(principal(cormat(d), nfactors = nfac, rotate = "varimax")$loadings)[, 1:nfac, drop = FALSE]
  }
  
  safe_load <- function(d) tryCatch(suppressWarnings(load_vm(d)), error = function(e) NULL)
  
  # congruence between two loading matrices (half A vs half B): the two halves
  # are independent, unlike a half compared with the full sample it belongs to
  match_pair <- function(LA, LB) {
    if (is.null(LA) || is.null(LB)) return(c(NA_real_, NA_real_))
    phi <- abs(factor.congruence(LA, LB))   # rows: RC of A, columns: RC of B
    if (phi[1, 1] + phi[2, 2] >= phi[1, 2] + phi[2, 1]) c(phi[1, 1], phi[2, 2])
    else c(phi[1, 2], phi[2, 1])
  }
  
  set.seed(seed)
  n <- nrow(dat)
  res <- matrix(NA_real_, B, nfac, dimnames = list(NULL, paste0("RC", 1:nfac)))
  for (b in seq_len(B)) {
    idx <- sample.int(n, floor(n / 2))
    res[b, ] <- match_pair(safe_load(dat[idx, , drop = FALSE]),
                           safe_load(dat[-idx, , drop = FALSE]))
  }
  
  summ <- t(apply(res, 2, function(x)
    c(median = median(x, na.rm = TRUE),
      q025   = unname(quantile(x, 0.025, na.rm = TRUE)),
      share_ge_.85 = mean(x >= 0.85, na.rm = TRUE),
      failed = sum(is.na(x)))))
  
  list(congruence = res, summary = round(summ, 2))
}

## Usage:
##   source("split_half_congruence.R")
##   sh <- lapply(groups[c("FIELD", "CORE", "DEM")],
##                function(cols) split_half_congruence(my.vars, cols))
##   lapply(sh, `[[`, "summary")