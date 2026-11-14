## regression_CORE.R
## Regression of the two CORE components (recognition of the crisis, severity in
## one's own field) on the components of the other item groups.
## Built on top of PCA_analysis.R, so numbering and signs of components match the report.
##
## Run after:  source("data_and_R_src/report_numbers.R")
## Needs in the session: my.vars, groups, res  (res <- PCA.analysis(my.vars, groups)).
## Optional: my (original labelled data), used only for the data check in step 1.
##
## Steps
##   1. data check: are codes 1-2 of vek, kar_titul, obor_tlak absent or lost as NA?
##   2. score weights fixed from the oblimin solutions (same numbering/signs as the report)
##   3. block-wise models on the original data, one predictor group at a time
##   4. joint models: multiple imputation of the items, scores recomputed in every
##      imputed data set, Rubin's rules; complete-case model for comparison
##   5. item-level models (logistic / proportional odds) as a robustness check
##
## Output: printed tables and CSV files in the folder "regression_CORE".

suppressPackageStartupMessages({
  library(psych); library(mice); library(sandwich); library(MASS)
})

## ---- 0. specification ------------------------------------------------------
if (!exists("m_imp")) m_imp <- 20   # number of imputed data sets (set m_imp <- 5 earlier for a quick test)
mid     <- TRUE               # drop rows whose outcome was originally missing
                              # (multiple imputation, then deletion; von Hippel 2007)
out_dir <- "regression_CORE"
dir.create(out_dir, showWarnings = FALSE)

if (!exists("res")) {
  source("data_and_R_src/PCA_analysis.R")
  res <- PCA.analysis(my.vars, groups)
}
stopifnot(all(c("FACTORS", "FIELD", "FREQ", "DEM", "NOREP", "CORE") %in% names(res)))

# predictor blocks: component scores (FAC1 ... NOR1) and single variables
blocks <- list(
  FACTORS = c("FAC1", "FAC2"),                  # oblimin RC1, RC2
  FIELD   = c("FLD1", "FLD2"),                  # oblimin RC1, RC2
  FREQ    = c("FRQ1", "frek_positiv"),          # oblimin RC1 + the single item
  DEM     = c("DEM1", "gender", "top_journal"), # seniority + two single variables
  NOREP   = c("NOR1", "nerep_jaCizi"))          # first component + one single item
# outcomes: oblimin components of CORE
outcomes   <- c("COR1", "COR2")
core_items <- colnames(res$CORE$R)

## ---- helpers -----------------------------------------------------------------
to_num <- function(col) {                       # same rule as in PCA.analysis()
  if (is.null(col)) return(NULL)
  if (is.factor(col)) as.integer(col) else suppressWarnings(as.numeric(col))
}
get_items <- function(g, data) {
  ctr <- res[[g]]$pca$center
  as.data.frame(lapply(data[, names(ctr), drop = FALSE], to_num))
}

## ---- 1. data check -------------------------------------------------------------
cat("\n=== 1. Data check: are codes 1-2 of vek, kar_titul, obor_tlak absent or lost? ===\n")
for (v in c("vek", "kar_titul", "obor_tlak")) {
  code <- to_num(my.vars[[v]])
  cat("\n", v, ": ", sum(!is.na(code)), " valid, ", sum(is.na(code)), " NA", sep = "")
  if (exists("my") && v %in% names(my)) {
    cat("; answered in the original data but NA after coding: ",
        sum(!is.na(my[[v]]) & is.na(code)), "\n", sep = "")
    print(table(label = my[[v]], code = code, useNA = "ifany"))
  } else cat("\n")
}
cat("\nYears in research (kar_praxe) by missingness of age (vek):\n")
print(table(age_missing = is.na(to_num(my.vars$vek)), years = to_num(my.vars$kar_praxe)))
cat("If age is missing mostly for 'less than 5 years', the missing ages are young\n",
    "respondents lost in coding. Fix this before imputing (MAR would not hold).\n", sep = "")

## ---- 2. fixed score weights ----------------------------------------------------
# oblimin solution of the group correlation matrix, columns matched to the varimax
# numbering and signed so that the largest |loading| is positive (as in the report)
oblimin_weights <- function(g) {
  R  <- res[[g]]$R
  vm <- unclass(principal(R, nfactors = 2, rotate = "varimax")$loadings)[, 1:2]
  ob <- principal(R, nfactors = 2, rotate = "oblimin")
  P  <- unclass(ob$loadings)[, 1:2]
  Wt <- ob$weights[, 1:2, drop = FALSE]
  rownames(Wt) <- rownames(P) <- colnames(R)      # weights rows follow the columns of R
  phi <- ob$Phi
  cg <- abs(factor.congruence(P, vm))
  if (cg[1, 2] + cg[2, 1] > cg[1, 1] + cg[2, 2]) {
    P <- P[, 2:1]; Wt <- Wt[, 2:1]; phi <- phi[2:1, 2:1]
  }
  for (j in 1:2) if (P[which.max(abs(P[, j])), j] < 0) {
    P[, j] <- -P[, j]; Wt[, j] <- -Wt[, j]
    phi[j, setdiff(1:2, j)] <- -phi[j, setdiff(1:2, j)]
    phi[setdiff(1:2, j), j] <- -phi[setdiff(1:2, j), j]
  }
  colnames(P) <- colnames(Wt) <- c("RC1", "RC2")
  list(P = P, W = Wt, phi = phi[1, 2])
}
ob <- lapply(setNames(nm = c("FACTORS", "FIELD", "DEM", "FREQ", "CORE")), oblimin_weights)

cat("\n=== 2. Oblimin pattern loadings used for the scores (compare with the report appendix) ===\n")
for (g in names(ob)) {
  cat("\n", g, "  (component correlation ", round(ob[[g]]$phi, 2), ")\n", sep = "")
  print(round(ob[[g]]$P, 2))
}

# NOREP: first unrotated component (general experience of failed replications)
nor <- principal(res$NOREP$R, nfactors = 1, rotate = "none")
nor_w <- nor$weights[, 1] * ifelse(nor$loadings[which.max(abs(nor$loadings)), 1] < 0, -1, 1)
names(nor_w) <- colnames(res$NOREP$R)

W <- list(FAC1 = ob$FACTORS$W[, 1], FAC2 = ob$FACTORS$W[, 2],
          FLD1 = ob$FIELD$W[, 1],   FLD2 = ob$FIELD$W[, 2],
          DEM1 = ob$DEM$W[, 1],     FRQ1 = ob$FREQ$W[, 1],
          NOR1 = nor_w,
          COR1 = ob$CORE$W[, 1],    COR2 = ob$CORE$W[, 2])
stopifnot(!anyNA(unlist(W)), all(vapply(W, function(w) !is.null(names(w)), logical(1))))
grp_of <- c(FAC1 = "FACTORS", FAC2 = "FACTORS", FLD1 = "FIELD", FLD2 = "FIELD",
            DEM1 = "DEM", FRQ1 = "FREQ", NOR1 = "NOREP", COR1 = "CORE", COR2 = "CORE")

# score = standardised items %*% weights; NA if an item of the group is missing
make_scores <- function(data, which = names(W)) {
  sapply(which, function(nm) {
    g <- grp_of[[nm]]; ctr <- res[[g]]$pca$center
    X <- as.matrix(get_items(g, data))[, names(ctr), drop = FALSE]
    drop(scale(X, center = ctr, scale = res[[g]]$pca$scale) %*% W[[nm]][names(ctr)])
  })
}
S0  <- make_scores(my.vars)
ref <- list(mean = colMeans(S0, na.rm = TRUE), sd = apply(S0, 2, stats::sd, na.rm = TRUE))
std <- function(S) sweep(sweep(S, 2, ref$mean[colnames(S)]), 2, ref$sd[colnames(S)], "/")
fp  <- to_num(my.vars$frek_positiv); fp_ref <- c(mean(fp, na.rm = TRUE), sd(fp, na.rm = TRUE))

# single variables; gender 0/1 (1 = woman), top_journal 0/1, nerep_jaCizi 0/1,
# frek_positiv standardised. Binary variables are not standardised.
singles <- function(data, which = c("gender", "top_journal", "frek_positiv", "nerep_jaCizi")) {
  s <- list(gender       = to_num(data$gender),
            top_journal  = to_num(data$kar_nejlepsi) - 1,
            frek_positiv = (to_num(data$frek_positiv) - fp_ref[1]) / fp_ref[2],
            nerep_jaCizi = to_num(data$nerep_jaCizi))
  as.data.frame(s[which])
}
D0 <- cbind(as.data.frame(std(S0)), singles(my.vars))     # original data, no imputation
D0[core_items] <- lapply(my.vars[core_items], to_num)
D0$krize_rep_noDK <- ifelse(D0$krize_rep == 2, NA, D0$krize_rep)  # sensitivity: "don't know" removed

cat("\n=== Scores on the original data (complete cases per group) ===\n")
print(colSums(!is.na(D0[, c(names(W), "gender", "top_journal", "frek_positiv", "nerep_jaCizi")])))
cat("Correlation of the two CORE scores: ",
    round(cor(D0$COR1, D0$COR2, use = "complete.obs"), 2), "\n", sep = "")

## ---- fitting and pooling functions ---------------------------------------------
fit_lm <- function(y, x, d) {                   # OLS with HC3 standard errors
  m <- lm(reformulate(x, y), data = d)
  list(est = coef(m), vcov = vcovHC(m, type = "HC3"), n = nobs(m),
       k = length(coef(m)), r2 = summary(m)$r.squared)
}
fit_item <- function(item, x, d) {              # logistic (binary) or proportional odds
  d2 <- d[, c(item, x), drop = FALSE]
  names(d2)[1] <- "y"
  n <- sum(complete.cases(d2))
  if (item == "rep_type") {
    m <- glm(reformulate(x, "y"), family = binomial, data = d2)
    est <- coef(m); V <- vcov(m)
  } else {
    d2$y <- factor(d2$y, ordered = TRUE)
    m <- polr(reformulate(x, "y"), data = d2, Hess = TRUE)
    est <- coef(m); V <- vcov(m)[names(est), names(est)]
  }
  list(est = est, vcov = V, n = n, k = length(est) + 1, r2 = NA_real_)
}
tidy_fit <- function(f) {                       # one fit, no imputation
  se <- sqrt(diag(f$vcov)[names(f$est)]); df <- f$n - f$k; t <- f$est / se
  if (df < 1) df <- NA_real_
  data.frame(term = names(f$est), estimate = unname(f$est), se = unname(se), df = df,
             t = unname(t), p = 2 * pt(-abs(t), df),
             lo = unname(f$est - qt(.975, df) * se), hi = unname(f$est + qt(.975, df) * se),
             fmi = NA_real_, n = f$n, r2 = f$r2, row.names = NULL)
}
pool_fits <- function(fits) {                   # Rubin's rules, Barnard-Rubin df
  m <- length(fits); terms <- names(fits[[1]]$est)
  Q <- t(sapply(fits, function(f) f$est[terms]))
  U <- t(sapply(fits, function(f) diag(f$vcov)[terms]))
  Qbar <- colMeans(Q); Ubar <- colMeans(U); B <- apply(Q, 2, stats::var)
  Tt <- Ubar + (1 + 1 / m) * B
  lambda <- (1 + 1 / m) * B / Tt
  dfc <- fits[[1]]$n - fits[[1]]$k
  df_old <- ifelse(lambda > 0, (m - 1) / lambda^2, Inf)
  df_obs <- (dfc + 1) / (dfc + 3) * dfc * (1 - lambda)
  df <- 1 / (1 / df_old + 1 / df_obs)
  se <- sqrt(Tt); t <- Qbar / se
  data.frame(term = terms, estimate = unname(Qbar), se = unname(se), df = unname(df),
             t = unname(t), p = unname(2 * pt(-abs(t), df)),
             lo = unname(Qbar - qt(.975, df) * se), hi = unname(Qbar + qt(.975, df) * se),
             fmi = unname(lambda), n = fits[[1]]$n, r2 = NA_real_, row.names = NULL)
}
show <- function(df) {
  op <- options(width = 200); on.exit(options(op))
  df <- df[df$term != "(Intercept)", setdiff(names(df), c("df", "t"))]
  df <- df[, colSums(!is.na(df)) > 0, drop = FALSE]            # hide all-NA columns
  num <- vapply(df, base::is.numeric, logical(1)); df[num] <- lapply(df[num], base::round, 3)
  if ("p" %in% names(df)) df$p <- round(df$p, 4)
  print(df, row.names = FALSE)
}

## ---- 3. block-wise models (original data, one group at a time) -----------------
cat("\n=== 3. Block-wise linear models (OLS, HC3); outcome and predictors in SD units ===\n")
dem_terms <- blocks$DEM
block_res <- list()
for (y in outcomes) for (b in names(blocks)) for (adj in c(FALSE, TRUE)) {
  if (adj && b == "DEM") next
  x <- c(blocks[[b]], if (adj) dem_terms)
  t <- tidy_fit(fit_lm(y, x, D0))
  t <- t[t$term %in% blocks[[b]], ]
  block_res[[length(block_res) + 1]] <-
    cbind(outcome = y, block = b, adjusted_for_DEM = adj, t)
}
block_res <- do.call(base::rbind, block_res)
show(block_res)
write.csv(block_res, file.path(out_dir, "block_models.csv"), row.names = FALSE)

## ---- 4. multiple imputation ------------------------------------------------------
# Items of FACTORS, FIELD, FREQ, DEM and CORE are imputed. NOREP is not: answers
# "I do not do experiments" / "I do not remember" are structurally missing, so all
# models with NOREP use only respondents with observed NOREP items.
imp_groups <- c("FACTORS", "FIELD", "FREQ", "DEM", "CORE")
imp_vars   <- unique(unlist(lapply(imp_groups, function(g) names(res[[g]]$pca$center))))
imp_dat    <- as.data.frame(lapply(my.vars[, imp_vars], to_num))
cat("\n=== 4. Multiple imputation: ", m_imp, " data sets, predictive mean matching, ",
    length(imp_vars), " items ===\n", sep = "")
imp <- mice(imp_dat, m = m_imp, method = "pmm", maxit = 10, seed = 1, printFlag = FALSE)
if (!is.null(imp$loggedEvents)) {
  cat("mice logged events (variables dropped as predictors, usually collinearity):\n")
  print(table(imp$loggedEvents$meth, imp$loggedEvents$out))
}
no_nor <- setdiff(names(W), "NOR1")
Dlist <- lapply(seq_len(m_imp), function(i) {
  L <- complete(imp, i)
  d <- cbind(as.data.frame(std(make_scores(L, no_nor))),
             singles(L, c("gender", "top_journal", "frek_positiv")),
             D0[, c("NOR1", "nerep_jaCizi")])             # NOREP stays as observed
  d[core_items] <- L[core_items]
  d$krize_rep_noDK <- ifelse(d$krize_rep == 2, NA, d$krize_rep)
  d
})
keep_for <- function(y) if (mid) !is.na(D0[[y]]) else rep(TRUE, nrow(D0))

x_all <- function(use_norep) unlist(blocks[c("FACTORS", "FIELD", "FREQ", "DEM",
                                             if (use_norep) "NOREP")], use.names = FALSE)

run_joint <- function(use_norep) {
  x <- x_all(use_norep); lab <- if (use_norep) "with NOREP" else "without NOREP"
  lin <- do.call(base::rbind, lapply(outcomes, function(y) {
    k <- keep_for(y)
    fits <- lapply(Dlist, function(d) fit_lm(y, x, d[k, ]))
    cbind(outcome = y, model = paste("MI,", lab), pool_fits(fits))
  }))
  cc <- do.call(base::rbind, lapply(outcomes, function(y)
    cbind(outcome = y, model = paste("complete case,", lab), tidy_fit(fit_lm(y, x, D0)))))
  list(lin = lin, cc = cc)
}

cat("\n--- Joint linear models, MI (without NOREP: all respondents) ---\n")
J1 <- run_joint(FALSE); show(J1$lin)
cat("\n--- Joint linear models, MI (with NOREP: respondents with observed NOREP items) ---\n")
J2 <- run_joint(TRUE);  show(J2$lin)
cat("\n--- Complete-case joint models, for comparison ---\n")
show(rbind(J1$cc, J2$cc))
joint_res <- rbind(J1$lin, J2$lin, J1$cc, J2$cc)
write.csv(joint_res, file.path(out_dir, "joint_models.csv"), row.names = FALSE)

cat("\n--- Variance inflation factors (first imputed data set) ---\n")
vif <- function(d, x) { X <- d[, x]; X <- X[complete.cases(X), ]
                        data.frame(term = x, VIF = round(diag(solve(cor(X))), 2)) }
print(vif(Dlist[[1]], x_all(TRUE)), row.names = FALSE)

## ---- 5. item-level models --------------------------------------------------------
cat("\n=== 5. Item-level models (MI): logistic for rep_type, proportional odds otherwise ===\n")
cat("Estimates are log-odds; OR = exp(estimate). Positive = higher/more serious category.\n")
items_out <- c(core_items, "krize_rep_noDK")
item_res <- list()
for (use_norep in c(FALSE, TRUE)) {
  x <- x_all(use_norep)
  for (it in items_out) {
    orig <- D0[[it]]; k <- if (mid) !is.na(orig) else rep(TRUE, nrow(D0))
    fits <- tryCatch(lapply(Dlist, function(d) fit_item(it, x, d[k, ])),
                     error = function(e) { message(it, ": ", conditionMessage(e)); NULL })
    if (is.null(fits)) next
    p <- pool_fits(fits)
    p$OR <- exp(p$estimate); p$OR_lo <- exp(p$lo); p$OR_hi <- exp(p$hi)
    item_res[[length(item_res) + 1]] <-
      cbind(outcome = it, model = if (use_norep) "with NOREP" else "without NOREP", p)
  }
}
item_res <- do.call(base::rbind, item_res)
show(item_res[, c("outcome", "model", "term", "estimate", "se", "p", "OR", "OR_lo", "OR_hi", "n")])
write.csv(item_res, file.path(out_dir, "item_models.csv"), row.names = FALSE)

# proportional-odds check on the first imputed data set (ordinal::clm, nominal_test)
if (requireNamespace("ordinal", quietly = TRUE)) {
  cat("\n--- Proportional-odds check (first imputed data set, small p = assumption doubtful) ---\n")
  for (it in setdiff(items_out, "rep_type")) {
    d <- Dlist[[1]][!is.na(D0[[it]]), c(it, x_all(FALSE))]
    names(d)[1] <- "y"; d$y <- factor(d$y, ordered = TRUE)
    cat("\n", it, "\n", sep = "")
    tryCatch(print(ordinal::nominal_test(ordinal::clm(reformulate(x_all(FALSE), "y"), data = d))),
             error = function(e) cat("  not available: ", conditionMessage(e), "\n", sep = ""))
  }
}

saveRDS(list(blocks = block_res, joint = joint_res, items = item_res,
             weights = W, ref = ref), file.path(out_dir, "regression_CORE.rds"))
cat("\nCSV files and regression_CORE.rds written to '", out_dir, "'.\n", sep = "")
