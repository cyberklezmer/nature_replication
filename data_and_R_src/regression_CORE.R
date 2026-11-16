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
##   1. data check: the corrected coding of vek, kar_titul, obor_tlak is in my.vars
##   2. score weights fixed from the oblimin solutions (same numbering/signs as the report);
##      optional comparison with C_oblique.csv
##   3. block-wise models on the original data, one predictor group at a time
##   4. joint models: multiple imputation of the items, scores recomputed in every
##      imputed data set, Rubin's rules; complete-case model for comparison;
##      sensitivity with a one-component FIELD
##   4b. Wald test of FLD1 = FLD2 per outcome (block-wise, MI, complete case); full covariance
##      matrix of the estimates, MI pooled by Rubin's rules -> wald_FLD.csv, wald_FLD_rows.tex
##   4c. alternative component solutions: (a) 3 and 4 oblimin components per predictor group
##      (outcomes COR1, COR2); (b) one unrotated component per group (outcomes COR0, COR1, COR2);
##      R^2 comparison with the 2-component model -> alt_models.csv, alt_fit_summary.csv, alt_loadings.csv
##   5. item-level models (logistic / proportional odds) as a robustness check
##
## Outcomes: COR1 (recognition of the crisis), COR2 (severity in one's own field) and,
## as a sensitivity analysis, COR0 = first unrotated component of CORE. Parallel analysis
## retained one component for CORE (2nd eigenvalue 0.87 vs 1.05 random), so the one- and the
## two-component description are both reported.
##
## Output: printed tables and CSV files in the folder "regression_CORE".

suppressPackageStartupMessages({
  library(psych); library(mice); library(sandwich); library(MASS)
  library(caret)   # for 10-fold cross-validation in step 4c
})

## ---- 0. specification ------------------------------------------------------
if (!exists("m_imp")) m_imp <- 20   # number of imputed data sets (set m_imp <- 5 earlier for a quick test)
mid     <- TRUE               # drop rows whose outcome was originally missing
                              # (multiple imputation, then deletion; von Hippel 2007)
out_dir <- "regression_CORE"
oblique_csv <- "C_oblique.csv"   # optional: oblique pattern loadings of the report (adjust path)
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
# extra block for a sensitivity analysis: FIELD as ONE component (first unrotated PC);
# parallel analysis gave one component for FIELD and its second component was the least
# stable one in the split-half check
blocks_extra <- list(FIELD_1comp = "FLD0")
# outcomes: oblimin components of CORE (+ one-component description as sensitivity)
outcomes   <- c("COR1", "COR2")
outcome_1c <- "COR0"
outcomes_all <- c(outcomes, outcome_1c)
core_items <- colnames(res$CORE$R)

## ---- helpers -----------------------------------------------------------------
to_num <- function(col) {                       # same rule as in PCA.analysis()
  if (is.null(col)) return(NULL)
  if (is.factor(col)) as.integer(col) else suppressWarnings(as.numeric(col))
}
bin01 <- function(x) {                          # binary variable as 0/1, whatever its coding
  if (is.null(x)) return(NULL)
  x <- to_num(x); r <- range(x, na.rm = TRUE)
  stopifnot(r[2] - r[1] == 1)
  x - r[1]
}
get_items <- function(g, data) {
  ctr <- res[[g]]$pca$center
  as.data.frame(lapply(data[, names(ctr), drop = FALSE], to_num))
}

## ---- 1. data check -------------------------------------------------------------
# After the corrected coding all codes of vek (1-5), kar_titul (1-5) and obor_tlak (1-4)
# are present and none of these variables has missing values.
cat("\n=== 1. Data check: corrected coding of vek, kar_titul, obor_tlak ===\n")
expected_codes <- list(vek = 1:5, kar_titul = 1:5, obor_tlak = 1:4)
for (v in names(expected_codes)) {
  code <- to_num(my.vars[[v]])
  cat("\n", v, ": ", sum(!is.na(code)), " valid, ", sum(is.na(code)), " NA\n", sep = "")
  print(table(code, useNA = "ifany"))
  if (anyNA(code) || !setequal(unique(code), expected_codes[[v]]))
    warning(v, ": codes differ from the expected ", min(expected_codes[[v]]), "-",
            max(expected_codes[[v]]), " without NA. Is my.vars the data with the corrected coding?")
}
expected_cc <- c(FACTORS = 343, NOREP = 451, FREQ = 582, FIELD = 878, DEM = 873, CORE = 700)
cc_now <- sapply(names(expected_cc), function(g) sum(complete.cases(get_items(g, my.vars))))
cat("\nComplete cases per group (now / in the report):\n")
print(rbind(now = cc_now, report = expected_cc))
if (any(cc_now != expected_cc)) warning("Complete cases differ from the report: check the data/coding.")

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

# comparison with the oblique pattern loadings of the report (C_oblique.csv). Only the sign of
# a component may differ (in the csv the second DEM component has gender as -0.83).
if (file.exists(oblique_csv)) {
  cat("\nMax. |difference| to the pattern loadings in ", oblique_csv, " (up to the sign of a component):\n", sep = "")
  cc <- read.csv(oblique_csv)
  for (g in names(ob)) for (j in 1:2) {
    cr <- cc[cc$group == g & cc$component == paste0("PC", j), ]
    p  <- ob[[g]]$P[cr$variable, j]
    d  <- min(max(abs(p - cr$pattern)), max(abs(p + cr$pattern)))
    cat(sprintf("  %-8s RC%d: %.4f%s\n", g, j, d, if (d > 0.01) "   <-- differs, check" else ""))
  }
} else cat("\n(", oblique_csv, " not found: comparison with the report skipped)\n", sep = "")

# first unrotated component of a group; largest |loading| positive.
# NOREP: general experience of failed replications; FIELD, CORE: one-component descriptions
pc1_weights <- function(g) {
  p <- principal(res[[g]]$R, nfactors = 1, rotate = "none")
  w <- p$weights[, 1] * ifelse(p$loadings[which.max(abs(p$loadings[, 1])), 1] < 0, -1, 1)
  names(w) <- colnames(res[[g]]$R)
  w
}
nor_w <- pc1_weights("NOREP")

W <- list(FAC1 = ob$FACTORS$W[, 1], FAC2 = ob$FACTORS$W[, 2],
          FLD1 = ob$FIELD$W[, 1],   FLD2 = ob$FIELD$W[, 2],
          DEM1 = ob$DEM$W[, 1],     FRQ1 = ob$FREQ$W[, 1],
          NOR1 = nor_w,
          COR1 = ob$CORE$W[, 1],    COR2 = ob$CORE$W[, 2],
          FLD0 = pc1_weights("FIELD"), COR0 = pc1_weights("CORE"))
stopifnot(!anyNA(unlist(W)), all(vapply(W, function(w) !is.null(names(w)), logical(1))))
grp_of <- c(FAC1 = "FACTORS", FAC2 = "FACTORS", FLD1 = "FIELD", FLD2 = "FIELD",
            DEM1 = "DEM", FRQ1 = "FREQ", NOR1 = "NOREP", COR1 = "CORE", COR2 = "CORE",
            FLD0 = "FIELD", COR0 = "CORE")

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
# frek_positiv standardised. Binary variables are not standardised. Only the requested
# variables are built (NOREP items are not in the imputed data).
singles <- function(data, which = c("gender", "top_journal", "frek_positiv", "nerep_jaCizi")) {
  f <- list(gender       = function() bin01(data$gender),
            top_journal  = function() bin01(data$kar_nejlepsi),
            frek_positiv = function() (to_num(data$frek_positiv) - fp_ref[1]) / fp_ref[2],
            nerep_jaCizi = function() bin01(data$nerep_jaCizi))
  as.data.frame(lapply(f[which], function(g) g()))
}
D0 <- cbind(as.data.frame(std(S0)), singles(my.vars))     # original data, no imputation
D0[core_items] <- lapply(my.vars[core_items], to_num)
D0$krize_rep_noDK <- ifelse(D0$krize_rep == 2, NA, D0$krize_rep)  # sensitivity: "don't know" removed

cat("\n=== Scores on the original data (complete cases per group) ===\n")
print(colSums(!is.na(D0[, c(names(W), "gender", "top_journal", "frek_positiv", "nerep_jaCizi")])))
cat("Correlation of the two CORE scores: ",
    round(cor(D0$COR1, D0$COR2, use = "complete.obs"), 2), "\n", sep = "")
cat("Correlation of the one-component CORE score with COR1 / COR2: ",
    round(cor(D0$COR0, D0$COR1, use = "complete.obs"), 2), " / ",
    round(cor(D0$COR0, D0$COR2, use = "complete.obs"), 2), "\n", sep = "")
cat("\nCoding of the binary predictors (0/1 as used in the models; 1 = woman, 1 = top journal):\n")
print(table(gender = D0$gender, useNA = "ifany")); print(table(top_journal = D0$top_journal, useNA = "ifany"))
cat("Gender 'other' (", sum(is.na(D0$gender)), " respondents) is set to missing and is not imputed.\n", sep = "")

## ---- fitting and pooling functions ---------------------------------------------
fit_lm <- function(y, x, d) {                   # OLS with HC3 standard errors
  m <- lm(reformulate(x, y), data = d)
  # AIC and BIC use the Gaussian log-likelihood at OLS residuals (same as stats::AIC/BIC).
  # HC3 corrects the standard errors but NOT the likelihood, so AIC/BIC remain valid
  # for comparing models on the same outcome and the same rows.
  list(est = coef(m), vcov = vcovHC(m, type = "HC3"), n = nobs(m),
       k = length(coef(m)), r2 = summary(m)$r.squared,
       aic = AIC(m), bic = BIC(m))
}
fit_item <- function(item, x, d) {              # logistic (binary) or proportional odds
  d2 <- d[, c(item, x), drop = FALSE]
  names(d2)[1] <- "y"
  n <- sum(complete.cases(d2))
  if (item == "rep_type") {
    d2$y <- d2$y - min(d2$y, na.rm = TRUE)         # 0/1 whatever the original coding
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
allb <- c(blocks, blocks_extra)
block_res <- list()
for (y in outcomes_all) for (b in names(allb)) for (adj in c(FALSE, TRUE)) {
  if (adj && b == "DEM") next
  x <- c(allb[[b]], if (adj) dem_terms)
  t <- tidy_fit(fit_lm(y, x, D0))
  t <- t[t$term %in% allb[[b]], ]
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
  d$gender[is.na(D0$gender)] <- NA                        # gender 'other': not imputed, excluded
  d[core_items] <- L[core_items]
  d$krize_rep_noDK <- ifelse(d$krize_rep == 2, NA, d$krize_rep)
  d
})
keep_for <- function(y) if (mid) !is.na(D0[[y]]) else rep(TRUE, nrow(D0))

x_all <- function(use_norep, field1 = FALSE) {
  b <- blocks
  if (field1) b$FIELD <- blocks_extra$FIELD_1comp        # FIELD as one component
  unlist(b[c("FACTORS", "FIELD", "FREQ", "DEM", if (use_norep) "NOREP")], use.names = FALSE)
}

run_joint <- function(use_norep, field1 = FALSE) {
  x <- x_all(use_norep, field1)
  lab <- paste0(if (use_norep) "with NOREP" else "without NOREP", if (field1) ", FIELD 1 comp." else "")
  lin <- do.call(base::rbind, lapply(outcomes_all, function(y) {
    k <- keep_for(y)
    fits <- lapply(Dlist, function(d) fit_lm(y, x, d[k, ]))
    cbind(outcome = y, model = paste("MI,", lab), pool_fits(fits))
  }))
  cc <- do.call(base::rbind, lapply(outcomes_all, function(y)
    cbind(outcome = y, model = paste("complete case,", lab), tidy_fit(fit_lm(y, x, D0)))))
  list(lin = lin, cc = cc)
}

cat("\n--- Joint linear models, MI (without NOREP: all respondents) ---\n")
J1 <- run_joint(FALSE); show(J1$lin)
cat("\n--- Joint linear models, MI (with NOREP: respondents with observed NOREP items) ---\n")
J2 <- run_joint(TRUE);  show(J2$lin)
cat("\n--- Sensitivity: FIELD as one component (MI, without NOREP) ---\n")
J3 <- run_joint(FALSE, field1 = TRUE); show(J3$lin)
cat("\n--- Complete-case joint models, for comparison ---\n")
show(rbind(J1$cc, J2$cc))
joint_res <- rbind(J1$lin, J2$lin, J3$lin, J1$cc, J2$cc)
write.csv(joint_res, file.path(out_dir, "joint_models.csv"), row.names = FALSE)

cat("\n--- Variance inflation factors (first imputed data set) ---\n")
vif <- function(d, x) { X <- d[, x]; X <- X[complete.cases(X), ]
                        data.frame(term = x, VIF = round(diag(solve(cor(X))), 2)) }
print(vif(Dlist[[1]], x_all(TRUE)), row.names = FALSE)

## ---- 4b. Wald test of FLD1 = FLD2 (per outcome) ---------------------------------
# H0: the effects of the two FIELD components on the outcome are equal (FLD1 - FLD2 = 0).
# Both scores are in SD units, so the coefficients are comparable. The test uses the full
# covariance matrix of the estimates (HC3 in every fit), i.e. it includes Cov(b_FLD1, b_FLD2),
# which is NOT available from the CSV tables (they hold only standard errors).
# One fit:  t = (b1 - b2) / sqrt(Var(b1) + Var(b2) - 2 Cov(b1, b2)),  Wald F = t^2 on 1 and n-k df.
# MI:       the contrast is computed in every imputed data set and pooled by Rubin's rules
#           (Barnard-Rubin df), the same way as the coefficients; F = t^2 on 1 and df.
contrast_fit <- function(f, a = "FLD1", b = "FLD2") {
  cv <- setNames(rep(0, length(f$est)), names(f$est)); cv[a] <- 1; cv[b] <- -1
  V  <- f$vcov[names(cv), names(cv)]
  list(b1 = unname(f$est[a]), b2 = unname(f$est[b]),
       diff = sum(cv * f$est), var = drop(crossprod(cv, V %*% cv)))
}
wald_one <- function(f) {
  k <- contrast_fit(f); df <- f$n - f$k; tt <- k$diff / sqrt(k$var)
  data.frame(FLD1 = k$b1, FLD2 = k$b2, diff = k$diff, se = sqrt(k$var), df = df,
             F = tt^2, p = 2 * pt(-abs(tt), df), fmi = NA_real_, n = f$n)
}
wald_pool <- function(fits) {
  m <- length(fits); ks <- lapply(fits, contrast_fit)
  Q <- vapply(ks, function(k) k$diff, numeric(1)); U <- vapply(ks, function(k) k$var, numeric(1))
  Qbar <- mean(Q); Ubar <- mean(U); B <- stats::var(Q)
  Tt <- Ubar + (1 + 1 / m) * B; lambda <- (1 + 1 / m) * B / Tt
  dfc <- fits[[1]]$n - fits[[1]]$k
  df_old <- if (lambda > 0) (m - 1) / lambda^2 else Inf
  df_obs <- (dfc + 1) / (dfc + 3) * dfc * (1 - lambda)
  df <- 1 / (1 / df_old + 1 / df_obs); tt <- Qbar / sqrt(Tt)
  data.frame(FLD1 = mean(vapply(ks, function(k) k$b1, numeric(1))),
             FLD2 = mean(vapply(ks, function(k) k$b2, numeric(1))),
             diff = Qbar, se = sqrt(Tt), df = df, F = tt^2,
             p = 2 * pt(-abs(tt), df), fmi = lambda, n = fits[[1]]$n)
}

cat("\n=== 4b. Wald test of FLD1 = FLD2 per outcome (F on 1 and df d.f.; diff = FLD1 - FLD2) ===\n")
wald_res <- list()
add_wald <- function(analysis, y, w) wald_res[[length(wald_res) + 1]] <<-
  cbind(analysis = analysis, outcome = y, w)
for (y in outcomes_all) {                                   # block-wise, original data
  for (adj in c(FALSE, TRUE)) {
    x <- c(blocks$FIELD, if (adj) dem_terms)
    add_wald(paste0("block FIELD", if (adj) ", adjusted for DEM" else ""), y,
             wald_one(fit_lm(y, x, D0)))
  }
}
for (use_norep in c(FALSE, TRUE)) {                         # joint models
  x <- x_all(use_norep); lab <- if (use_norep) "with NOREP" else "without NOREP"
  for (y in outcomes_all) {
    k <- keep_for(y)
    add_wald(paste("MI,", lab), y, wald_pool(lapply(Dlist, function(d) fit_lm(y, x, d[k, ]))))
    add_wald(paste("complete case,", lab), y, wald_one(fit_lm(y, x, D0)))
  }
}
wald_res <- do.call(base::rbind, wald_res)
op <- options(width = 200)
print(within(wald_res, { FLD1 <- round(FLD1, 3); FLD2 <- round(FLD2, 3); diff <- round(diff, 3)
                         se <- round(se, 3); df <- round(df, 1); F <- round(F, 2)
                         p <- round(p, 4); fmi <- round(fmi, 3) }), row.names = FALSE)
options(op)
write.csv(wald_res, file.path(out_dir, "wald_FLD.csv"), row.names = FALSE)

# rows for the LaTeX table of the report (MI without NOREP, MI with NOREP, complete case without NOREP)
fm <- function(x, d = 2) { s <- formatC(x, format = "f", digits = d); s[as.numeric(s) == 0] <- formatC(0, format = "f", digits = d)
                           sub("^-", "$-$", s) }
fmt_p <- function(p) ifelse(p < 0.001, "$<0.001$", formatC(p, format = "f", digits = 3))
tab <- wald_res[wald_res$analysis %in% c("MI, without NOREP", "MI, with NOREP",
                                          "complete case, without NOREP"), ]
tab <- tab[order(match(tab$analysis, c("MI, without NOREP", "MI, with NOREP", "complete case, without NOREP")),
                 match(tab$outcome, outcomes_all)), ]
lab <- c("MI, without NOREP" = "MI, all respondents", "MI, with NOREP" = "MI, observed NOREP",
         "complete case, without NOREP" = "Complete cases")
rows <- sprintf("%s & %s & %d & %s & %s & %s & %.3f & %.2f & %s \\\\", lab[tab$analysis], tab$outcome, tab$n,
                fm(tab$FLD1), fm(tab$FLD2), fm(tab$diff), tab$se, tab$F, fmt_p(tab$p))
writeLines(rows, file.path(out_dir, "wald_FLD_rows.tex"))

## ---- 4c. Alternative component solutions ------------------------------------------
# (a) 3 and 4 oblimin components per predictor group (FACTORS, FIELD, FREQ, DEM) instead of 2;
#     outcomes COR1 and COR2 (two oblimin components of CORE, as in the main analysis).
#     The components replace the 2-component scores AND the single items of these groups
#     (frek_positiv, gender, top_journal are items of FREQ / DEM and are part of the components).
#     NOREP has only 3 items and keeps NOR1 + nerep_jaCizi.
# (a2) FACTORS-only 3-component model: 3 oblimin components for FACTORS (FAC3c1--FAC3c3),
#     all other groups exactly as in the benchmark (FLD1, FLD2, FRQ1, frek_positiv,
#     DEM1, gender, top_journal; + NOR1, nerep_jaCizi when use_norep = TRUE).
#     Outcomes COR1 and COR2.
# (a3) FACTORS-only 4-component model: 4 oblimin components for FACTORS (FAC4c1--FAC4c4),
#     all other groups exactly as in the benchmark.
#     Outcomes COR1 and COR2. This is the focused alternative motivated by the
#     substantively interpretable four-way split of FACTORS items.
# (b) ONE unrotated component per predictor group (FAC0, FLD0, FRQ0, DEM0; NOR1 is already the
#     first component of NOREP), no single items; outcomes COR0 (the one-component CORE
#     description) and, for comparison, COR1 and COR2.
# Everything else is as in the main joint models: the same imputed items, scores recomputed in
# every imputed data set (weights fixed from the full-sample solution), rows with originally
# missing outcome dropped, gender "other" excluded, OLS with HC3, Rubin's rules; models with
# NOREP use only respondents with observed NOREP items. The 2-component model of the report is
# refitted with the same code as a benchmark, so that R^2 can be compared on the same rows.
# The numbering of the components follows the order returned by psych::principal() (largest
# variance first); check the printed pattern loadings before interpreting them.
alt_k      <- c(3, 4)
alt_groups <- c(FACTORS = "FAC", FIELD = "FLD", FREQ = "FRQ", DEM = "DEM")

oblimin_k <- function(g, k) {                   # k oblimin components, largest |loading| positive
  R   <- res[[g]]$R
  ob  <- principal(R, nfactors = k, rotate = "oblimin")
  P   <- unclass(ob$loadings)[, 1:k, drop = FALSE]
  Wt  <- ob$weights[, 1:k, drop = FALSE]
  phi <- ob$Phi
  rownames(P) <- rownames(Wt) <- colnames(R)
  sg  <- ifelse(apply(P, 2, function(p) p[which.max(abs(p))]) < 0, -1, 1)
  P   <- sweep(P, 2, sg, "*"); Wt <- sweep(Wt, 2, sg, "*"); phi <- phi * outer(sg, sg)
  cn  <- paste0("RC", seq_len(k))
  colnames(P) <- colnames(Wt) <- cn; dimnames(phi) <- list(cn, cn)
  list(P = P, W = Wt, phi = phi)
}

cat("\n=== 4c. Alternative component solutions ===\n")
W_alt <- list(); grp_alt <- character(0); alt_load <- list()
for (k in alt_k) for (g in names(alt_groups)) {
  o <- oblimin_k(g, k)
  cat("\n", g, ": ", k, " oblimin components, pattern loadings (columns = components)\n", sep = "")
  print(round(o$P, 2))
  cat("Correlations between the components:\n"); print(round(o$phi, 2))
  for (j in seq_len(k)) {
    nm <- paste0(alt_groups[[g]], k, "c", j)
    W_alt[[nm]] <- o$W[, j]; grp_alt[nm] <- g
  }
  alt_load[[length(alt_load) + 1]] <-
    data.frame(group = g, k = k, component = rep(colnames(o$P), each = nrow(o$P)),
               variable = rep(rownames(o$P), times = k), pattern = as.vector(o$P),
               row.names = NULL)
}
one_comp <- c(FACTORS = "FAC0", FREQ = "FRQ0", DEM = "DEM0")      # FIELD: FLD0 and NOREP: NOR1 exist
for (gg in names(one_comp)) {                                       # one unrotated component
  nm <- one_comp[[gg]]
  W_alt[[nm]] <- pc1_weights(gg); grp_alt[nm] <- gg
  pl <- principal(res[[gg]]$R, nfactors = 1, rotate = "none")$loadings[, 1]
  pl <- pl * ifelse(pl[which.max(abs(pl))] < 0, -1, 1)
  cat("\n", gg, ": first unrotated component (", nm, "), loadings\n", sep = ""); print(round(pl, 2))
}
alt_load <- do.call(base::rbind, alt_load)
write.csv(alt_load, file.path(out_dir, "alt_loadings.csv"), row.names = FALSE)

make_scores_alt <- function(data, which = names(W_alt)) {
  sapply(which, function(nm) {
    g <- grp_alt[[nm]]; ctr <- res[[g]]$pca$center
    X <- as.matrix(get_items(g, data))[, names(ctr), drop = FALSE]
    drop(scale(X, center = ctr, scale = res[[g]]$pca$scale) %*% W_alt[[nm]][names(ctr)])
  })
}
S0a  <- make_scores_alt(my.vars)
ref2 <- list(mean = colMeans(S0a, na.rm = TRUE), sd = apply(S0a, 2, stats::sd, na.rm = TRUE))
std2 <- function(S) sweep(sweep(S, 2, ref2$mean[colnames(S)]), 2, ref2$sd[colnames(S)], "/")
D02    <- cbind(D0, as.data.frame(std2(S0a)))
Dlist2 <- lapply(seq_len(m_imp), function(i)
  cbind(Dlist[[i]], as.data.frame(std2(make_scores_alt(complete(imp, i))))))
ok_g     <- !is.na(D0$gender)                    # gender "other" is excluded everywhere
keep_alt <- function(y) keep_for(y) & ok_g
stopifnot(!anyNA(W_alt), all(c("FAC0", "FRQ0", "DEM0") %in% names(D02)))

x_k <- function(k, use_norep)                    # all k components of FACTORS, FIELD, FREQ, DEM
  c(unlist(lapply(alt_groups, function(p) paste0(p, k, "c", seq_len(k))), use.names = FALSE),
    if (use_norep) c("NOR1", "nerep_jaCizi"))
x_1 <- function(use_norep) c("FAC0", "FLD0", "FRQ0", "DEM0", if (use_norep) "NOR1")
# FACTORS-only k-component models: k oblimin components for FACTORS, benchmark spec elsewhere.
# FIELD: FLD1 + FLD2; FREQ: FRQ1 + frek_positiv; DEM: DEM1 + gender + top_journal.
x_fac3 <- function(use_norep)
  c("FAC3c1", "FAC3c2", "FAC3c3",               # 3 oblimin components of FACTORS
    blocks$FIELD,                                 # FLD1, FLD2
    blocks$FREQ,                                  # FRQ1, frek_positiv
    blocks$DEM,                                   # DEM1, gender, top_journal
    if (use_norep) blocks$NOREP)                  # NOR1, nerep_jaCizi
x_fac4 <- function(use_norep)
  c("FAC4c1", "FAC4c2", "FAC4c3", "FAC4c4",     # 4 oblimin components of FACTORS
    blocks$FIELD,                                 # FLD1, FLD2
    blocks$FREQ,                                  # FRQ1, frek_positiv
    blocks$DEM,                                   # DEM1, gender, top_journal
    if (use_norep) blocks$NOREP)                  # NOR1, nerep_jaCizi

adj_r2 <- function(r2, n, k) 1 - (1 - r2) * (n - 1) / (n - k)   # k = number of coefficients incl. intercept

## 10-fold cross-validated R^2 for one data frame.
## Uses base R (no caret dependency at prediction time) to keep the loop fast.
## Returns the mean and SD of fold-level R^2 (pred vs. obs on the held-out fold).
cv10_r2 <- function(y, x, d, seed = 42) {
  set.seed(seed)
  d_complete <- d[complete.cases(d[, c(y, x)]), ]
  n   <- nrow(d_complete)
  idx <- sample(seq_len(n))
  folds <- cut(seq_along(idx), breaks = 10, labels = FALSE)
  r2_folds <- vapply(1:10, function(f) {
    train <- d_complete[idx[folds != f], , drop = FALSE]
    test  <- d_complete[idx[folds == f], , drop = FALSE]
    m     <- lm(reformulate(x, y), data = train)
    pred  <- predict(m, newdata = test)
    obs   <- test[[y]]
    ss_res <- sum((obs - pred)^2)
    ss_tot <- sum((obs - mean(obs))^2)
    if (ss_tot == 0) return(NA_real_)
    1 - ss_res / ss_tot
  }, numeric(1))
  list(mean = mean(r2_folds, na.rm = TRUE), sd = stats::sd(r2_folds, na.rm = TRUE))
}

## Pool AIC/BIC across imputed data sets (arithmetic mean; Lydersen et al. 2010).
## Note: this is an approximation — no between-imputation variance is propagated.
## It is used here only for model ranking, not for inference.
pool_aic_bic <- function(fits)
  c(AIC_MI = mean(vapply(fits, function(f) f$aic, numeric(1))),
    BIC_MI = mean(vapply(fits, function(f) f$bic, numeric(1))))

## Pool 10-fold CV R^2 across imputed data sets (mean of per-imputation means).
pool_cv <- function(y, x, dlist, kp, seed = 42) {
  cv_vals <- vapply(dlist, function(d) cv10_r2(y, x, d[kp, ], seed = seed)$mean, numeric(1))
  c(cvR2_MI_mean = mean(cv_vals), cvR2_MI_sd = stats::sd(cv_vals))
}

run_alt <- function(x, lab, ys) {
  coefs <- list(); summ <- list()
  for (y in ys) {
    kp   <- keep_alt(y)
    fits <- lapply(Dlist2, function(d) fit_lm(y, x, d[kp, ]))
    pm   <- pool_fits(fits)
    r2m  <- mean(vapply(fits, function(f) f$r2, numeric(1)))     # mean R^2 over the imputed data sets
    pm$r2 <- r2m
    ab   <- pool_aic_bic(fits)                                    # pooled AIC and BIC (mean over imputations)
    cv   <- pool_cv(y, x, Dlist2, kp)                            # 10-fold CV R^2 (mean and SD over imputations)
    cc   <- fit_lm(y, x, D02); tc <- tidy_fit(cc)
    cv_cc <- cv10_r2(y, x, D02[kp, ])                            # CV on complete cases (all rows with y observed)
    coefs[[length(coefs) + 1]] <- base::rbind(cbind(outcome = y, model = paste("MI,", lab), pm),
                                              cbind(outcome = y, model = paste("complete case,", lab), tc))
    summ[[length(summ) + 1]] <- data.frame(
      outcome = y, model = lab, predictors = length(x),
      n_MI = fits[[1]]$n, R2_MI = r2m, adjR2_MI = adj_r2(r2m, fits[[1]]$n, fits[[1]]$k),
      AIC_MI = ab[["AIC_MI"]], BIC_MI = ab[["BIC_MI"]],
      cvR2_MI = cv[["cvR2_MI_mean"]], cvR2_MI_sd = cv[["cvR2_MI_sd"]],
      n_cc = cc$n, R2_cc = cc$r2, adjR2_cc = adj_r2(cc$r2, cc$n, cc$k),
      AIC_cc = cc$aic, BIC_cc = cc$bic,
      cvR2_cc = cv_cc$mean, cvR2_cc_sd = cv_cc$sd)
  }
  list(coef = do.call(base::rbind, coefs), summ = do.call(base::rbind, summ))
}

alt_out <- list()
for (use_norep in c(FALSE, TRUE)) {
  nl <- if (use_norep) "with NOREP" else "without NOREP"
  specs <- c(list(list(x_all(use_norep), paste("2 oblimin + single items (report),", nl), outcomes_all)),
             lapply(alt_k, function(k)
               list(x_k(k, use_norep), paste(k, "oblimin components,", nl), outcomes)),
             list(list(x_fac3(use_norep), paste("3 oblimin FACTORS + benchmark rest,", nl), outcomes)),
             list(list(x_fac4(use_norep), paste("4 oblimin FACTORS + benchmark rest,", nl), outcomes)),
             list(list(x_1(use_norep), paste("1 unrotated component,", nl), outcomes_all)))
  for (s in specs) {
    r <- run_alt(s[[1]], s[[2]], s[[3]])
    alt_out[[length(alt_out) + 1]] <- r
    cat("\n--- ", s[[2]], " (MI) ---\n", sep = "")
    show(r$coef[grepl("^MI", r$coef$model), ])
  }
}
alt_res  <- do.call(base::rbind, lapply(alt_out, function(r) r$coef))
alt_summ <- do.call(base::rbind, lapply(alt_out, function(r) r$summ))
cat("\n--- Fit comparison (MI: mean over imputed data sets; CC: complete cases) ---\n")
cat("    adjR2: adjusted R^2 (complexity-penalised)\n")
cat("    AIC/BIC: Gaussian-likelihood criteria averaged over imputed data sets (MI)\n")
cat("             or on the complete-case model; lower is better\n")
cat("    cvR2: 10-fold cross-validated R^2 (mean over imputed data sets and its SD)\n\n")
op <- options(width = 240)
num <- vapply(alt_summ, base::is.numeric, logical(1))
alt_print <- alt_summ; alt_print[num] <- lapply(alt_print[num], base::round, 3)
print(alt_print, row.names = FALSE)
options(op)
write.csv(alt_res,  file.path(out_dir, "alt_models.csv"),      row.names = FALSE)
write.csv(alt_summ, file.path(out_dir, "alt_fit_summary.csv"), row.names = FALSE)

cat("\n--- Variance inflation factors, models with NOREP (first imputed data set) ---\n")
for (k in alt_k) { cat("\n", k, " oblimin components\n", sep = ""); print(vif(Dlist2[[1]], x_k(k, TRUE)), row.names = FALSE) }
cat("\n3 oblimin FACTORS + benchmark rest\n"); print(vif(Dlist2[[1]], x_fac3(TRUE)), row.names = FALSE)
cat("\n4 oblimin FACTORS + benchmark rest\n"); print(vif(Dlist2[[1]], x_fac4(TRUE)), row.names = FALSE)
cat("\n1 unrotated component\n"); print(vif(Dlist2[[1]], x_1(TRUE)), row.names = FALSE)

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

saveRDS(list(blocks = block_res, joint = joint_res, wald = wald_res, alt = alt_res, alt_fit = alt_summ, alt_weights = W_alt, items = item_res,
             weights = W, ref = ref), file.path(out_dir, "regression_CORE.rds"))
cat("\nCSV files and regression_CORE.rds written to '", out_dir, "'.\n", sep = "")
