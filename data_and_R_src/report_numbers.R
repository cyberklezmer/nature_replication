## report_numbers.R
## Collects every number used in the LaTeX report, so nothing has to be read
## from plots, and builds the table of scales.
##
## Needs (already in your session or sourced below):
##   my        original data (factors with labels)
##   my.vars   numeric coded data used in the PCA
##   groups    named list of variable names per group
##   PCA_analysis.R (updated version) and split_half_congruence.R
##
## Run it and paste the printed output (or attach the CSV files written to
## the folder "report_numbers").

library(psych)

source("PCA_analysis.R")
source("split_half_congruence.R")

dir.create("report_numbers", showWarnings = FALSE)

## ---- run the analysis for all groups (also rewrites the PDF figures) -------
res <- PCA.analysis(my.vars, groups)

## ---- A. retention numbers ---------------------------------------------------
retention <- do.call(rbind, lapply(names(res), function(g) {
  r   <- res[[g]]
  eig <- r$eigenvalues
  ran <- r$random_eigenvalues
  rn  <- function(i) if (is.null(ran)) NA_real_ else ran[i]
  # number of components read from the plotted values: retain while observed > random
  pa_plot <- if (is.null(ran)) NA else {
    k <- which(!(eig > ran[seq_along(eig)]))
    if (length(k)) k[1] - 1 else length(eig)
  }
  data.frame(
    group = g, items = r$n_items, n_total = r$n_total,
    n_complete = r$n_complete, n_obs_pa = r$n_obs_pa,
    correlation = r$correlation,
    pa_ncomp = r$ncp_parallel, pa_from_plotted_values = pa_plot,
    eig1 = eig[1], eig2 = eig[2], eig3 = eig[3],
    rand1 = rn(1), rand2 = rn(2), rand3 = rn(3),
    pc1_pct = 100 * eig[1] / sum(eig), pc2_pct = 100 * eig[2] / sum(eig),
    var_2comp_pct = 100 * sum(eig[1:2]) / sum(eig),
    rc1_pct = unname(r$varimax_pct[1]), rc2_pct = unname(r$varimax_pct[2]),
    ncomp_for_90pct = which(r$cumvar >= 0.9 - 1e-8)[1],
    row.names = NULL
  )
}))
cat("\n=== A. Retention ===\n")
print(retention, digits = 3, row.names = FALSE)
write.csv(retention, "report_numbers/A_retention.csv", row.names = FALSE)

## ---- B. varimax loadings ----------------------------------------------------
varimax_tab <- do.call(rbind, lapply(names(res), function(g) {
  L <- unclass(res[[g]]$varimax$loadings)[, 1:2, drop = FALSE]
  data.frame(group = g, variable = rownames(L), RC1 = L[, 1], RC2 = L[, 2],
             row.names = NULL)
}))
cat("\n=== B. Varimax loadings ===\n")
print(varimax_tab, digits = 3, row.names = FALSE)
write.csv(varimax_tab, "report_numbers/B_varimax.csv", row.names = FALSE)

## ---- C. oblique (oblimin) pattern and structure, component correlation ------
ob <- lapply(res, function(r) principal(r$R, nfactors = 2, rotate = "oblimin"))
oblique_tab <- do.call(rbind, lapply(names(ob), function(g) {
  o <- ob[[g]]
  P <- unclass(o$loadings)
  S <- unclass(o$Structure)
  do.call(rbind, lapply(colnames(P), function(cn)
    data.frame(group = g, variable = rownames(P), component = cn,
               pattern = P[, cn], structure = S[, cn], phi = o$Phi[1, 2],
               row.names = NULL)))
}))
cat("\n=== C. Oblique (oblimin): pattern, structure, component correlation phi ===\n")
print(oblique_tab, digits = 3, row.names = FALSE)
write.csv(oblique_tab, "report_numbers/C_oblique.csv", row.names = FALSE)

## ---- D. split-half stability for all groups (may take a few minutes) --------
sh <- lapply(groups, function(cols) split_half_congruence(my.vars, cols))
split_tab <- do.call(rbind, lapply(names(sh), function(g)
  data.frame(group = g, component = rownames(sh[[g]]$summary),
             sh[[g]]$summary, row.names = NULL)))
cat("\n=== D. Split-half congruence (half A vs half B, 500 splits) ===\n")
print(split_tab, row.names = FALSE)
write.csv(split_tab, "report_numbers/D_splithalf.csv", row.names = FALSE)

## ---- E. table of scales ------------------------------------------------------
## For every variable: how the original answers (my) map to the codes used in
## the analysis (my.vars), number of valid and missing answers.
scale_table <- function(my, my.vars, groups) {
  do.call(rbind, lapply(names(groups), function(g) {
    do.call(rbind, lapply(groups[[g]], function(v) {
      lab  <- as.character(my[[v]])
      code <- my.vars[[v]]
      tab  <- table(label = lab, code = code, useNA = "ifany")
      cds  <- colnames(tab)
      map  <- vapply(seq_len(ncol(tab)), function(j) {
        labs <- rownames(tab)[tab[, j] > 0]
        paste(ifelse(is.na(labs), "NA", labs), collapse = " / ")
      }, character(1))
      cds[is.na(cds)] <- "NA"
      ok <- cds != "NA"
      ord <- order(suppressWarnings(as.numeric(cds[ok])))
      data.frame(
        group = g, variable = v,
        n_valid = sum(!is.na(code)), n_missing = sum(is.na(code)),
        n_categories = sum(ok),
        lowest  = paste0(cds[ok][ord][1], " = ", map[ok][ord][1]),
        highest = paste0(cds[ok][ord][length(ord)], " = ", map[ok][ord][length(ord)]),
        all_codes = paste0(cds, " = ", map, collapse = "; "),
        stringsAsFactors = FALSE)
    }))
  }))
}
scales <- scale_table(my, my.vars, groups)
rownames(scales) <- NULL
cat("\n=== E. Scales ===\n")
print(scales[, c("group", "variable", "n_valid", "n_missing", "n_categories",
                 "lowest", "highest")], row.names = FALSE)
cat("\nFull code-to-label mapping:\n")
print(scales[, c("variable", "all_codes")], row.names = FALSE)
write.csv(scales, "report_numbers/E_scales.csv", row.names = FALSE)
