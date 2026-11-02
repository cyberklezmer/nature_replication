# Convert all columns used in any group to numeric ordinal codes
encode_ordinal <- function(w, groups, levels_list = NULL) {
  cols <- unique(unlist(groups))
  for (cn in cols) {
    x <- w[[cn]]
    
    if (!is.null(levels_list) && !is.null(levels_list[[cn]])) {
      # explicit order given for this column
      x <- as.integer(factor(x, levels = levels_list[[cn]], ordered = TRUE))
      
    } else if (is.factor(x)) {
      # factor / ordered factor: position of the level (check order first!)
      x <- as.integer(x)
      
    } else if (is.character(x)) {
      num <- suppressWarnings(as.numeric(x))
      if (all(is.na(num[!is.na(x)]))) {
        stop("Column '", cn, "' is text without numeric codes. ",
             "Give its level order in levels_list.")
      }
      x <- num
    }
    w[[cn]] <- x
  }
  w
}

## PCA.analysis(): PCA of ordinal survey items based on polychoric correlations
## (missing values handled by pairwise deletion, no imputation), with PDF output:
##   1. cumulative-variance scree plot
##   2. eigenvalue scree plot (Kaiser line + parallel-analysis line)
##   3. loadings plot PC1 vs PC2 (unrotated)
##   4. loadings plot RC1 vs RC2 (varimax rotation)
##
## Arguments
##   w          data frame with the items (ordered factors or numeric codes)
##   groups     named list of character vectors with column names of w,
##              e.g. list(FACTORS = c("fakt_podvod", "fakt_karier", ...))
##              (a single unnamed character vector is treated as one group)
##   out_dir    folder for the PDF files
##
## The recommended number of components comes from parallel analysis.
##
## Value: invisibly, a list with one element per group (pca, eigenvalues,
##        cumvar, ncp_parallel, loadings, varimax, scores).

library(psych)
library(ggplot2)
library(ggrepel)

PCA.analysis <- function(w, groups, out_dir = "figures") {
  
  if (is.character(groups)) groups <- list(group = groups)
  if (is.null(names(groups))) names(groups) <- paste0("group", seq_along(groups))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  ## ---- helpers ------------------------------------------------------------
  to_num <- function(col) {
    if (is.factor(col)) return(as.integer(col))
    suppressWarnings(as.numeric(col))   # character/free text becomes all NA
  }
  
  # prcomp-like object from a polychoric correlation matrix
  polychoric_prcomp <- function(dat) {
    rng <- diff(range(as.matrix(dat), na.rm = TRUE)) + 1
    if (rng > 8) {
      # polychoric() refuses items with more than 8 categories; use Spearman
      message("  range of codes = ", rng, " (> 8): using Spearman correlations ",
              "instead of polychoric for this group")
      R <- cor(dat, method = "spearman", use = "pairwise.complete.obs")
    } else {
      R <- suppressWarnings(polychoric(dat, global = FALSE)$rho)
    }
    R   <- cor.smooth(R)
    eig <- eigen(R, symmetric = TRUE)
    rot <- eig$vectors
    dimnames(rot) <- list(colnames(dat), paste0("PC", seq_len(ncol(dat))))
    ctr <- colMeans(dat, na.rm = TRUE)
    sds <- apply(dat, 2, sd, na.rm = TRUE)
    x   <- scale(as.matrix(dat), center = ctr, scale = sds) %*% rot
    rownames(x) <- rownames(dat)
    list(pca = structure(list(sdev = sqrt(pmax(eig$values, 0)), rotation = rot,
                              center = ctr, scale = sds, x = x),
                         class = "prcomp"),
         R = R)
  }
  
  # make the largest absolute loading of every component positive
  fix_pca_signs <- function(p) {
    for (j in seq_len(ncol(p$rotation))) {
      idx <- which.max(abs(p$rotation[, j]))
      if (p$rotation[idx, j] < 0) {
        p$rotation[, j] <- -p$rotation[, j]
        if (!is.null(p$x)) p$x[, j] <- -p$x[, j]
      }
    }
    p
  }
  
  theme_pdf <- theme_bw(base_size = 14) +
    theme(plot.title = element_text(face = "bold"))
  
  plot_loadings <- function(L, xlab, ylab, title, file) {
    df <- data.frame(Dim1 = L[, 1], Dim2 = L[, 2], variable = rownames(L))
    p <- ggplot(df, aes(Dim1, Dim2, label = variable)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70") +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70") +
      geom_segment(aes(xend = 0, yend = 0), alpha = 0.4,
                   arrow = arrow(length = unit(0.2, "cm"))) +
      geom_point(size = 4, colour = "#377EB8") +
      geom_text_repel(size = 4.5, max.overlaps = 30, box.padding = 0.5,
                      point.padding = 0.3, segment.alpha = 0.6) +
      labs(x = xlab, y = ylab, title = title) +
      theme_pdf
    ggsave(file, p, width = 10, height = 8)
  }
  
  ## ---- main loop ----------------------------------------------------------
  results <- list()
  
  for (g in names(groups)) {
    cols <- groups[[g]]
    if (length(cols) == 0) {
      warning(g, ": group has no columns, skipped."); next
    }
    if (!all(cols %in% names(w))) {
      stop(g, ": columns not found in w: ",
           paste(setdiff(cols, names(w)), collapse = ", "))
    }
    dat  <- as.data.frame(lapply(w[, cols, drop = FALSE], to_num))
    rownames(dat) <- rownames(w)
    
    keep <- vapply(dat, function(x) isTRUE(sum(!is.na(x)) > 1 && sd(x, na.rm = TRUE) > 0),
                   logical(1))
    if (any(!keep)) {
      warning(g, ": dropped non-numeric, all-NA or constant columns: ",
              paste(names(dat)[!keep], collapse = ", "))
      dat <- dat[, keep, drop = FALSE]
    }
    p_items <- ncol(dat)
    if (p_items < 3) { warning(g, ": fewer than 3 usable items, skipped."); next }
    
    fit <- polychoric_prcomp(dat)
    pca <- fix_pca_signs(fit$pca)
    
    eig    <- pca$sdev^2
    pct    <- 100 * eig / sum(eig)
    cumvar <- cumsum(eig) / sum(eig)
    
    # parallel analysis (random-data eigenvalues from the same matrix size)
    pa <- NULL
    invisible(capture.output(
      pa <- tryCatch(
        suppressWarnings(suppressMessages(
          fa.parallel(fit$R, n.obs = max(sum(complete.cases(dat)), nrow(dat) / 2),
                      fa = "pc", n.iter = 50, plot = FALSE))),
        error = function(e) NULL)))
    ncp_pa <- if (!is.null(pa)) pa$ncomp else NA
    pa_eig <- if (!is.null(pa) && !is.null(pa$pc.sim)) pa$pc.sim else NULL
    
    tag <- g
    
    ## 1. cumulative variance scree plot
    d1 <- data.frame(dim = 0:p_items, cumvar = c(0, cumvar))
    p1 <- ggplot(d1, aes(dim, cumvar)) +
      geom_line(colour = "#1F4E79", linewidth = 0.9) +
      geom_point(colour = "#1F4E79", size = 3) +
      scale_x_continuous(breaks = d1$dim) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
      labs(x = "Component", y = "Explained variance ratio",
           title = paste("PCA scree plot:", g)) +
      theme_classic(base_size = 14) +
      theme(panel.grid.major = element_line(colour = "grey90"),
            plot.title = element_text(face = "bold"))
    ggsave(file.path(out_dir, paste0("PCA-scree-cumvar-", tag, ".pdf")),
           p1, width = 10, height = 6)
    
    ## 2. eigenvalue scree plot with Kaiser line and parallel analysis
    d2 <- data.frame(dim = seq_len(p_items), eigenvalue = eig)
    lab_obs <- "Observed data"
    lab_ran <- "Random data (parallel analysis)"
    p2 <- ggplot(d2, aes(dim)) +
      geom_line(aes(y = eigenvalue, colour = lab_obs)) +
      geom_point(aes(y = eigenvalue, colour = lab_obs), size = 3) +
      scale_x_continuous(breaks = d2$dim) +
      labs(x = "Component", y = "Eigenvalue", colour = NULL,
           title = paste("Eigenvalue scree plot:", g),
           subtitle = paste0("Parallel analysis suggests ", ncp_pa,
                             " component(s)")) +
      theme_pdf + theme(legend.position = "bottom")
    cols_leg <- setNames("#377EB8", lab_obs)
    if (!is.null(pa_eig)) {
      d2$random <- pa_eig[seq_len(p_items)]
      p2 <- p2 + geom_line(data = d2, aes(y = random, colour = lab_ran),
                           linetype = "dotted") +
        geom_point(data = d2, aes(y = random, colour = lab_ran),
                   shape = 4, size = 3)
      cols_leg <- c(cols_leg, setNames("firebrick", lab_ran))
    } else {
      # fallback only when parallel analysis failed: show the Kaiser line
      p2 <- p2 + geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
        labs(caption = "Dashed line: Kaiser criterion (eigenvalue = 1)")
    }
    p2 <- p2 + scale_colour_manual(values = cols_leg)
    ggsave(file.path(out_dir, paste0("PCA-scree-eigen-", tag, ".pdf")),
           p2, width = 8, height = 6)
    
    ## 3. unrotated loadings, PC1 vs PC2
    plot_loadings(pca$rotation[, 1:2, drop = FALSE],
                  paste0("PC1 (", round(pct[1], 1), "%)"),
                  paste0("PC2 (", round(pct[2], 1), "%)"),
                  paste("PCA variable loadings:", g),
                  file.path(out_dir, paste0("PCA-loadings-", tag, ".pdf")))
    
    ## 4. varimax-rotated two-component solution
    vm <- principal(fit$R, nfactors = 2, rotate = "varimax")
    L  <- unclass(vm$loadings)[, 1:2, drop = FALSE]
    vm_pct <- 100 * colSums(L^2) / p_items
    plot_loadings(L,
                  paste0("RC1 (", round(vm_pct[1], 1), "%)"),
                  paste0("RC2 (", round(vm_pct[2], 1), "%)"),
                  paste("Varimax-rotated loadings:", g),
                  file.path(out_dir, paste0("PCA-varimax-", tag, ".pdf")))
    
    ## Scores for the recommended number of components (NA if an item is missing)
    # (falls back to the Kaiser count only if parallel analysis failed)
    ncp_use <- if (is.na(ncp_pa)) max(sum(eig > 1), 1) else ncp_pa
    
    results[[g]] <- list(
      pca          = pca,
      eigenvalues  = eig,
      cumvar       = cumvar,
      ncp_parallel = ncp_pa,
      loadings     = pca$rotation[, 1:max(ncp_use, 2), drop = FALSE],
      varimax      = vm,
      scores       = pca$x[, 1:max(ncp_use, 1), drop = FALSE]
    )
    
    message(sprintf("%s: %d items, parallel analysis -> %s comp.",
                    g, p_items, ncp_pa))
  }
  
  invisible(results)
}