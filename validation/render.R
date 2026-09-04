#!/usr/bin/env Rscript
# One-stop renderer for the TargetedBeast operator-validation artefacts, so the README can embed
# them as static figures that regenerate on their own:
#   1. detailed_balance.md  (written by the detailed-balance test suites)  -> detailed_balance.png
#   2. rateprior_logs/*.log (written by RatePriorMCMCTest)                  -> rateprior_<stat>.png
#
# Usage:  Rscript validation/render.R [validationDir] [ratePriorLogDir]
suppressMessages({library(gridExtra); library(grid); library(data.table)})

args   <- commandArgs(trailingOnly = TRUE)
VALDIR <- if (length(args) >= 1) args[1] else "validation"
LOGDIR <- if (length(args) >= 2) args[2] else {
  a <- file.path(VALDIR, "rateprior_logs"); b <- "/tmp/rateprior_logs"
  if (dir.exists(a) && length(list.files(a, "\\.log$"))) a else b
}

## ---------- 1. detailed-balance matrix (test x bucket-type -> pass/FAIL/vacuous) -> PNG + md ----------
# Built from the per-pair CSV: each mapper (bucket type) is its own column, so you can see per test
# which bucket dimensions balanced, which never moved (vacuous), and which are broken (NaN/FAIL).
render_matrix <- function(csv, png, mdout = NULL) {
  if (!file.exists(csv)) { cat("skip table: no", csv, "\n"); return(invisible()) }
  d <- fread(csv); d[, compared := as.logical(compared)]
  res <- d[, {
    cmp <- .SD[compared == TRUE]
    r <- if (nrow(cmp) == 0) "vacuous"
         else if (any(!is.finite(cmp$forward) | !is.finite(cmp$backward) | !is.finite(cmp$tol))) "NaN"
         else if (any(abs(cmp$forward - cmp$backward) > cmp$tol)) "FAIL"
         else "pass"
    list(result = r)
  }, by = .(suite, test, mapper)]
  w <- dcast(res, suite + test ~ mapper, value.var = "result", fill = "-")
  mcols <- setdiff(names(w), c("suite", "test"))
  pref  <- c("TreeLength", "TreeImbalance", "RootFirstChildHeight", "RateLogSD", "Rate0", "Down", "Imbalance")
  mcols <- c(intersect(pref, mcols), setdiff(mcols, pref))
  w <- w[, c("suite", "test", mcols), with = FALSE]; setorder(w, suite, test)
  d2 <- as.data.frame(w, check.names = FALSE)

  cellcol <- function(v) ifelse(v == "pass", "#d9f0d3",
                         ifelse(v == "FAIL", "#f4a9a0",
                         ifelse(v %in% c("vacuous", "NaN"), "#ffe6b3", "#f5f5f5")))
  fill <- matrix("white", nrow(d2), ncol(d2))
  fill[, 1] <- "#ececec"; fill[, 2] <- "#f7f7f7"
  for (j in seq_along(mcols)) fill[, 2 + j] <- cellcol(d2[[2 + j]])

  th <- ttheme_minimal(base_size = 9,
    core    = list(bg_params = list(fill = fill, col = "grey80", lwd = 0.6),
                   fg_params = list(hjust = 0, x = 0.08, fontsize = 8)),
    colhead = list(bg_params = list(fill = "grey25", col = "grey25"),
                   fg_params = list(col = "white", hjust = 0, x = 0.05, fontface = "bold", fontsize = 8)))
  g <- tableGrob(d2, rows = NULL, theme = th)
  wi <- convertWidth(sum(g$widths), "in", valueOnly = TRUE) + 0.3
  hi <- convertHeight(sum(g$heights), "in", valueOnly = TRUE) + 0.3
  png(png, width = wi, height = hi, units = "in", res = 220, type = "cairo", bg = "white")
  grid.newpage(); grid.draw(g); invisible(dev.off())
  cat("wrote", png, "\n")

  if (!is.null(mdout)) {
    ln <- c(paste0("| ", paste(names(d2), collapse = " | "), " |"),
            paste0("|", paste(rep("---", ncol(d2)), collapse = "|"), "|"),
            apply(d2, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
    writeLines(ln, mdout); cat("wrote", mdout, "\n")
  }
}

## ---------- 2. prior-invariance sampled distributions -> PNG per stat ----------
render_distributions <- function(logdir, outdir) {
  base <- file.path(logdir, "ratePrior_BASELINE.log")
  if (!file.exists(base)) { cat("skip distributions: no logs in", logdir, "\n"); return(invisible()) }
  read_log <- function(f) tryCatch(fread(f, skip = "Sample", header = TRUE), error = function(e) NULL)
  b <- read_log(base)
  STATS <- intersect(c("ucldStdev", "ucldMean", "tree.height", "tree.treeLength", "tree.imbalance"), names(b))
  files <- list.files(logdir, "^ratePrior_.*\\.log$", full.names = TRUE)
  runs <- list()
  for (f in files) {
    k <- sub("^ratePrior_", "", sub("\\.log$", "", basename(f)))
    if (k == "BASELINE") next
    d <- read_log(f); if (is.null(d) || nrow(d) < 10) next
    runs[[k]] <- d[floor(0.2 * nrow(d)):nrow(d)]
  }
  nm <- names(runs); ord <- c(intersect("STD", nm), sort(setdiff(nm, c("STD", "ALL"))), intersect("ALL", nm))
  runs <- runs[ord]
  for (st in STATS) {
    bx <- b[[st]]; rng <- quantile(bx, c(0.002, 0.998), na.rm = TRUE)
    bd <- density(bx, from = rng[1], to = rng[2], n = 512)
    n <- length(runs); nc <- 5; nr <- ceiling(n / nc)
    png(file.path(outdir, paste0("rateprior_", gsub("\\.", "_", st), ".png")),
        width = 15, height = 2.4 * nr + 0.6, units = "in", res = 220, type = "cairo", bg = "white")
    par(mfrow = c(nr, nc), mar = c(2.4, 2.2, 2, 0.6), mgp = c(1.3, 0.5, 0), oma = c(0, 0, 2.2, 0))
    for (k in names(runs)) {
      x <- runs[[k]][[st]]; x <- x[is.finite(x)]
      dd <- density(x, from = rng[1], to = rng[2], n = 512)
      reld <- abs(mean(x) - mean(bx)) / abs(mean(bx))
      col <- if (reld > 0.05) "firebrick" else if (k == "ALL") "purple" else if (k == "STD") "grey30" else "navy"
      plot(bd, xlim = rng, ylim = c(0, max(bd$y, dd$y) * 1.05), main = "", xlab = "", ylab = "", col = NA)
      polygon(c(bd$x[1], bd$x, bd$x[512]), c(0, bd$y, 0), col = "grey85", border = "grey70")
      lines(dd, col = col, lwd = 2)
      title(main = sprintf("%s (%.0f%%)", k, 100 * reld), col.main = col, cex.main = 0.85)
      abline(v = mean(bx), col = "grey55", lty = 3); abline(v = mean(x), col = col, lty = 2)
    }
    mtext(sprintf("%s  -  grey = DirectSimulator truth; line = operator MCMC (red = |Delta mean|>5%%)", st),
          outer = TRUE, cex = 1.0, font = 2, line = 0.4)
    invisible(dev.off())
    cat("wrote", file.path(outdir, paste0("rateprior_", gsub("\\.", "_", st), ".png")), "\n")
  }
}

## ---------- 3. all detailed-balance buckets -> forward-vs-backward scatter per test ----------
render_balance_scatter <- function(csv, png) {
  if (!file.exists(csv)) { cat("skip buckets: no", csv, "\n"); return(invisible()) }
  d <- fread(csv); d[, compared := as.logical(compared)]

  # Stable colour per bucket type (mapper), consistent across every panel and shared by one legend.
  pref   <- c("TreeLength", "TreeImbalance", "RootFirstChildHeight", "RateLogSD", "Rate0", "Down", "Imbalance")
  mtypes <- unique(d$mapper); mtypes <- c(intersect(pref, mtypes), sort(setdiff(mtypes, pref)))
  okabe  <- c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#56B4E9", "#D55E00", "#F0E442", "#666666", "#000000")
  cols   <- if (length(mtypes) <= length(okabe)) okabe[seq_along(mtypes)] else hcl.colors(length(mtypes), "Dark3")
  pal    <- setNames(cols, mtypes)

  tests <- unique(d$test)
  n <- length(tests); nc <- 5; nr <- ceiling((n + 1) / nc)   # +1 cell reserved for the shared legend
  png(png, width = 15, height = 2.9 * nr + 0.4, units = "in", res = 220, type = "cairo", bg = "white")
  par(mfrow = c(nr, nc), mar = c(3, 3, 2, 0.6), mgp = c(1.7, 0.5, 0))
  for (t in tests) {
    s <- d[test == t]
    cmp <- s[compared == TRUE]; vac <- s[compared == FALSE]
    fin <- cmp[is.finite(forward) & is.finite(backward)]
    rng <- if (nrow(fin)) range(c(fin$forward, fin$backward, 0)) else c(0, 1)
    plot(NA, xlim = rng, ylim = rng, xlab = "forward flow", ylab = "backward flow", main = t, cex.main = 0.8)
    abline(0, 1, col = "grey60", lty = 2)
    points(vac$forward, vac$backward, col = "grey75", pch = 1, cex = 0.6)          # vacuous buckets
    points(fin$forward, fin$backward, col = pal[fin$mapper], pch = 19, cex = 0.7)  # hit + balanced, coloured by type
    bad <- fin[abs(forward - backward) > tol]                                      # off-diagonal beyond 3σ
    if (nrow(bad)) points(bad$forward, bad$backward, col = "firebrick", pch = 1, cex = 1.7, lwd = 2)  # ring keeps the type colour visible
    note <- if (nrow(fin) == 0) sprintf("%d hit (NaN flow!)", nrow(cmp))
            else sprintf("%d hit / %d vacuous%s", nrow(fin), nrow(vac), if (nrow(bad)) paste0(" / ", nrow(bad), " FAIL") else "")
    legend("topleft", note, bty = "n", cex = 0.72, text.col = if (nrow(fin) == 0) "firebrick" else "grey30")
  }
  # Shared legend in the reserved cell: bucket types, plus the vacuous / FAIL glyphs.
  plot.new()
  legend("center", title = "bucket type (colour)",
         legend = c(names(pal), "vacuous (low count)", "FAIL (|Δflow| > 3σ)"),
         col    = c(unname(pal), "grey75", "firebrick"),
         pch    = c(rep(19, length(pal)), 1, 1),
         pt.cex = c(rep(1.1, length(pal)), 0.9, 1.5), pt.lwd = 2,
         bty = "n", cex = 0.95)
  invisible(dev.off())
  cat("wrote", png, "\n")
}

dir.create(VALDIR, showWarnings = FALSE, recursive = TRUE)
render_matrix(file.path(VALDIR, "detailed_balance_pairs.csv"), file.path(VALDIR, "detailed_balance.png"),
              file.path(VALDIR, "detailed_balance.md"))
render_balance_scatter(file.path(VALDIR, "detailed_balance_pairs.csv"), file.path(VALDIR, "detailed_balance_buckets.png"))
render_distributions(LOGDIR, VALDIR)
