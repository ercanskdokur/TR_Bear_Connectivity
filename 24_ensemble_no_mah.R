## ============================================================================
## 24_ensemble_no_mah.R
## Project: TR_Bear_Connectivity
## Purpose: Drop the MAH algorithm from the ensemble (TSS = 0.16, essentially
##   random) and recompute Ensemble/W_MEAN and Ensemble/MEAN for the bear SDM
##   AND the conflict SDM, for present and all 18 future projections.
##
## Strategy:
##   1. Back up the original 8-algorithm ensembles (suffix _orig8alg.tif) once
##   2. Re-aggregate per-algorithm rasters across the kept 7 algorithms
##   3. W_MEAN uses TSS-weighted mean (weights from Evaluation_Table.txt)
##   4. MEAN uses unweighted mean
##   5. Compute per-pixel difference vs. original ensemble; save summary
##
## Outputs:
##   24_ensemble_diff_summary.csv  — per-scenario diff stats
##   24_eval_no_mah.csv            — 7-algorithm eval consensus
##   (Ensemble/W_MEAN and Ensemble/MEAN overwritten in both ENMTML result dirs)
## ============================================================================

suppressPackageStartupMessages({
  library(terra); library(dplyr); library(readr)
})

.tb_find_paths_R <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- a[grepl("--file=", a)]
  if (length(f)) {
    d <- dirname(normalizePath(sub("--file=", "", f[1]), mustWork = FALSE))
    if (file.exists(file.path(d, "00_paths.R"))) return(d)
  }
  env_dir <- Sys.getenv("TB_SCRIPTS", unset = "")
  if (nzchar(env_dir) && file.exists(file.path(env_dir, "00_paths.R"))) return(env_dir)
  if (file.exists("00_paths.R")) return(getwd())
  stop("Cannot find 00_paths.R")
}
setwd(.tb_find_paths_R())
source("00_paths.R"); source("00_helpers.R")
tb_log_init("24_ensemble_no_mah")

DROP_ALG    <- "MAH"
ALL_ALGS    <- c("BIO","BRT","GAM","GLM","MAH","MXD","RDF","SVM")
KEEP_ALGS   <- setdiff(ALL_ALGS, DROP_ALG)

## ---- Helpers ---------------------------------------------------------------
.read_eval <- function(enm_dir) {
  fn <- file.path(enm_dir, "Evaluation_Table.txt")
  if (!file.exists(fn)) { tb_log(sprintf("MISSING eval: %s", fn), "WARN"); return(NULL) }
  read.table(fn, header = TRUE, sep = "\t", check.names = FALSE)
}

.tss_weights <- function(eval_df, keep = KEEP_ALGS) {
  e <- eval_df[eval_df$Algorithm %in% keep, c("Algorithm", "TSS")]
  e <- e[!duplicated(e$Algorithm), ]
  w <- pmax(e$TSS, 0)
  names(w) <- e$Algorithm
  if (sum(w) > 0) w <- w / sum(w)
  w
}

## Locate algorithm raster file for present (Algorithm/<ALG>/<sp>.tif) or
## future (Projection/<scen>/<ALG>/<sp>.tif).
.find_alg_tif <- function(parent_dir, alg, sp) {
  paths <- c(
    file.path(parent_dir, "Algorithm", alg, paste0(sp, ".tif")),
    file.path(parent_dir, alg, paste0(sp, ".tif"))
  )
  for (p in paths) if (file.exists(p)) return(p)
  hits <- list.files(file.path(parent_dir, alg),
                     pattern = paste0("^", sp, "\\.tif$"),
                     full.names = TRUE, recursive = FALSE)
  if (length(hits)) return(hits[1])
  hits <- list.files(file.path(parent_dir, "Algorithm", alg),
                     pattern = paste0("^", sp, "\\.tif$"),
                     full.names = TRUE, recursive = FALSE)
  if (length(hits)) return(hits[1])
  NA_character_
}

## Locate ensemble file (Ensemble/<kind>/<sp>.tif).
.find_ens_tif <- function(parent_dir, kind, sp) {
  p <- file.path(parent_dir, "Ensemble", kind, paste0(sp, ".tif"))
  if (file.exists(p)) return(p)
  hits <- list.files(file.path(parent_dir, "Ensemble", kind),
                     pattern = paste0("^", sp, "\\.tif$"),
                     full.names = TRUE, recursive = FALSE)
  if (length(hits)) return(hits[1])
  NA_character_
}

## Recompute one scenario directory.
.recompute_one <- function(scen_dir, sp, tss_w, scen_label) {
  algo_rasters <- list()
  for (a in names(tss_w)) {
    f <- .find_alg_tif(scen_dir, a, sp)
    if (is.na(f)) { tb_log(sprintf("missing alg=%s for %s", a, scen_label), "WARN"); next }
    algo_rasters[[a]] <- terra::rast(f)
  }
  if (!length(algo_rasters)) return(NULL)
  stk <- terra::rast(algo_rasters)

  ## Weighted mean (W_MEAN)
  w_used <- tss_w[names(algo_rasters)]
  w_used <- w_used / sum(w_used)
  wmean <- sum(stk * w_used)

  ## Simple mean (MEAN)
  mn <- terra::app(stk, mean, na.rm = TRUE)

  ## Backup original ensembles once, then overwrite
  out <- list()
  for (kind in c("W_MEAN", "MEAN")) {
    orig <- .find_ens_tif(scen_dir, kind, sp)
    if (is.na(orig)) {
      tb_log(sprintf("orig %s ensemble missing for %s — skipping diff",
                      kind, scen_label), "WARN")
      next
    }
    bak <- sub("\\.tif$", "_orig8alg.tif", orig)
    if (!file.exists(bak)) file.copy(orig, bak)
    new_r <- if (kind == "W_MEAN") wmean else mn
    terra::writeRaster(new_r, orig, overwrite = TRUE)
    old_r <- terra::rast(bak)
    d <- new_r - old_r
    out[[kind]] <- data.frame(
      scenario = scen_label,
      ensemble = kind,
      mean_abs_diff = mean(abs(terra::values(d)), na.rm = TRUE),
      max_abs_diff  = max(abs(terra::values(d)), na.rm = TRUE),
      rmse          = sqrt(mean(terra::values(d)^2, na.rm = TRUE)),
      cor_old_new   = unname(cor(terra::values(old_r), terra::values(new_r),
                                  use = "pairwise.complete.obs"))
    )
  }
  do.call(rbind, out)
}

## Process one ENMTML output tree (bear or conflict).
.process_tree <- function(enm_root, sp, label) {
  tb_log_section(sprintf("%s — %s", label, enm_root))
  eval_df <- .read_eval(enm_root)
  if (is.null(eval_df)) return(NULL)
  tss_w <- .tss_weights(eval_df)
  tb_log(sprintf("[%s] TSS weights (no MAH): %s", label,
                  paste(sprintf("%s=%.3f", names(tss_w), tss_w), collapse = " ")))

  diffs <- list()
  ## Present
  diffs[["present"]] <- .recompute_one(enm_root, sp, tss_w, "present")

  ## Future scenarios
  proj_root <- file.path(enm_root, "Projection")
  if (dir.exists(proj_root)) {
    scens <- list.dirs(proj_root, recursive = FALSE)
    for (sd in scens) {
      sl <- basename(sd)
      diffs[[sl]] <- .recompute_one(sd, sp, tss_w, sl)
    }
  }
  diff_df <- do.call(rbind, diffs)
  diff_df$tree <- label
  diff_df
}

## ---- BEAR SDM --------------------------------------------------------------
diff_bear    <- .process_tree(TB_OUT_ENMTML,          "Ursus_arctos",          "bear")

## ---- CONFLICT SDM ----------------------------------------------------------
diff_conf    <- .process_tree(TB_OUT_ENMTML_CONFLICT, "Ursus_arctos_conflict", "conflict")

## ---- Save summaries --------------------------------------------------------
diff_all <- rbind(diff_bear, diff_conf)
tb_save_table(diff_all, "24_ensemble_diff_summary")

## Consensus 7-algorithm eval (mean per metric across kept algos)
eval_bear <- .read_eval(TB_OUT_ENMTML)
eval_conf <- .read_eval(TB_OUT_ENMTML_CONFLICT)
.eval_consensus <- function(eval_df, keep = KEEP_ALGS) {
  e <- eval_df[eval_df$Algorithm %in% keep, ]
  num_cols <- setdiff(names(e),
                      c("Sp","Algorithm","Threshold","Partition"))
  ag <- aggregate(e[, num_cols], by = list(Algorithm = e$Algorithm), FUN = mean,
                  na.rm = TRUE)
  ag
}
eval_no_mah <- rbind(
  cbind(tree = "bear",     .eval_consensus(eval_bear)),
  cbind(tree = "conflict", .eval_consensus(eval_conf)))
tb_save_table(eval_no_mah, "24_eval_no_mah")

tb_log_section("Summary")
print(diff_all)
tb_log_session()
tb_log("24_ensemble_no_mah DONE")
