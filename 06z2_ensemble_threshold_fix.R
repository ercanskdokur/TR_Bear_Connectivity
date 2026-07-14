## ============================================================================
## 06z2_ensemble_threshold_fix.R
## Project: TR_Bear_Connectivity — recompute Ensemble MAX_TSS thresholds + binaries
## Background: 06z_ensemble_manual.R computed continuous Ensembles successfully
##   but binary thresholds came back NA because Occurrences_Fitting.txt contains
##   only presences (ENMTML did not persist pseudo-absences separately for BOOT).
## Fix: sample random background points from the species' accessible-area mask
##   (Extent_Masks/<sp>.tif) and use them as pseudo-absences to compute the
##   MAX_TSS threshold for each ensemble flavor. Then write the binary rasters
##   (overwriting the NA placeholders the previous pass produced).
## ============================================================================

suppressPackageStartupMessages({
  library(terra)
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
tb_log_init("06z2_ensemble_threshold_fix")

SP <- "Ursus_arctos"
N_BG <- 10000L
set.seed(42L)

mask_file <- file.path(TB_OUT_ENMTML, "Extent_Masks", paste0(SP, ".tif"))
occ_file  <- file.path(TB_OUT_ENMTML, "Occurrences_Cleaned.txt")

tb_log(sprintf("mask  = %s", mask_file))
tb_log(sprintf("occ   = %s", occ_file))

mask_r <- terra::rast(mask_file)
occ    <- read.table(occ_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
pres_mat <- as.matrix(occ[, c("x", "y")])
tb_log(sprintf("presences cleaned: %d", nrow(pres_mat)))

## Sample N_BG random background points from non-NA cells in mask
## terra::spatSample() with raster: returns data.frame with x,y if xy=TRUE
bg_xy <- terra::spatSample(mask_r, size = N_BG, method = "random",
                           na.rm = TRUE, xy = TRUE, values = FALSE)
bg_mat <- as.matrix(bg_xy[, c("x", "y")])
tb_log(sprintf("background sampled: %d (from accessible-area mask)", nrow(bg_mat)))

## ---- MAX_TSS threshold helper -----------------------------------------------
.max_tss_threshold <- function(r, pres_mat, bg_mat) {
  vp <- terra::extract(r, pres_mat)
  va <- terra::extract(r, bg_mat)
  ## extract() with matrix returns vector OR single-column data.frame depending on version
  if (is.data.frame(vp)) vp <- vp[, 1]
  if (is.data.frame(va)) va <- va[, 1]
  vp <- vp[!is.na(vp)]; va <- va[!is.na(va)]
  if (length(vp) < 5 || length(va) < 5) return(NA_real_)
  cand <- unique(c(vp, va))
  if (length(cand) > 500) cand <- quantile(cand, probs = seq(0.005, 0.995, length.out = 500), names = FALSE)
  cand <- sort(unique(cand))
  ts <- vapply(cand, function(thr) {
    tpr <- mean(vp >= thr); tnr <- mean(va <  thr)
    tpr + tnr - 1
  }, numeric(1))
  cand[which.max(ts)]
}

## ---- Present ensembles -------------------------------------------------------
tb_log_section("PRESENT")
ens_mean  <- terra::rast(file.path(TB_OUT_ENMTML, "Ensemble", "MEAN",   paste0(SP, ".tif")))
ens_wmean <- terra::rast(file.path(TB_OUT_ENMTML, "Ensemble", "W_MEAN", paste0(SP, ".tif")))

thr_m <- .max_tss_threshold(ens_mean,  pres_mat, bg_mat)
thr_w <- .max_tss_threshold(ens_wmean, pres_mat, bg_mat)
tb_log(sprintf("MEAN   MAX_TSS = %.6f", thr_m))
tb_log(sprintf("W_MEAN MAX_TSS = %.6f", thr_w))

terra::writeRaster(ens_mean  >= thr_m,
                   file.path(TB_OUT_ENMTML, "Ensemble", "MEAN",   "MAX_TSS", paste0(SP, ".tif")),
                   overwrite = TRUE)
terra::writeRaster(ens_wmean >= thr_w,
                   file.path(TB_OUT_ENMTML, "Ensemble", "W_MEAN", "MAX_TSS", paste0(SP, ".tif")),
                   overwrite = TRUE)
tb_log("present binary written")

## ---- Future scenarios --------------------------------------------------------
tb_log_section("FUTURE 18 SCENARIOS")
proj_dir <- file.path(TB_OUT_ENMTML, "Projection")
scenarios <- list.dirs(proj_dir, recursive = FALSE)
tb_tic()
for (scn in scenarios) {
  f_m <- file.path(scn, "Ensemble", "MEAN",   paste0(SP, ".tif"))
  f_w <- file.path(scn, "Ensemble", "W_MEAN", paste0(SP, ".tif"))
  if (!file.exists(f_m) || !file.exists(f_w)) {
    tb_log(sprintf("WARN: missing ensemble files in %s", basename(scn)), "WARN")
    next
  }
  r_m <- terra::rast(f_m); r_w <- terra::rast(f_w)
  terra::writeRaster(r_m >= thr_m,
                     file.path(scn, "Ensemble", "MEAN",   "MAX_TSS", paste0(SP, ".tif")),
                     overwrite = TRUE)
  terra::writeRaster(r_w >= thr_w,
                     file.path(scn, "Ensemble", "W_MEAN", "MAX_TSS", paste0(SP, ".tif")),
                     overwrite = TRUE)
}
tb_toc(sprintf("future binaries (%d scn x 2 methods)", length(scenarios)))

## ---- Overwrite Thresholds_Ensemble.txt --------------------------------------
thr_table <- data.frame(
  Sp        = SP,
  Ensemble  = c("MEAN", "W_MEAN"),
  THR       = "MAX_TSS",
  THR_VALUE = c(thr_m, thr_w),
  stringsAsFactors = FALSE
)
write.table(thr_table, file.path(TB_OUT_ENMTML, "Thresholds_Ensemble.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)
tb_log(sprintf("wrote: Thresholds_Ensemble.txt (%.4f / %.4f)", thr_m, thr_w))

tb_log_session()
tb_log("06z2_ensemble_threshold_fix DONE")
