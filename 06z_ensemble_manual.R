## ============================================================================
## 06z_ensemble_manual.R
## Project: TR_Bear_Connectivity — manual post-hoc ENSEMBLE step
## Purpose: ENMTML 06_enmtml_run.R completed all 10 BOOT replicates + 8 algos +
##   18 future projections, but the Ensemble step was silently skipped because
##   `ensemble = list(...)` argument was passed (ENMTML expects named vector).
##   Pipeline produced per-algorithm rasters already averaged across replicates
##   (Algorithm/<algo>/Ursus_arctos.tif) and per-scenario per-algo projections
##   (Projection/<scn>/<algo>/Ursus_arctos.tif). This script recomputes the
##   two ensemble flavors ENMTML would have produced:
##     - MEAN: simple pixel-wise mean across the 8 algorithms
##     - W_MEAN: TSS-weighted pixel-wise mean (weights from Evaluation_Table.txt)
##   Continuous rasters + binary rasters (thresholded with ensemble's own MAX_TSS).
##
##   Mathematical equivalence: ENMTML's Ensemble_TMLA does exactly this
##   weighted averaging — manual implementation is identical and defensible.
## Outputs:
##   Ensemble/MEAN/Ursus_arctos.tif         (continuous, present)
##   Ensemble/W_MEAN/Ursus_arctos.tif       (continuous, present)
##   Ensemble/MEAN/MAX_TSS/Ursus_arctos.tif (binary, present)
##   Ensemble/W_MEAN/MAX_TSS/Ursus_arctos.tif
##   Projection/<scn>/Ensemble/MEAN/Ursus_arctos.tif        (continuous, future)
##   Projection/<scn>/Ensemble/W_MEAN/Ursus_arctos.tif
##   Projection/<scn>/Ensemble/MEAN/MAX_TSS/Ursus_arctos.tif (binary, future)
##   Projection/<scn>/Ensemble/W_MEAN/MAX_TSS/Ursus_arctos.tif
##   Thresholds_Ensemble.txt (ensemble-level MAX_TSS thresholds)
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
cat(sprintf("[bootstrap] wd = %s\n", getwd()))
source("00_paths.R"); source("00_helpers.R")
tb_log_init("06z_ensemble_manual")

SP <- "Ursus_arctos"
ALGOS <- TB_ENM_ALGORITHMS
tb_log(sprintf("species  = %s", SP))
tb_log(sprintf("algos    = %s", paste(ALGOS, collapse = ", ")))
tb_log(sprintf("result   = %s", TB_OUT_ENMTML))

## ============================================================================
## 1. Load TSS weights from Evaluation_Table.txt
## ============================================================================
tb_log_section("1. LOAD TSS WEIGHTS")

eval_tab <- read.table(file.path(TB_OUT_ENMTML, "Evaluation_Table.txt"),
                       header = TRUE, sep = "\t", stringsAsFactors = FALSE)
tb_log(sprintf("Evaluation_Table.txt rows=%d", nrow(eval_tab)))

tss_w <- setNames(eval_tab$TSS[match(ALGOS, eval_tab$Algorithm)], ALGOS)
if (any(is.na(tss_w))) {
  miss <- ALGOS[is.na(tss_w)]
  tb_log(sprintf("WARN: algos missing TSS in Eval_Table: %s", paste(miss, collapse=",")), "WARN")
}
tss_w <- pmax(tss_w, 0, na.rm = TRUE)   # clamp negative to 0
tb_log(sprintf("TSS weights: %s", paste(sprintf("%s=%.3f", names(tss_w), tss_w), collapse=", ")))

## ============================================================================
## 2. Ensemble helper functions
## ============================================================================
.ensemble_mean <- function(file_paths) {
  stk <- terra::rast(file_paths)
  terra::app(stk, fun = mean, na.rm = TRUE)
}

.ensemble_wmean <- function(file_paths, weights) {
  stk <- terra::rast(file_paths)
  if (length(weights) != terra::nlyr(stk)) stop("weights length mismatch")
  ws <- sum(weights)
  if (ws <= 0) stop("zero weight sum")
  # Per-pixel: sum(w_i * x_i) / sum(w_i)  (na ignored per layer)
  weighted <- stk * weights
  num <- terra::app(weighted, fun = sum, na.rm = TRUE)
  num / ws
}

## Compute MAX_TSS threshold for ensemble: given continuous raster + presences
## + pseudo-absences (from Occurrences_Cleaned + ENMTML's PA), choose threshold
## maximizing (sensitivity + specificity - 1)
.ensemble_threshold <- function(rast_ens, pres_xy, abs_xy) {
  vals_p <- terra::extract(rast_ens, pres_xy)[, 2]
  vals_a <- terra::extract(rast_ens, abs_xy)[, 2]
  vals_p <- vals_p[!is.na(vals_p)]
  vals_a <- vals_a[!is.na(vals_a)]
  if (length(vals_p) < 5 || length(vals_a) < 5) return(NA_real_)
  # Candidate thresholds: unique sorted continuous values
  cand <- sort(unique(c(vals_p, vals_a)))
  if (length(cand) > 200) cand <- quantile(cand, probs = seq(0.01, 0.99, length.out = 200), names = FALSE)
  tss_vals <- vapply(cand, function(thr) {
    tpr <- mean(vals_p >= thr)            # sensitivity
    tnr <- mean(vals_a <  thr)            # specificity
    tpr + tnr - 1
  }, numeric(1))
  cand[which.max(tss_vals)]
}

## ============================================================================
## 3. Present ensembles
## ============================================================================
tb_log_section("3. PRESENT ENSEMBLES")

present_files <- vapply(ALGOS, function(a) {
  file.path(TB_OUT_ENMTML, "Algorithm", a, paste0(SP, ".tif"))
}, character(1))
missing <- present_files[!file.exists(present_files)]
if (length(missing) > 0) stop(sprintf("Missing present files: %s", paste(missing, collapse=", ")))
tb_log(sprintf("present files: %d/%d found", sum(file.exists(present_files)), length(present_files)))

dir.create(file.path(TB_OUT_ENMTML, "Ensemble", "MEAN"),   recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(TB_OUT_ENMTML, "Ensemble", "W_MEAN"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(TB_OUT_ENMTML, "Ensemble", "MEAN",   "MAX_TSS"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(TB_OUT_ENMTML, "Ensemble", "W_MEAN", "MAX_TSS"), recursive = TRUE, showWarnings = FALSE)

tb_tic()
ens_mean_present  <- .ensemble_mean(present_files)
ens_wmean_present <- .ensemble_wmean(present_files, weights = unname(tss_w))
out_mean_p  <- file.path(TB_OUT_ENMTML, "Ensemble", "MEAN",   paste0(SP, ".tif"))
out_wmean_p <- file.path(TB_OUT_ENMTML, "Ensemble", "W_MEAN", paste0(SP, ".tif"))
terra::writeRaster(ens_mean_present,  out_mean_p,  overwrite = TRUE)
terra::writeRaster(ens_wmean_present, out_wmean_p, overwrite = TRUE)
tb_toc("present ensemble (MEAN + W_MEAN)")
tb_log(sprintf("wrote: %s", out_mean_p))
tb_log(sprintf("wrote: %s", out_wmean_p))

## ============================================================================
## 4. Future projection ensembles (18 scenarios)
## ============================================================================
tb_log_section("4. FUTURE PROJECTION ENSEMBLES")

proj_dir <- file.path(TB_OUT_ENMTML, "Projection")
scenarios <- list.dirs(proj_dir, recursive = FALSE)
tb_log(sprintf("scenarios: %d", length(scenarios)))

tb_tic()
for (scn in scenarios) {
  scn_name <- basename(scn)
  scn_files <- vapply(ALGOS, function(a) {
    file.path(scn, a, paste0(SP, ".tif"))
  }, character(1))
  exist_ok <- file.exists(scn_files)
  if (!all(exist_ok)) {
    tb_log(sprintf("WARN %s: missing %d algo files (skipping)", scn_name, sum(!exist_ok)), "WARN")
    next
  }
  out_dir_m  <- file.path(scn, "Ensemble", "MEAN")
  out_dir_w  <- file.path(scn, "Ensemble", "W_MEAN")
  out_dir_mb <- file.path(out_dir_m, "MAX_TSS")
  out_dir_wb <- file.path(out_dir_w, "MAX_TSS")
  for (d in c(out_dir_m, out_dir_w, out_dir_mb, out_dir_wb)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
  ens_m  <- .ensemble_mean(scn_files)
  ens_w  <- .ensemble_wmean(scn_files, weights = unname(tss_w))
  terra::writeRaster(ens_m, file.path(out_dir_m, paste0(SP, ".tif")), overwrite = TRUE)
  terra::writeRaster(ens_w, file.path(out_dir_w, paste0(SP, ".tif")), overwrite = TRUE)
}
tb_toc(sprintf("future ensemble (%d scenarios x 2 methods)", length(scenarios)))

## ============================================================================
## 5. Ensemble-level MAX_TSS thresholds
## ============================================================================
tb_log_section("5. ENSEMBLE THRESHOLDS")

## Recover presence + pseudo-absence used in modelling.
## ENMTML wrote Occurrences_Fitting.txt = presences + PA actually used.
occ_fit_file <- file.path(TB_OUT_ENMTML, "Occurrences_Fitting.txt")
if (!file.exists(occ_fit_file)) {
  tb_log(sprintf("WARN: %s missing — skip ensemble threshold computation", occ_fit_file), "WARN")
  thr_table <- data.frame()
} else {
  occ_fit <- read.table(occ_fit_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  tb_log(sprintf("Occurrences_Fitting rows=%d cols=%s", nrow(occ_fit), paste(names(occ_fit), collapse=",")))
  ## detect PresAbse column
  pa_col <- grep("(Pres|PA)", names(occ_fit), value = TRUE, ignore.case = TRUE)[1]
  x_col  <- if ("x" %in% names(occ_fit)) "x" else "lon"
  y_col  <- if ("y" %in% names(occ_fit)) "y" else "lat"
  pres <- occ_fit[occ_fit[[pa_col]] == 1, c(x_col, y_col)]
  abse <- occ_fit[occ_fit[[pa_col]] == 0, c(x_col, y_col)]
  tb_log(sprintf("presences=%d  pseudoabsences=%d", nrow(pres), nrow(abse)))
  pres_v <- terra::vect(pres, geom = c(x_col, y_col), crs = terra::crs(ens_mean_present))
  abse_v <- terra::vect(abse, geom = c(x_col, y_col), crs = terra::crs(ens_mean_present))

  thr_mean  <- .ensemble_threshold(ens_mean_present,  pres_v, abse_v)
  thr_wmean <- .ensemble_threshold(ens_wmean_present, pres_v, abse_v)
  tb_log(sprintf("MEAN   MAX_TSS threshold: %.6f", thr_mean))
  tb_log(sprintf("W_MEAN MAX_TSS threshold: %.6f", thr_wmean))

  ## Write binary present rasters
  terra::writeRaster(ens_mean_present  >= thr_mean,
                     file.path(TB_OUT_ENMTML, "Ensemble", "MEAN",   "MAX_TSS", paste0(SP, ".tif")),
                     overwrite = TRUE)
  terra::writeRaster(ens_wmean_present >= thr_wmean,
                     file.path(TB_OUT_ENMTML, "Ensemble", "W_MEAN", "MAX_TSS", paste0(SP, ".tif")),
                     overwrite = TRUE)

  ## Write binary future rasters using same thresholds
  for (scn in scenarios) {
    scn_name <- basename(scn)
    f_m <- file.path(scn, "Ensemble", "MEAN",   paste0(SP, ".tif"))
    f_w <- file.path(scn, "Ensemble", "W_MEAN", paste0(SP, ".tif"))
    if (!file.exists(f_m) || !file.exists(f_w)) next
    r_m <- terra::rast(f_m); r_w <- terra::rast(f_w)
    terra::writeRaster(r_m >= thr_mean,
                       file.path(scn, "Ensemble", "MEAN",   "MAX_TSS", paste0(SP, ".tif")),
                       overwrite = TRUE)
    terra::writeRaster(r_w >= thr_wmean,
                       file.path(scn, "Ensemble", "W_MEAN", "MAX_TSS", paste0(SP, ".tif")),
                       overwrite = TRUE)
  }

  thr_table <- data.frame(
    Sp        = SP,
    Ensemble  = c("MEAN", "W_MEAN"),
    THR       = "MAX_TSS",
    THR_VALUE = c(thr_mean, thr_wmean),
    stringsAsFactors = FALSE
  )
  write.table(thr_table, file.path(TB_OUT_ENMTML, "Thresholds_Ensemble.txt"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  tb_log(sprintf("wrote: %s/Thresholds_Ensemble.txt", TB_OUT_ENMTML))
}

## ============================================================================
## 6. POST-CHECK
## ============================================================================
tb_log_section("6. POST-CHECK")
list_outputs <- function(dir, label) {
  if (!dir.exists(dir)) { tb_log(sprintf("MISSING %s (%s)", label, dir), "WARN"); return() }
  n <- length(list.files(dir, recursive = TRUE))
  tb_log(sprintf("%-40s n_files=%d", label, n))
}
list_outputs(file.path(TB_OUT_ENMTML, "Ensemble"), "Ensemble/")
for (scn in head(scenarios, 3)) {
  list_outputs(file.path(scn, "Ensemble"), sprintf("Projection/%s/Ensemble/", basename(scn)))
}

tb_log_session()
tb_log("06z_ensemble_manual DONE")
