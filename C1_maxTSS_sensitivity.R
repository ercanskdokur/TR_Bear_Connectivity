#!/usr/bin/env Rscript
# ============================================================================
# C1_maxTSS_sensitivity.R
# TR_Bear_Connectivity: sensitivity of habitat delineation to the ensemble
#   suitability threshold (MAX_TSS), evaluated on BOTH ensemble variants.
#
# Thresholds are re-derived here with the same presence/background MAX_TSS
# procedure used to threshold the ensembles upstream, rather than hard-coded.
#
# Outputs:
#   C1_maxTSS_sensitivity.csv   — one row per (ensemble, threshold)
#   C1_maxTSS_n_patches.png, C1_maxTSS_habitat_area.png,
#   C1_maxTSS_largest_patch.png
#
# Earlier fixes retained: patch delineation uses terra::patches()+terra::freq()
#   (the idiom used by 35/36/37/C3) rather than rasterToPolygons(dissolve=FALSE),
#   which vectorised every cell separately and left the result table empty.
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse); library(terra); library(ggplot2)
})

.tb_find_paths_R <- function() {
  a <- commandArgs(trailingOnly = FALSE); f <- a[grepl("--file=", a)]
  if (length(f)) { d <- dirname(normalizePath(sub("--file=", "", f[1]), mustWork = FALSE))
    if (file.exists(file.path(d, "00_paths.R"))) return(d) }
  env_dir <- Sys.getenv("TB_SCRIPTS", unset = "")
  if (nzchar(env_dir) && file.exists(file.path(env_dir, "00_paths.R"))) return(env_dir)
  if (file.exists("00_paths.R")) return(getwd()); stop("Cannot find 00_paths.R")
}
setwd(.tb_find_paths_R()); source("00_paths.R"); source("00_helpers.R")

cat("\n=== STRAND C1: MAX_TSS Sensitivity (both ensemble variants) ===\n")

SP       <- "Ursus_arctos"
min_size <- TB_PATCH_MIN_KM2
N_BG     <- 10000L

# ---- 1) LOAD BOTH CONTINUOUS ENSEMBLE SURFACES -----------------------------
f7 <- file.path(TB_OUT_ENMTML, "Ensemble", "W_MEAN", paste0(SP, ".tif"))
f8 <- file.path(TB_OUT_ENMTML, "Ensemble", "W_MEAN", paste0(SP, "_orig8alg.tif"))

if (!file.exists(f7)) stop("Ensemble not found; run the SDM first (06_enmtml_run.R)")

surfaces <- list()
surfaces[["7-algorithm (MAH excluded)"]] <- f7
if (file.exists(f8)) {
  surfaces[["8-algorithm (all)"]] <- f8
} else {
  cat("NOTE: <sp>_orig8alg.tif absent - 24_ensemble_no_mah.R has not run,\n")
  cat("      so only the current ensemble is swept.\n")
}

# ---- 2) PRESENCE + BACKGROUND FOR THE MAX_TSS CRITERION --------------------
set.seed(TB_SEED)
mask_r <- terra::rast(file.path(TB_OUT_ENMTML, "Extent_Masks", paste0(SP, ".tif")))
occ    <- read.table(file.path(TB_OUT_ENMTML, "Occurrences_Cleaned.txt"),
                     header = TRUE, sep = "\t", stringsAsFactors = FALSE)
pres_mat <- as.matrix(occ[, c("x", "y")])
bg_xy <- terra::spatSample(mask_r, size = N_BG, method = "random",
                           na.rm = TRUE, xy = TRUE, values = FALSE)
bg_mat <- as.matrix(bg_xy[, c("x", "y")])
cat(sprintf("presences = %d | background = %d\n", nrow(pres_mat), nrow(bg_mat)))

.max_tss_threshold <- function(r, pres_mat, bg_mat) {
  vp <- terra::extract(r, pres_mat); va <- terra::extract(r, bg_mat)
  if (is.data.frame(vp)) vp <- vp[, 1]
  if (is.data.frame(va)) va <- va[, 1]
  vp <- vp[!is.na(vp)]; va <- va[!is.na(va)]
  if (length(vp) < 5 || length(va) < 5) return(NA_real_)
  cand <- unique(c(vp, va))
  if (length(cand) > 500)
    cand <- quantile(cand, probs = seq(0.005, 0.995, length.out = 500), names = FALSE)
  cand <- sort(unique(cand))
  ts <- vapply(cand, function(thr) mean(vp >= thr) + mean(va < thr) - 1, numeric(1))
  cand[which.max(ts)]
}

optima <- vapply(surfaces, function(p) .max_tss_threshold(terra::rast(p), pres_mat, bg_mat),
                 numeric(1))
for (nm in names(optima)) cat(sprintf("own MAX_TSS optimum | %-28s = %.6f\n", nm, optima[[nm]]))

# ---- 3) COMMON THRESHOLD GRID ----------------------------------------------
# Spans every optimum found, padded by 0.03 either side, so each surface is
# evaluated both at its own optimum and at the other's.
grid_lo  <- floor((min(optima, na.rm = TRUE) - 0.03) * 100) / 100
grid_hi  <- ceiling((max(optima, na.rm = TRUE) + 0.03) * 100) / 100
tss_grid <- sort(unique(c(seq(grid_lo, grid_hi, by = 0.01), unname(optima))))
cat(sprintf("\nthreshold grid: %.2f to %.2f by 0.01, plus each own optimum (%d values)\n\n",
            grid_lo, grid_hi, length(tss_grid)))

# ---- 4) SWEEP ---------------------------------------------------------------
rows <- list()
for (nm in names(surfaces)) {
  suitability <- terra::rast(surfaces[[nm]])
  if (terra::nlyr(suitability) > 1) suitability <- suitability[[1]]
  cell_km2 <- prod(terra::res(suitability)) / 1e6
  own_opt  <- optima[[nm]]
  cat(sprintf("--- %s ---\n", nm))

  for (tss in tss_grid) {
    binary <- terra::ifel(suitability >= tss, 1, NA)
    pat <- terra::patches(binary, directions = 8, zeroAsNA = TRUE)
    fr  <- terra::freq(pat) |> as.data.frame() |>
      dplyr::rename(patch_id = value, n_cells = count) |>
      dplyr::mutate(area_km2 = n_cells * cell_km2)

    fr_filt <- fr |> dplyr::filter(area_km2 >= min_size)
    if (!nrow(fr_filt)) {
      cat(sprintf("  thr %.4f: no patch passes the %d km2 filter\n", tss, min_size)); next
    }
    rows[[paste(nm, tss)]] <- tibble::tibble(
      ensemble            = nm,
      threshold           = tss,
      is_own_optimum      = isTRUE(all.equal(tss, own_opt)),
      offset_from_own_opt = tss - own_opt,
      n_patches           = nrow(fr),
      n_patches_filtered  = nrow(fr_filt),
      total_area_km2      = sum(fr$area_km2),
      filtered_area_km2   = sum(fr_filt$area_km2),
      largest_patch_km2   = max(fr_filt$area_km2),
      mean_patch_km2      = mean(fr_filt$area_km2))
    cat(sprintf("  thr %.4f%s: %d patches, %.0f km2 (>=%d km2: %d patches, %.0f km2)\n",
                tss, if (isTRUE(all.equal(tss, own_opt))) " *own optimum*" else "",
                nrow(fr), sum(fr$area_km2), min_size, nrow(fr_filt), sum(fr_filt$area_km2)))
  }
}
sensitivity_results <- dplyr::bind_rows(rows)
if (!nrow(sensitivity_results))
  stop("No threshold in the tested range produced any source patch -- check TB_PATCH_MIN_KM2")

# ---- 5) SUMMARY -------------------------------------------------------------
cat("\n=== VALUES AT EACH ENSEMBLE'S OWN OPTIMUM ===\n")
print(sensitivity_results |>
        dplyr::filter(is_own_optimum) |>
        dplyr::select(ensemble, threshold, n_patches_filtered,
                      filtered_area_km2, largest_patch_km2))

cat("\n=== RANGE ACROSS THE WHOLE GRID, PER ENSEMBLE ===\n")
print(sensitivity_results |>
        dplyr::group_by(ensemble) |>
        dplyr::summarise(
          thr_min      = min(threshold), thr_max = max(threshold),
          cores_min    = min(n_patches_filtered), cores_max = max(n_patches_filtered),
          area_min_km2 = min(filtered_area_km2), area_max_km2 = max(filtered_area_km2),
          cv_cores     = sd(n_patches_filtered) / mean(n_patches_filtered),
          .groups = "drop"))

# ---- 6) FIGURES -------------------------------------------------------------
opt_df <- tibble::tibble(ensemble = names(optima), opt = unname(optima))
pal <- c("7-algorithm (MAH excluded)" = "#2E5F8A", "8-algorithm (all)" = "#C0562A")

.mk <- function(yvar, ylab, ttl) {
  ggplot(sensitivity_results, aes(x = threshold, y = .data[[yvar]], colour = ensemble)) +
    geom_vline(data = opt_df, aes(xintercept = opt, colour = ensemble),
               linetype = "dashed", show.legend = FALSE) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.4) +
    geom_point(data = dplyr::filter(sensitivity_results, is_own_optimum),
               size = 4.2, shape = 21, fill = "white", stroke = 1.3) +
    scale_colour_manual(values = pal, name = "Ensemble") +
    labs(title = ttl, subtitle = "Dashed line and ringed point mark each ensemble's own MAX_TSS optimum",
         x = "Suitability threshold", y = ylab) +
    theme_minimal(base_size = 12) + theme(legend.position = "bottom")
}

ggsave(file.path(TB_OUT_FIGURES, "C1_maxTSS_n_patches.png"),
       .mk("n_patches_filtered", sprintf("Source patches (>= %d km2)", min_size),
           "Threshold sensitivity: number of source patches"), width = 9, height = 6.5, dpi = 300)
ggsave(file.path(TB_OUT_FIGURES, "C1_maxTSS_habitat_area.png"),
       .mk("filtered_area_km2", "Source-patch area (km2)",
           "Threshold sensitivity: total source-patch area"), width = 9, height = 6.5, dpi = 300)
ggsave(file.path(TB_OUT_FIGURES, "C1_maxTSS_largest_patch.png"),
       .mk("largest_patch_km2", "Largest patch (km2)",
           "Threshold sensitivity: largest source patch"), width = 9, height = 6.5, dpi = 300)

# ---- 7) SAVE ----------------------------------------------------------------
write_csv(sensitivity_results, file.path(TB_OUT_TABLES, "C1_maxTSS_sensitivity.csv"))
cat("\n✓ Saved: C1_maxTSS_sensitivity.csv\n")
cat("✓ Figures: C1_maxTSS_*.png\n")

cv <- sensitivity_results |>
  dplyr::group_by(ensemble) |>
  dplyr::summarise(cv = sd(n_patches_filtered) / mean(n_patches_filtered), .groups = "drop")
cat("\nRobustness of patch count to the threshold (CV within each ensemble):\n")
for (i in seq_len(nrow(cv)))
  cat(sprintf("  %-28s CV = %.3f -> %s\n", cv$ensemble[i], cv$cv[i],
              if (cv$cv[i] < 0.10) "stable" else if (cv$cv[i] < 0.25) "moderate" else "sensitive"))

cat("\n=== STRAND C1 COMPLETE ===\n")
