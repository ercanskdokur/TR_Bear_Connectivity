#!/usr/bin/env Rscript
# ============================================================================
# C2_pseudoabsence_sensitivity.R  (fixed 2026-08-19)
# TR_Bear_Connectivity: Sensitivity to Pseudo-Absence Ratio -- eval S4.8/S4.9
#
# Bugs fixed on top of the earlier library(enmtml) typo:
#  1) The ENMTML() call used a fictional argument set (occurrence=, predictors=,
#     output=, PA=, PA.ratio=, ...) that doesn't match the real function
#     signature. Rewritten to match the working call in 06_enmtml_run.R
#     (pred_dir=, occ_file=, pres_abs_ratio=, ...). proj_dir is omitted here
#     (NULL) since this strand only needs present-day patch topology, not the
#     18 future projections -- this also cuts runtime substantially.
#  2) rasterToPolygons(dissolve = FALSE) vectorised one polygon per raster
#     cell rather than per patch, so no "patch" ever reached the 83 km^2
#     filter (identical bug to C1). Replaced with terra::patches()+freq(),
#     matching the rest of the pipeline.
#  3) Each re-run's own MAX_TSS threshold (from its Thresholds_Ensemble.txt)
#     is now used for binarisation instead of a hardcoded 0.50, since a
#     different PA ratio genuinely shifts the optimal threshold.
# ============================================================================

suppressPackageStartupMessages({
  library(tidyverse); library(terra); library(ENMTML); library(ggplot2)
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

cat("\n=== STRAND C2: Pseudo-Absence Ratio Sensitivity ===\n")

occ_data <- read_delim(TB_PRESENCE_TXT, delim = "\t", show_col_types = FALSE)
n_presence <- nrow(occ_data)
cat(sprintf("\u2713 Loaded %d presence points\n", n_presence))

mask_fixed_shp <- file.path(TB_DATA_ROOT, "TR_mask", "TR_mask.shp")
min_size <- TB_PATCH_MIN_KM2

## NOTE 2026-08-20: the main production run (06_enmtml_run.R) never called
## set.seed() before ENMTML::ENMTML(), so its PA sampling and BOOT partition
## draws are not reproducible -- there is no seed to "match" exactly. What we
## CAN do is seed this sweep so its own three ratios are drawn from a common,
## reproducible baseline and are fairly comparable to each other. The absolute
## area at PA=1:1 here is still expected to differ from the archived
## production figure (141,049 km2); that run-to-run stochastic spread is
## itself an honest finding, reported alongside the PA-ratio effect.
set.seed(TB_SEED)

pa_ratios <- c(1, 10, 15)
rows <- list()

for (pa_ratio in pa_ratios) {
  cat(sprintf("\n  Testing PA ratio %d:1 (%d presences x %d = %d PAs)...\n",
              pa_ratio, n_presence, pa_ratio, n_presence * pa_ratio))
  output_dir <- file.path(TB_OUT_ENMTML, sprintf("present_pa_ratio_%d", pa_ratio))
  if (dir.exists(output_dir)) unlink(output_dir, recursive = TRUE, force = TRUE)
  dir.create(output_dir, recursive = TRUE)

  {
  ENMTML::ENMTML(
    pred_dir            = TB_PRED_ENMTML_PRESENT,
    proj_dir            = NULL,
    result_dir          = output_dir,
    occ_file            = TB_OCC_ENMTML,
    sp                  = "species",
    x                   = "x",
    y                   = "y",
    min_occ             = TB_ENM_MIN_OCC,
    thin_occ            = TB_ENM_THIN_OCC,
    eval_occ            = NULL,
    colin_var           = TB_ENM_COLIN_VAR,
    imp_var             = FALSE,
    sp_accessible_area  = c(method = "MASK", filepath = mask_fixed_shp),
    pseudoabs_method    = TB_ENM_PA_METHOD,
    pres_abs_ratio      = pa_ratio,
    part                = TB_ENM_PART,
    save_part           = FALSE,
    save_final          = TB_ENM_SAVE_FINAL,
    algorithm           = TB_ENM_ALGORITHMS,
    thr                 = TB_ENM_THR,
    msdm                = TB_ENM_MSDM,
    ensemble            = TB_ENM_ENSEMBLE,
    extrapolation       = FALSE,
    cores               = TB_ENM_CORES
  )
  cat("    \u2713 ENMTML run complete\n")
  }

  ## ---- mean TSS across algorithms ----
  eval_file <- file.path(output_dir, "Evaluation_Table.txt")
  mean_tss <- NA_real_
  if (file.exists(eval_file)) {
    ev <- read_delim(eval_file, delim = "\t", show_col_types = FALSE)
    if ("TSS" %in% names(ev)) mean_tss <- mean(ev$TSS, na.rm = TRUE)
    cat(sprintf("    Mean TSS: %s\n", if (is.na(mean_tss)) "NA" else sprintf("%.3f", mean_tss)))
  }

  ## ---- this run's own MAX_TSS threshold ----
  ## Raw ENMTML output uses a transposed layout with no real header: row1 =
  ## species (repeated), row2 = ensemble type (WMEA/MEA), row3 = threshold
  ## type, row4 = value, one column per (ensemble x threshold) combination.
  thr_file <- file.path(output_dir, "Thresholds_Ensemble.txt")
  run_thr <- 0.5
  if (file.exists(thr_file)) {
    ## file's own first line is a throwaway "V1".."V12" header written by
    ## write.table(); default col_names=TRUE consumes it as the header, so
    ## the remaining tibble rows are: 1=species(repeated), 2=ensemble type,
    ## 3=threshold type, 4=value.
    th_raw <- read_delim(thr_file, delim = "\t", show_col_types = FALSE)
    ens_row <- as.character(th_raw[2, ])
    thr_row <- as.character(th_raw[3, ])
    val_row <- as.numeric(th_raw[4, ])
    hit <- which(ens_row == "WMEA" & thr_row == "MAX_TSS")
    if (length(hit) == 1) run_thr <- val_row[hit]
  }
  cat(sprintf("    Threshold (W_MEAN MAX_TSS): %.4f\n", run_thr))

  ensemble_file <- file.path(output_dir, "Ensemble", "W_MEAN", "Ursus_arctos.tif")
  if (!file.exists(ensemble_file)) {
    cat(sprintf("    WARNING: Ensemble not found at %s\n", ensemble_file))
    next
  }
  suitability <- terra::rast(ensemble_file)
  if (terra::nlyr(suitability) > 1) suitability <- suitability[[1]]
  sv <- terra::values(suitability, na.rm = TRUE)
  if (max(sv) > 1) suitability <- suitability / max(sv)
  cell_km2 <- prod(terra::res(suitability)) / 1e6

  binary <- terra::ifel(suitability > run_thr, 1, NA)
  pat <- terra::patches(binary, directions = 8, zeroAsNA = TRUE)
  fr  <- terra::freq(pat) |> as.data.frame() |>
    dplyr::rename(patch_id = value, n_cells = count) |>
    dplyr::mutate(area_km2 = n_cells * cell_km2)
  n_patches_all <- nrow(fr)
  total_area_km2 <- sum(fr$area_km2)
  fr_filt <- fr |> dplyr::filter(area_km2 >= min_size)

  if (nrow(fr_filt) > 0) {
    rows[[as.character(pa_ratio)]] <- tibble::tibble(
      PA_ratio = pa_ratio, n_presences = n_presence, n_pseudoabsences = n_presence * pa_ratio,
      mean_TSS = mean_tss, threshold_used = run_thr,
      n_patches = n_patches_all, n_patches_filtered = nrow(fr_filt),
      total_area_km2 = total_area_km2, filtered_area_km2 = sum(fr_filt$area_km2),
      largest_patch_km2 = max(fr_filt$area_km2))
    cat(sprintf("    Habitat: %.0f km2 (%d patches >=%d km2)\n",
                sum(fr_filt$area_km2), nrow(fr_filt), min_size))
  } else {
    cat(sprintf("    WARNING: no patch passed the %d km2 filter at PA ratio %d\n", min_size, pa_ratio))
  }
}

sensitivity_results <- dplyr::bind_rows(rows)
if (!nrow(sensitivity_results)) stop("No PA ratio produced any source patch -- check ENMTML output above")

cat("\n=== SENSITIVITY RESULTS ===\n")
print(sensitivity_results)

cat("\nInterpretation:\n")
cat("- If TSS and habitat area stable across PA ratios: ROBUST to PA sampling\n")
cat("- If TSS/area changes dramatically: SENSITIVE to PA sampling; use balanced or higher PA ratio\n")

p1 <- sensitivity_results %>%
  ggplot(aes(x = PA_ratio, y = mean_TSS)) +
  geom_line(linewidth = 1.2, color = "steelblue") + geom_point(size = 4, color = "steelblue") +
  labs(title = "PA Ratio Sensitivity: Model Performance (TSS)",
       x = "Pseudo-Absence : Presence Ratio", y = "Mean TSS") + theme_minimal()
ggsave(file.path(TB_OUT_FIGURES, "C2_PA_ratio_TSS.png"), p1, width = 8, height = 6)

p2 <- sensitivity_results %>%
  ggplot(aes(x = PA_ratio, y = filtered_area_km2)) +
  geom_line(linewidth = 1.2, color = "darkgreen") + geom_point(size = 4, color = "darkgreen") +
  labs(title = "PA Ratio Sensitivity: Habitat Area",
       x = "Pseudo-Absence : Presence Ratio", y = "Total Habitat Area (km2)") + theme_minimal()
ggsave(file.path(TB_OUT_FIGURES, "C2_PA_ratio_habitat_area.png"), p2, width = 8, height = 6)

p3 <- sensitivity_results %>%
  ggplot(aes(x = PA_ratio, y = n_patches_filtered)) +
  geom_line(linewidth = 1.2, color = "darkorange") + geom_point(size = 4, color = "darkorange") +
  labs(title = "PA Ratio Sensitivity: Number of Patches",
       x = "Pseudo-Absence : Presence Ratio", y = sprintf("Number of Patches (>=%d km2)", min_size)) +
  theme_minimal()
ggsave(file.path(TB_OUT_FIGURES, "C2_PA_ratio_n_patches.png"), p3, width = 8, height = 6)

write_csv(sensitivity_results, file.path(TB_OUT_TABLES, "C2_PA_ratio_sensitivity.csv"))
cat("\n\u2713 Saved: C2_PA_ratio_sensitivity.csv\n")
cat("\n=== STRAND C2 COMPLETE ===\n")
