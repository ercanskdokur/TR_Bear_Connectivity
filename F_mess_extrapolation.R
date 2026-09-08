#!/usr/bin/env Rscript
# ============================================================================
# F_mess_extrapolation.R  (rewritten 2026-08-19)
# TR_Bear_Connectivity: Non-analog (extrapolated) climate area per future scenario
# (eval item 9 / S4.4 -- niche truncation concern)
#
# This computes it directly rather than depending on ENMTML's internal MOP:
# Multivariate Environmental Similarity Surface (MESS), via dismo::mess()
# (terra 1.7.71 on this cluster predates terra::mess(), so each target stack
# is written to a temp GeoTIFF and re-read as a raster::stack() for the dismo
# call). Reference envelope = present-day environmental values at
# occupied/accessible land cells (the same space the model was trained on).
# Target = each of the 18 future GCM x SSP x period predictor stacks, using
# the SAME retained predictor set as the main SDM (Table S1, 17 variables;
# hardcoded below to match, since ENMTML's internal VIF-reduced set is not
# re-exported as a plain list anywhere in the pipeline outputs).
#
# MESS < 0 flags a cell whose environment falls outside the training range on
# at least one variable (non-analog / true extrapolation). Per-GCM non-analog
# fraction is computed first, then averaged across
# the 3 GCMs within each of the 6 SSP x period scenarios -- consistent with
# how every other climate summary in this pipeline is GCM-averaged as a
# derived quantity, not by averaging raw layers before the metric is computed.
#
# Outputs:
#   tables/F_mess_extrapolation.csv        scenario x gcm -> pct_nonanalog_national, pct_nonanalog_suitable
#   tables/F_mess_extrapolation_summary.csv scenario (GCM-averaged) -> mean/min/max pct_nonanalog
#   figures/F_mess_extrapolation/fig_F_mess_maps.png      one present-vs-worst-case MESS map per SSP (near/far x worst SSP)
#   figures/F_mess_extrapolation/fig_F_nonanalog_bar.png  GCM-mean +/- range, by scenario
# ============================================================================

suppressPackageStartupMessages({ library(terra); library(raster); library(dismo); library(dplyr); library(ggplot2) })
select <- dplyr::select  # guard: raster masks dplyr::select

.tb_find_paths_R <- function() {
  a <- commandArgs(trailingOnly = FALSE); f <- a[grepl("--file=", a)]
  if (length(f)) { d <- dirname(normalizePath(sub("--file=", "", f[1]), mustWork = FALSE))
    if (file.exists(file.path(d, "00_paths.R"))) return(d) }
  env_dir <- Sys.getenv("TB_SCRIPTS", unset = "")
  if (nzchar(env_dir) && file.exists(file.path(env_dir, "00_paths.R"))) return(env_dir)
  if (file.exists("00_paths.R")) return(getwd()); stop("Cannot find 00_paths.R")
}
setwd(.tb_find_paths_R()); source("00_paths.R"); source("00_helpers.R")
tb_log_init("F_mess_extrapolation")

FIG_SUBDIR <- "F_mess_extrapolation"

## Table S1 retained predictor set (VIF-reduced, used in the main SDM)
RETAINED_VARS <- c("Slope", "Bio09", "d2Forest", "GHMTS", "Bio08", "Bio03",
                    "Elevation", "Bio14", "Bio15", "Bio04", "PopDen",
                    "d2ArtificialSurfaces", "d2Roads", "Bio19", "d2Water",
                    "d2Agricultural", "Aspect")

TB_GCMS <- c("GFDL_ESM4", "IPSL_CM6A_LR", "MPI_ESM1_2_HR")

present_dir <- file.path(TB_PRED_ENMTML_DIR, "present")
future_root <- file.path(TB_PRED_ENMTML_DIR, "future")

load_stack <- function(dir) {
  files <- file.path(dir, paste0(RETAINED_VARS, ".tif"))
  miss <- !file.exists(files)
  if (any(miss)) stop("Missing predictor(s) in ", dir, ": ", paste(RETAINED_VARS[miss], collapse=", "))
  r <- terra::rast(files); names(r) <- RETAINED_VARS; r
}

tr_mask <- terra::vect(file.path(TB_DATA_ROOT, "TR_mask", "TR_mask.shp")) |>
  terra::project(TB_CRS_PROJ)

tb_log("Loading present (reference) predictor stack...")
present_stack <- load_stack(present_dir)
present_stack <- terra::mask(present_stack, tr_mask)
ref_vals <- terra::values(present_stack, na.rm = TRUE)
ref_vals <- as.data.frame(ref_vals)
tb_log(sprintf("Reference envelope: %d cells x %d variables", nrow(ref_vals), ncol(ref_vals)))

bp <- file.path(TB_OUT_HS_BINARY, "present_wmean.tif")
present_suit_mask <- if (file.exists(bp)) terra::rast(bp) else NULL

## ---- per-GCM MESS for all 18 future projections ----------------------------
future_dirs <- list.dirs(future_root, recursive = FALSE, full.names = TRUE)
tb_log(sprintf("Found %d future predictor folders", length(future_dirs)))

rows <- list()
mess_cache <- list()
for (fd in future_dirs) {
  scen_full <- basename(fd)   # e.g. 2041_2070_GFDL_ESM4_ssp126
  ## parse period_gcm_ssp -> period, gcm, ssp
  gcm_hit <- TB_GCMS[sapply(TB_GCMS, function(g) grepl(g, scen_full, fixed = TRUE))]
  if (length(gcm_hit) != 1) { tb_log(sprintf("  [skip] cannot parse GCM from %s", scen_full)); next }
  ssp_hit <- TB_SSPS[sapply(TB_SSPS, function(s) grepl(s, scen_full, fixed = TRUE))]
  period_hit <- TB_PERIODS[sapply(TB_PERIODS, function(p) grepl(p, scen_full, fixed = TRUE))]
  if (length(ssp_hit) != 1 || length(period_hit) != 1) { tb_log(sprintf("  [skip] cannot parse period/ssp from %s", scen_full)); next }
  scenario <- paste(period_hit, ssp_hit, sep = "_")

  fut_stack <- load_stack(fd)
  fut_stack <- terra::mask(terra::resample(fut_stack, present_stack), tr_mask)
  tmp_tif <- tempfile(fileext = ".tif")
  terra::writeRaster(fut_stack, tmp_tif, overwrite = TRUE)
  fut_rs <- raster::stack(tmp_tif); names(fut_rs) <- RETAINED_VARS
  m_raster <- dismo::mess(fut_rs, ref_vals, full = FALSE)
  m <- terra::rast(m_raster)
  file.remove(tmp_tif)
  mess_cache[[scen_full]] <- m

  vals <- terra::values(m, na.rm = TRUE)
  n_national <- length(vals); n_nonanalog <- sum(vals < 0, na.rm = TRUE)
  pct_national <- 100 * n_nonanalog / n_national

  pct_suitable <- NA_real_
  if (!is.null(present_suit_mask)) {
    sm <- terra::resample(present_suit_mask, m, method = "near")
    vs <- terra::values(m)[terra::values(sm) == 1]; vs <- vs[!is.na(vs)]
    if (length(vs)) pct_suitable <- 100 * sum(vs < 0) / length(vs)
  }

  rows[[scen_full]] <- data.frame(scenario = scenario, gcm = gcm_hit,
    pct_nonanalog_national = pct_national, pct_nonanalog_suitable = pct_suitable)
  tb_log(sprintf("  %-40s national=%.2f%%  suitable=%s", scen_full, pct_national,
                 if (is.na(pct_suitable)) "NA" else sprintf("%.2f%%", pct_suitable)))
}
out <- do.call(rbind, rows)
tb_save_table(out, "F_mess_extrapolation")

## present self-check (should be ~0)
tmp_tif0 <- tempfile(fileext = ".tif")
terra::writeRaster(present_stack, tmp_tif0, overwrite = TRUE)
present_rs <- raster::stack(tmp_tif0); names(present_rs) <- RETAINED_VARS
present_mess_r <- dismo::mess(present_rs, ref_vals, full = FALSE)
file.remove(tmp_tif0)
pv <- raster::values(present_mess_r); pv <- pv[!is.na(pv)]
pct_present <- 100 * sum(pv < 0) / length(pv)
tb_log(sprintf("Present self-comparison non-analog: %.3f%% (sanity check, should be ~0)", pct_present))

## ---- GCM-averaged summary per scenario --------------------------------------
scen_label <- c(present = "Present",
  "2041_2070_ssp126" = "2070s SSP126", "2041_2070_ssp370" = "2070s SSP370",
  "2041_2070_ssp585" = "2070s SSP585", "2071_2100_ssp126" = "2100s SSP126",
  "2071_2100_ssp370" = "2100s SSP370", "2071_2100_ssp585" = "2100s SSP585")

summary_df <- out |>
  group_by(scenario) |>
  summarise(mean_pct_nonanalog_national = mean(pct_nonanalog_national),
            min_pct_nonanalog_national  = min(pct_nonanalog_national),
            max_pct_nonanalog_national  = max(pct_nonanalog_national),
            mean_pct_nonanalog_suitable = mean(pct_nonanalog_suitable, na.rm = TRUE),
            .groups = "drop") |>
  mutate(label = scen_label[scenario])
tb_save_table(summary_df, "F_mess_extrapolation_summary")
tb_log("GCM-averaged non-analog area by scenario:")
for (i in seq_len(nrow(summary_df))) {
  r <- summary_df[i, ]
  tb_log(sprintf("  %-16s national=%.1f%% [%.1f-%.1f]  suitable-habitat=%.1f%%",
    r$label, r$mean_pct_nonanalog_national, r$min_pct_nonanalog_national,
    r$max_pct_nonanalog_national, r$mean_pct_nonanalog_suitable))
}

## ---- figures -----------------------------------------------------------------
tb_log_section("Figures")

## one representative MESS map per SSP x period (first GCM available) to keep the panel readable
show_scen <- future_dirs[sapply(names(mess_cache), function(n) grepl(TB_GCMS[1], n, fixed = TRUE))]
if (length(mess_cache)) {
  keep <- names(mess_cache)[grepl(TB_GCMS[1], names(mess_cache), fixed = TRUE)]
  df_maps <- lapply(keep, function(n) {
    d <- as.data.frame(mess_cache[[n]], xy = TRUE); names(d)[3] <- "mess"
    scen <- sub(paste0("_", TB_GCMS[1]), "", n)
    d$scenario <- scen_label[[scen]]; d
  }) |> bind_rows()
  p1 <- ggplot(df_maps, aes(x, y, fill = mess)) +
    geom_raster() + coord_equal() + facet_wrap(~ scenario) +
    scale_fill_gradient2(low = "#D55E00", mid = "#F2F2F2", high = "#0072B2",
                          midpoint = 0, name = "MESS\n(<0 = non-analog)") +
    labs(title = sprintf("MESS extrapolation (%s), by scenario", TB_GCMS[1])) +
    theme_trbear(base_size = 11) + theme(axis.text = element_blank(), axis.ticks = element_blank())
  tb_save_fig(p1, "fig_F_mess_maps", w = 12, h = 9, subdir = FIG_SUBDIR)
}

p2 <- ggplot(summary_df, aes(reorder(label, mean_pct_nonanalog_national), mean_pct_nonanalog_national)) +
  geom_col(fill = "#D55E00") +
  geom_errorbar(aes(ymin = min_pct_nonanalog_national, ymax = max_pct_nonanalog_national), width = 0.25) +
  geom_text(aes(label = sprintf("%.1f%%", mean_pct_nonanalog_national)), hjust = -0.15, size = 3.2) +
  coord_flip(clip = "off") +
  labs(title = "Non-analog (extrapolated) area by scenario, GCM-averaged",
       subtitle = "MESS < 0 = non-analog; error bars = range across 3 GCMs",
       x = NULL, y = "% of national area") +
  theme_trbear_bar(base_size = 12)
tb_save_fig(p2, "fig_F_nonanalog_bar", w = 9, h = 6, subdir = FIG_SUBDIR)

tb_log_session(); tb_log("F_mess_extrapolation DONE")
