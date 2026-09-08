#!/usr/bin/env Rscript
# ============================================================================
# C3_c_sensitivity_full.R
# TR_Bear_Connectivity: FULL c-sensitivity across ALL climate scenarios (9 / eval S4.8)
#
#
# Design: identical resistance/cost-distance/PC/IIC machinery as scripts 35 and
# 36. Patches per scenario come from the existing binary MAX_TSS maps (fixed,
# independent of c). Resistance is regenerated from the scenario's CONTINUOUS
# suitability raster (hs_present/wmean.tif; hs_gcm_avg/{scenario}.tif) for each
# c, because c changes the suitability->resistance transform. R_eff is fixed at
# its PRESENT, c=4 value (Reff = 5.99, from 35_costdist_calibration.csv) for
# ALL scenarios and ALL c, exactly as in the main manuscript pipeline, so the
# only thing that varies across the sweep is the resistance surface itself.
#
# Outputs:
#   tables/C3_c_sensitivity_full.csv   scenario x c x d_km -> PC, IIC, ECA, n50_PC
#   tables/C3_manuscript_summary.csv   scenario x c (d=100 focal) -> PC_retention_pct,
#                                        ECA_retention_pct, artifact_ratio (=ECA_ret/area_ret)
#   figures/C3_c_sensitivity/fig_C3_PC_loss.png
#   figures/C3_c_sensitivity/fig_C3_PC_vs_ECA.png
#   figures/C3_c_sensitivity/fig_C3_heatmap.png
# ============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("igraph", quietly = TRUE))
    install.packages("igraph", repos = "https://cloud.r-project.org")
  library(terra); library(sf); library(ggplot2); library(dplyr); library(tidyr); library(igraph)
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
tb_log_init("C3_c_sensitivity_full")

FIG_SUBDIR     <- "C3_c_sensitivity"
C_VALUES       <- c(2, 4, 6, 8)
DIST_KM        <- c(50, 100, 200)
FOCAL_D        <- 100
PROB_AT_THRESH <- 0.5

scenarios <- c("present",
               sprintf("%s_%s", rep(TB_PERIODS, each = length(TB_SSPS)),
                                 rep(TB_SSPS,    times = length(TB_PERIODS))))

.bin_path <- function(s) {
  if (s == "present") file.path(TB_OUT_HS_BINARY, "present_wmean.tif")
  else                file.path(TB_OUT_HS_BINARY, sprintf("future_%s.tif", s))
}
.hs_path <- function(s) {
  if (s == "present") file.path(TB_OUT_HS_PRESENT, "wmean.tif")
  else                file.path(TB_OUT_HS_AVG, sprintf("%s.tif", s))
}
.hs_to_R <- function(h, c, rmin = TB_RESIST_MIN, rmax = TB_RESIST_MAX)
  rmax - (rmax - rmin) * (1 - exp(-c * h)) / (1 - exp(-c))

tr_mask_sf <- sf::st_read(file.path(TB_DATA_ROOT, "TR_mask", "TR_mask.shp"),
                          quiet = TRUE) |> sf::st_transform(TB_CRS_PROJ)
AL_km2 <- as.numeric(sum(sf::st_area(tr_mask_sf))) / 1e6
tb_log(sprintf("A_L = %.0f km2", AL_km2))

## R_eff fixed at present, c=4 (manuscript convention: script 35 calibration)
cal <- tryCatch(read.csv(file.path(TB_OUT_TABLES, "35_costdist_calibration.csv")),
                error = function(e) NULL)
if (is.null(cal)) stop("35_costdist_calibration.csv not found -- run script 35 first")
R_eff_fixed <- cal$R_eff[1]
tb_log(sprintf("R_eff (fixed, present c=4) = %.3f", R_eff_fixed))

## ---- patch delineation per scenario (binary map; c-independent) -----------
delineate <- function(s, min_km2 = TB_PATCH_MIN_KM2) {
  pp <- .bin_path(s); if (!file.exists(pp)) return(NULL)
  r <- terra::rast(pp); names(r) <- "suit"
  cell_km2 <- prod(terra::res(r)) / 1e6
  pat <- terra::patches(r, directions = 8, zeroAsNA = TRUE)
  fr  <- terra::freq(pat) |> as.data.frame() |>
    dplyr::rename(patch_id = value, n_cells = count) |>
    dplyr::mutate(area_km2 = n_cells * cell_km2) |>
    dplyr::filter(!is.na(patch_id), area_km2 >= min_km2)
  if (!nrow(fr)) return(NULL)
  keep <- terra::ifel(pat %in% fr$patch_id, pat, NA); names(keep) <- "patch"
  cen <- lapply(fr$patch_id, function(i) {
    crd <- terra::crds(terra::as.points(terra::ifel(keep == i, 1, NA)))
    data.frame(patch_id = i, x = mean(crd[, 1]), y = mean(crd[, 2]))
  }) |> do.call(what = rbind)
  list(ids = fr$patch_id, area = fr$area_km2, cen = cen, template = r)
}

cost_matrix <- function(res_r, cen_xy) {
  n <- nrow(cen_xy); cells <- terra::cellFromXY(res_r, cen_xy)
  D <- matrix(NA_real_, n, n)
  for (k in seq_len(n)) {
    orig <- res_r; orig[cells[k]] <- 0
    cd <- terra::costDist(orig, target = 0, scale = 1, maxiter = 3000)
    D[k, ] <- terra::extract(cd, cen_xy)[, 1]
  }
  D <- (D + t(D)) / 2
  fin <- D[is.finite(D)]; big <- if (length(fin)) max(fin) * 10 else 1e9
  D[!is.finite(D)] <- big; diag(D) <- 0; D / 1000
}
PC_of <- function(D, area, keepv, d_eff) {
  ar <- area[keepv]; if (length(ar) < 2) return(if (length(ar)) ar^2 / AL_km2^2 else 0)
  alpha <- -log(PROB_AT_THRESH) / d_eff
  W <- -log(pmax(exp(-alpha * D[keepv, keepv]), .Machine$double.eps)); diag(W) <- 0
  g <- igraph::graph_from_adjacency_matrix(W, mode = "undirected", weighted = TRUE)
  sum(outer(ar, ar) * exp(-igraph::distances(g, weights = igraph::E(g)$weight))) / AL_km2^2
}
IIC_of <- function(D, area, d_eff) {
  Adj <- (D <= d_eff) & (D > 0)
  g <- igraph::graph_from_adjacency_matrix(Adj, mode = "undirected")
  sum(outer(area, area) / (1 + igraph::distances(g)), na.rm = TRUE) / AL_km2^2
}
dPC_vec <- function(D, area, d_eff) {
  n <- length(area); pc0 <- PC_of(D, area, seq_len(n), d_eff)
  vapply(seq_len(n), function(k)
    if (pc0 > 0) 100 * (pc0 - PC_of(D, area, setdiff(seq_len(n), k), d_eff)) / pc0 else 0, numeric(1))
}
n50_targeted <- function(D, area, d_eff) {
  n <- length(area); pc0 <- PC_of(D, area, seq_len(n), d_eff)
  ord <- order(dPC_vec(D, area, d_eff), decreasing = TRUE)
  remaining <- seq_len(n)
  for (i in seq_len(n - 1L)) {
    remaining <- setdiff(remaining, ord[i])
    if (PC_of(D, area, remaining, d_eff) < 0.5 * pc0) return(i)
  }
  NA_integer_
}

## ---- main sweep: scenario x c ----------------------------------------------
rows <- list()
for (s in scenarios) {
  tb_log_section(sprintf("scenario = %s", s))
  pd <- delineate(s)
  if (is.null(pd)) { tb_log(sprintf("  [skip] no patches for %s", s)); next }
  hs_path <- .hs_path(s)
  if (!file.exists(hs_path)) { tb_log(sprintf("  [skip] HS raster missing: %s", hs_path)); next }
  hs <- terra::rast(hs_path); if (terra::nlyr(hs) > 1) hs <- hs[[1]]
  hs <- terra::resample(hs, pd$template, method = "bilinear")
  cen_xy <- as.matrix(pd$cen[, c("x", "y")])
  tb_log(sprintf("  n_patches = %d, hab_km2 = %.0f", length(pd$ids), sum(pd$area)))

  for (cc in C_VALUES) {
    res_c <- terra::app(hs, function(h) .hs_to_R(h, cc)); names(res_c) <- "R"
    D <- cost_matrix(res_c, cen_xy)
    for (d in DIST_KM) {
      d_eff <- d * R_eff_fixed
      pc  <- PC_of(D, pd$area, seq_len(length(pd$ids)), d_eff)
      iic <- IIC_of(D, pd$area, d_eff)
      eca <- sqrt(pc) * AL_km2
      n50 <- if (d == FOCAL_D) n50_targeted(D, pd$area, d_eff) else NA_integer_
      rows[[paste(s, cc, d)]] <- data.frame(
        scenario = s, c = cc, d_km = d, n_patches = length(pd$ids),
        hab_km2 = sum(pd$area), PC = pc, IIC = iic, ECA_km2 = eca, n50_PC = n50)
    }
    tb_log(sprintf("    c=%g done", cc))
  }
}
out <- do.call(rbind, rows)
tb_save_table(out, "C3_c_sensitivity_full")

## ---- manuscript summary: retention relative to present, at focal d --------
foc <- out |> dplyr::filter(d_km == FOCAL_D)
pres <- foc |> dplyr::filter(scenario == "present") |>
  dplyr::select(c, PC_present = PC, ECA_present = ECA_km2, hab_present = hab_km2)
summary_df <- foc |>
  dplyr::left_join(pres, by = "c") |>
  dplyr::mutate(
    area_retention_pct = 100 * hab_km2 / hab_present,
    PC_retention_pct    = 100 * PC / PC_present,
    ECA_retention_pct   = 100 * ECA_km2 / ECA_present,
    artifact_ratio       = ECA_retention_pct / area_retention_pct) |>
  dplyr::select(scenario, c, n_patches, hab_km2, area_retention_pct,
                PC, PC_retention_pct, ECA_km2, ECA_retention_pct, artifact_ratio, n50_PC)
tb_save_table(summary_df, "C3_manuscript_summary")

## ---- figures ----------------------------------------------------------------
tb_log_section("Figures")
scen_label <- c(present = "Present",
  "2041_2070_ssp126" = "2070s SSP126", "2041_2070_ssp370" = "2070s SSP370",
  "2041_2070_ssp585" = "2070s SSP585", "2071_2100_ssp126" = "2100s SSP126",
  "2071_2100_ssp370" = "2100s SSP370", "2071_2100_ssp585" = "2100s SSP585")
summary_df$lbl <- factor(scen_label[summary_df$scenario], levels = scen_label[scenarios])

## Two-hue structural palette: period (2070s = blue family, 2100s = orange
## family) x severity (light -> dark = SSP126 -> SSP370 -> SSP585). Pastel and
## colourblind-safe (blue-orange axis is preserved under all CVD types).
scen_pal <- c(
  "2070s SSP126" = "#A8D0E6", "2070s SSP370" = "#5B9BD5", "2070s SSP585" = "#2E5F8A",
  "2100s SSP126" = "#F5C89A", "2100s SSP370" = "#E8935A", "2100s SSP585" = "#C0562A"
)

p1 <- ggplot(summary_df |> dplyr::filter(scenario != "present"),
             aes(factor(c), PC_retention_pct, fill = lbl)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = scen_pal) +
  labs(title = "PC retention (% of present) by resistance shape constant c",
       x = "c", y = "PC retention (%)", fill = "Scenario") +
  theme_trbear_bar(base_size = 12)
tb_save_fig(p1, "fig_C3_PC_loss", w = 10, h = 6, subdir = FIG_SUBDIR)

p2 <- ggplot(summary_df |> dplyr::filter(scenario != "present"),
             aes(area_retention_pct, ECA_retention_pct, color = lbl, shape = factor(c))) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 3) +
  scale_color_manual(values = scen_pal) +
  labs(title = "ECA retention vs. area retention (1:1 line = no artefact)",
       x = "Area retention (%)", y = "ECA retention (%)",
       color = "Scenario", shape = "c") +
  theme_trbear_bar(base_size = 12)
tb_save_fig(p2, "fig_C3_PC_vs_ECA", w = 8, h = 7, subdir = FIG_SUBDIR)

p3 <- ggplot(summary_df |> dplyr::filter(scenario != "present"),
             aes(factor(c), lbl, fill = artifact_ratio)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", artifact_ratio)), size = 3.2) +
  scale_fill_gradient2(midpoint = 1, low = "#D55E00", mid = "white", high = "#0072B2",
                        name = "ECA_ret /\narea_ret") +
  labs(title = "Artefact ratio (ECA retention / area retention) across c and scenario",
       x = "c", y = NULL) +
  theme_trbear_bar(base_size = 12)
tb_save_fig(p3, "fig_C3_heatmap", w = 9, h = 6, subdir = FIG_SUBDIR)

tb_log_session(); tb_log("C3_c_sensitivity_full DONE")
