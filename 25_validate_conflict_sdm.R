## ============================================================================
## 25_validate_conflict_sdm.R
## Project: TR_Bear_Connectivity
## Purpose: Diagnose whether the conflict SDM is mostly re-predicting bear
##   habitat (circular-predictor risk). Quantifies the spatial agreement
##   between bear HS and conflict HS rasters via:
##     (i)   Spearman rho across all cells (and within TR mask only)
##     (ii)  Pearson r across all cells
##     (iii) hexbin scatter + density overlay
##     (iv)  Difference map (conflict_HS - bear_HS) with quantile breaks
##   If rho > 0.7, the conflict SDM does not provide information beyond what
##   the bear SDM already encodes, and is treated as encounter-conditioned conflict risk.
##
## Inputs:
##   bear     : <TB_OUT_ENMTML>/Ensemble/W_MEAN/Ursus_arctos.tif
##   conflict : <TB_OUT_ENMTML_CONFLICT>/Ensemble/W_MEAN/Ursus_arctos_conflict.tif
##   TR_mask  : <TB_DATA_ROOT>/TR_mask/TR_mask.shp
##
## Outputs (tables/):
##   25_conflict_sdm_correlation.csv
## Outputs (figures/25_validate/):
##   fig25a_scatter_hexbin.png
##   fig25b_difference_map.png
## ============================================================================

suppressPackageStartupMessages({
  library(terra); library(sf); library(ggplot2); library(dplyr); library(tidyterra)
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
tb_log_init("25_validate_conflict_sdm")

bear_tif <- file.path(TB_OUT_ENMTML,          "Ensemble", "W_MEAN",
                       "Ursus_arctos.tif")
conf_tif <- file.path(TB_OUT_ENMTML_CONFLICT, "Ensemble", "W_MEAN",
                       "Ursus_arctos_conflict.tif")
stopifnot(file.exists(bear_tif), file.exists(conf_tif))

bear <- terra::rast(bear_tif)
conf <- terra::rast(conf_tif)

if (!terra::compareGeom(bear, conf, stopOnError = FALSE))
  conf <- terra::resample(conf, bear, method = "bilinear")

tr_shp <- file.path(TB_DATA_ROOT, "TR_mask", "TR_mask.shp")
tr_sf  <- sf::st_read(tr_shp, quiet = TRUE) |> sf::st_transform(terra::crs(bear))
bear_m <- terra::mask(bear, terra::vect(tr_sf))
conf_m <- terra::mask(conf, terra::vect(tr_sf))
names(bear_m) <- "bear_hs"
names(conf_m) <- "conf_hs"

vals <- as.data.frame(c(bear_m, conf_m), na.rm = TRUE)
names(vals) <- c("bear_hs", "conf_hs")
tb_log(sprintf("n cells inside TR with both rasters = %d", nrow(vals)))

## ---- Correlations ----------------------------------------------------------
tb_log_section("Correlations")
rho_s   <- cor(vals$bear_hs, vals$conf_hs, method = "spearman")
rho_p   <- cor(vals$bear_hs, vals$conf_hs, method = "pearson")
tb_log(sprintf("Spearman rho = %.4f   Pearson r = %.4f", rho_s, rho_p))

## Subset: cells where bear_hs > 0.5 (suspected presence)
sub <- vals[vals$bear_hs > 0.5, ]
rho_s_pres <- if (nrow(sub) > 50) cor(sub$bear_hs, sub$conf_hs, method = "spearman") else NA_real_
rho_p_pres <- if (nrow(sub) > 50) cor(sub$bear_hs, sub$conf_hs, method = "pearson")  else NA_real_

decision <- dplyr::case_when(
  rho_s <  0.4              ~ "DEFENSIBLE — conflict SDM provides distinct information; report as 'conflict risk'",
  rho_s <= 0.7              ~ "CAVEATED — interpret as 'encounter-conditioned conflict risk'",
  TRUE                       ~ "REDUNDANT — conflict SDM mirrors bear HS; consider dropping or using conflict-specific predictors")

cor_df <- data.frame(
  metric = c("spearman_rho_all", "pearson_r_all",
              "spearman_rho_presence", "pearson_r_presence",
              "n_cells_all", "n_cells_presence",
              "decision"),
  value  = c(round(rho_s, 4), round(rho_p, 4),
              round(rho_s_pres, 4), round(rho_p_pres, 4),
              nrow(vals), nrow(sub),
              decision))
tb_save_table(cor_df, "25_conflict_sdm_correlation")
tb_log(sprintf("DECISION: %s", decision))

## ---- Fig25a — hexbin scatter -----------------------------------------------
tb_log_section("Fig25a hexbin scatter")
samp <- vals
if (nrow(samp) > 2e5) samp <- samp[sample(seq_len(nrow(samp)), 2e5), ]

p_sc <- ggplot(samp, aes(bear_hs, conf_hs)) +
  geom_hex(bins = 60) +
  scale_fill_viridis_c(option = "viridis", trans = "log10",
                        name = "Cell count (log10)") +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 6),
              se = FALSE, color = "#9E2A2B", linewidth = 0.8) +
  geom_abline(slope = 1, intercept = 0, color = "white",
              linewidth = 0.4, linetype = 2) +
  labs(title    = "Conflict SDM vs. bear habitat suitability (present)",
       subtitle = sprintf("Spearman rho = %.3f (all cells, n = %d) | %.3f (bear HS > 0.5, n = %d).  Red = GAM fit;  dashed = 1:1.",
                          rho_s, nrow(vals), rho_s_pres, nrow(sub)),
       x = "Bear habitat suitability (W_MEAN)",
       y = "Conflict-trained SDM (W_MEAN)") +
  theme_trbear_bar(base_size = 11)
tb_save_fig(p_sc, "fig25a_scatter_hexbin", w = 9, h = 7, subdir = "25_validate")

## ---- Fig25b — difference map -----------------------------------------------
tb_log_section("Fig25b difference map")
diff_r <- conf_m - bear_m
names(diff_r) <- "diff"

p_d <- ggplot()
world_sf <- tryCatch(
  rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") |>
    sf::st_transform(terra::crs(diff_r)),
  error = function(e) NULL)
if (!is.null(world_sf))
  p_d <- p_d + geom_sf(data = world_sf, fill = "#E8E8E8",
                       color = "#7C8A93", linewidth = 0.3)
tr_bb <- sf::st_bbox(tr_sf)
pad   <- 30000
p_d <- p_d +
  tidyterra::geom_spatraster(data = diff_r, na.rm = TRUE) +
  scale_fill_gradient2(low = "#0072B2", mid = "#F2F2F2", high = "#9E2A2B",
                       midpoint = 0, limits = c(-1, 1),
                       name = "Conflict − Bear\nsuitability",
                       na.value = NA) +
  geom_sf(data = tr_sf, fill = NA, color = TB_COLOR_FRAME, linewidth = 0.5) +
  coord_sf(xlim = c(tr_bb["xmin"] - pad, tr_bb["xmax"] + pad),
           ylim = c(tr_bb["ymin"] - pad, tr_bb["ymax"] + pad),
           datum = sf::st_crs(4326), expand = FALSE) +
  labs(title    = "Where the conflict risk deviates from the bear distribution",
       subtitle = "Blue = bear-suitable but low conflict; red = high conflict above what bear suitability predicts.") +
  theme_trbear(base_size = 11)
tb_save_fig(p_d, "fig25b_difference_map", w = 11, h = 7, subdir = "25_validate")

tb_log_session()
tb_log("25_validate_conflict_sdm DONE")
