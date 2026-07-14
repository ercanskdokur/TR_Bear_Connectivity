## ============================================================================
## 10_gcm_average.R
## Project: TR_Bear_Connectivity
## Purpose: Reduce 18 future scenarios → 6 GCM-averaged scenarios
##   (2 periods × 3 SSPs), plus per-cell CV across GCMs as uncertainty.
##   Outputs HS_avg + HS_cv rasters + final 6-panel figure pair.
## ============================================================================

suppressPackageStartupMessages({
  library(terra); library(sf); library(ggplot2); library(dplyr); library(tidyr)
  library(tidyterra); library(patchwork); library(rnaturalearth); library(scales)
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
tb_log_init("10_gcm_average")

SP <- "Ursus_arctos"
FIG_SUBDIR <- "10_gcm_avg"

## ----------------------------------------------------------------------------
## 1) Combine each (period, ssp) over 3 GCMs
## ----------------------------------------------------------------------------
tb_log_section("Build GCM averages")

thr_tbl   <- read.table(file.path(TB_OUT_ENMTML, "Thresholds_Ensemble.txt"),
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE)
thr_wmean <- thr_tbl$THR_VALUE[thr_tbl$Ensemble == "W_MEAN"]
tb_log(sprintf("threshold (W_MEAN MAX_TSS) = %.4f", thr_wmean))

scen_avg <- expand.grid(period = TB_PERIODS, ssp = TB_SSPS,
                        stringsAsFactors = FALSE)
scen_avg <- scen_avg[order(scen_avg$period, scen_avg$ssp), ]
rownames(scen_avg) <- NULL
scen_avg$scen_key <- sprintf("%s_%s", scen_avg$period, scen_avg$ssp)

avg_rasters <- vector("list", nrow(scen_avg))
cv_rasters  <- vector("list", nrow(scen_avg))
summary_rows <- vector("list", nrow(scen_avg))

cell_area_km2 <- NULL

tb_tic()
for (i in seq_len(nrow(scen_avg))) {
  p <- scen_avg$period[i]; s <- scen_avg$ssp[i]
  paths <- file.path(TB_OUT_HS_FUTURE,
                     sprintf("%s_%s_%s.tif", p, TB_GCMS, s))
  ok <- file.exists(paths)
  if (!all(ok)) {
    tb_log(sprintf("missing GCM(s) for %s_%s: %s",
                   p, s, paste(TB_GCMS[!ok], collapse = ",")), "WARN")
    if (!any(ok)) next
  }
  stk <- terra::rast(paths[ok])
  names(stk) <- TB_GCMS[ok]
  if (is.null(cell_area_km2)) cell_area_km2 <- prod(terra::res(stk)) / 1e6

  r_avg <- terra::app(stk, fun = mean, na.rm = TRUE); names(r_avg) <- "HS"
  r_sd  <- terra::app(stk, fun = sd,   na.rm = TRUE); names(r_sd)  <- "SD"
  r_cv  <- r_sd / r_avg
  r_cv  <- terra::ifel(r_avg < 1e-6, NA, r_cv)        ## avoid divide-by-near-zero
  names(r_cv) <- "CV"

  ## Save
  fn_avg <- file.path(TB_OUT_HS_AVG, sprintf("%s.tif", scen_avg$scen_key[i]))
  fn_cv  <- file.path(TB_OUT_HS_AVG,
                      sprintf("%s_cv.tif", scen_avg$scen_key[i]))
  terra::writeRaster(r_avg, fn_avg, overwrite = TRUE, datatype = "FLT4S",
                     gdal = c("COMPRESS=DEFLATE","PREDICTOR=2","TILED=YES"))
  terra::writeRaster(r_cv,  fn_cv,  overwrite = TRUE, datatype = "FLT4S",
                     gdal = c("COMPRESS=DEFLATE","PREDICTOR=2","TILED=YES"))

  ## Binary at present threshold
  r_bin <- r_avg >= thr_wmean
  fn_bin <- file.path(TB_OUT_HS_BINARY,
                      sprintf("future_%s.tif", scen_avg$scen_key[i]))
  terra::writeRaster(r_bin, fn_bin, overwrite = TRUE, datatype = "INT1U",
                     gdal = c("COMPRESS=DEFLATE","TILED=YES"))

  avg_rasters[[i]] <- r_avg
  cv_rasters[[i]]  <- r_cv

  vals  <- terra::values(r_avg, mat = FALSE); vals <- vals[!is.na(vals)]
  bvals <- terra::values(r_bin, mat = FALSE); bvals <- bvals[!is.na(bvals)]
  cvals <- terra::values(r_cv,  mat = FALSE); cvals <- cvals[!is.na(cvals)]
  n_suit <- sum(bvals == 1, na.rm = TRUE); n_tot <- length(bvals)

  summary_rows[[i]] <- data.frame(
    scenario     = scen_avg$scen_key[i],
    period       = p, ssp = s,
    n_gcms       = sum(ok),
    HS_mean      = mean(vals),
    HS_median    = median(vals),
    pct_suit     = 100 * n_suit / n_tot,
    area_suit_km2 = n_suit * cell_area_km2,
    CV_mean      = mean(cvals),
    CV_median    = median(cvals),
    CV_p90       = quantile(cvals, 0.90, names = FALSE)
  )
}
tb_toc("6 GCM averages")
summary_df <- do.call(rbind, summary_rows)
tb_save_table(summary_df, "10_gcm_avg_summary")

## Disagreement detail: per-scenario fraction of cells where GCMs disagree about
## suitability (|≥thr| varies across GCMs)
tb_log_section("GCM disagreement table")

disagree_rows <- vector("list", nrow(scen_avg))
for (i in seq_len(nrow(scen_avg))) {
  p <- scen_avg$period[i]; s <- scen_avg$ssp[i]
  paths <- file.path(TB_OUT_HS_FUTURE, sprintf("%s_%s_%s.tif", p, TB_GCMS, s))
  ok <- file.exists(paths)
  if (sum(ok) < 2) next
  stk <- terra::rast(paths[ok])
  bin <- stk >= thr_wmean
  n_agree_all <- terra::app(bin, fun = function(v) all(v == v[1]))
  vals <- terra::values(n_agree_all, mat = FALSE); vals <- vals[!is.na(vals)]
  disagree_rows[[i]] <- data.frame(
    scenario        = scen_avg$scen_key[i],
    period          = p, ssp = s,
    pct_cells_agree = 100 * mean(vals == 1),
    pct_cells_disagree = 100 * mean(vals == 0)
  )
}
disagree_df <- do.call(rbind, disagree_rows)
tb_save_table(disagree_df, "10_gcm_disagreement")

## ----------------------------------------------------------------------------
## FIGURES
## ----------------------------------------------------------------------------

## Common map config
world_sf <- tryCatch(
  rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") |>
    sf::st_transform(TB_CRS_PROJ),
  error = function(e) NULL)
tr_mask_shp <- file.path(TB_DATA_ROOT, "TR_mask", "TR_mask.shp")
tr_mask_sf <- if (file.exists(tr_mask_shp))
  sf::st_read(tr_mask_shp, quiet = TRUE) |> sf::st_transform(TB_CRS_PROJ) else NULL

template <- avg_rasters[[which(!sapply(avg_rasters, is.null))[1]]]
e <- terra::ext(template); ext_pad <- 30000
xlim_p <- c(e$xmin - ext_pad, e$xmax + ext_pad)
ylim_p <- c(e$ymin - ext_pad, e$ymax + ext_pad)

.tb_panel <- function(r, title, fill_palette = "mako", fill_dir = 1,
                      fill_label, limits, show_legend = FALSE, oob = scales::squish) {
  p <- ggplot()
  if (!is.null(world_sf)) p <- p + geom_sf(data = world_sf,
                                            fill = TB_FILL_LAND,
                                            color = TB_COLOR_LAND,
                                            linewidth = 0.25)
  p <- p +
    tidyterra::geom_spatraster(data = r, na.rm = TRUE) +
    scale_fill_viridis_c(option = fill_palette, direction = fill_dir,
                         na.value = "transparent", name = fill_label,
                         limits = limits, oob = oob)
  if (!is.null(tr_mask_sf)) p <- p +
    geom_sf(data = tr_mask_sf, fill = NA, color = TB_COLOR_FRAME,
            linewidth = 0.35)
  p <- p +
    coord_sf(xlim = xlim_p, ylim = ylim_p, datum = sf::st_crs(4326),
             expand = FALSE) +
    labs(title = title) +
    theme_trbear(base_size = 11) +
    theme(plot.title = element_text(size = 12, hjust = 0.5,
                                    margin = ggplot2::margin(b = 3)),
          plot.margin = ggplot2::margin(3, 3, 3, 3))
  if (!show_legend) p <- p + theme(legend.position = "none")
  p
}

## ---- fig10a: 6-panel GCM-averaged HS ----------------------------------------
tb_log_section("fig10a GCM-avg HS 6-panel")

## Order: rows = period (near, far) × cols = SSP (126, 370, 585)
plots_a <- vector("list", nrow(scen_avg))
for (i in seq_len(nrow(scen_avg))) {
  r <- avg_rasters[[i]]
  if (is.null(r)) { plots_a[[i]] <- patchwork::plot_spacer(); next }
  ttl <- sprintf("%s — %s",
                 TB_PERIOD_LABELS[scen_avg$period[i]],
                 TB_SSP_LABELS[scen_avg$ssp[i]])
  plots_a[[i]] <- .tb_panel(r, ttl,
                            fill_label = "HS",
                            limits = c(0, 1),
                            show_legend = (i == nrow(scen_avg)))
}
ord_a <- with(scen_avg, order(period, ssp))
plots_a <- plots_a[ord_a]

p10a <- patchwork::wrap_plots(plots_a, ncol = 3, nrow = 2, guides = "collect") +
  patchwork::plot_annotation(
    title    = "Habitat suitability – GCM-averaged future projections",
    theme = theme(plot.title    = element_text(face = "bold", size = 18,
                                               color = TB_COLOR_FRAME)))
tb_save_fig(p10a, "fig10a_gcm_avg_6panel", w = 22, h = 13, subdir = FIG_SUBDIR)

## ---- fig10b: GCM uncertainty (CV) 6-panel -----------------------------------
tb_log_section("fig10b GCM CV 6-panel")

cv_max <- 0.6  ## cap colour scale for legibility
plots_b <- vector("list", nrow(scen_avg))
for (i in seq_len(nrow(scen_avg))) {
  r <- cv_rasters[[i]]
  if (is.null(r)) { plots_b[[i]] <- patchwork::plot_spacer(); next }
  ttl <- sprintf("%s — %s",
                 TB_PERIOD_LABELS[scen_avg$period[i]],
                 TB_SSP_LABELS[scen_avg$ssp[i]])
  plots_b[[i]] <- .tb_panel(r, ttl,
                            fill_palette = "rocket", fill_dir = -1,
                            fill_label = "CV\n(SD/mean)",
                            limits = c(0, cv_max),
                            show_legend = (i == nrow(scen_avg)))
}
plots_b <- plots_b[ord_a]

p10b <- patchwork::wrap_plots(plots_b, ncol = 3, nrow = 2, guides = "collect") +
  patchwork::plot_annotation(
    title    = "Inter-GCM uncertainty — coefficient of variation",
    theme = theme(plot.title    = element_text(face = "bold", size = 18,
                                               color = TB_COLOR_FRAME)))
tb_save_fig(p10b, "fig10b_gcm_uncertainty_6panel", w = 22, h = 13, subdir = FIG_SUBDIR)

tb_log_session()
tb_log("10_gcm_average DONE")
