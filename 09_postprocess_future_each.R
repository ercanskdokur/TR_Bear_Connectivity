## ============================================================================
## 09_postprocess_future_each.R
## Project: TR_Bear_Connectivity
## Purpose: All 18 future HS rasters (2 periods × 3 SSPs × 3 GCMs) — W_MEAN.
##   - Copy to derived/hs_future/
##   - Per-scenario summary table
##   - Two figure panels: raw HS (18 maps) and delta vs present (18 maps)
## ============================================================================

suppressPackageStartupMessages({
  library(terra); library(sf); library(ggplot2); library(dplyr); library(tidyr)
  library(tidyterra); library(patchwork); library(rnaturalearth)
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
tb_log_init("09_postprocess_future_each")

SP <- "Ursus_arctos"
FIG_SUBDIR <- "09_future"

## ----------------------------------------------------------------------------
## 1) Inputs
## ----------------------------------------------------------------------------
tb_log_section("Load present + 18 future")

## Prefer derived copy if 08 has run, else fall back directly to ENMTML output
.pres_path1 <- file.path(TB_OUT_HS_PRESENT, "wmean.tif")
.pres_path2 <- file.path(TB_OUT_ENMTML, "Ensemble", "W_MEAN", paste0(SP, ".tif"))
present_wmean <- terra::rast(if (file.exists(.pres_path1)) .pres_path1 else .pres_path2)
if (terra::nlyr(present_wmean) > 1) present_wmean <- present_wmean[[1]]
names(present_wmean) <- "HS"
tb_log(sprintf("present source: %s", terra::sources(present_wmean)[1]))

scen_grid <- expand.grid(
  period = TB_PERIODS, gcm = TB_GCMS, ssp = TB_SSPS,
  stringsAsFactors = FALSE
) |>
  mutate(scen = sprintf("%s_%s_%s", period, gcm, ssp))

scen_grid <- scen_grid[order(scen_grid$period, scen_grid$ssp, scen_grid$gcm), ]
rownames(scen_grid) <- NULL

scen_grid$src_path <- file.path(
  TB_OUT_ENMTML, "Projection", scen_grid$scen,
  "Ensemble", "W_MEAN", paste0(SP, ".tif")
)
missing <- !file.exists(scen_grid$src_path)
if (any(missing)) {
  tb_log(sprintf("MISSING %d scenarios: %s", sum(missing),
                 paste(scen_grid$scen[missing], collapse = ", ")), "WARN")
}

tb_log(sprintf("scenarios: %d total, %d present, %d missing",
               nrow(scen_grid), sum(!missing), sum(missing)))

## TR mask + world for basemap
world_sf <- tryCatch(
  rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") |>
    sf::st_transform(TB_CRS_PROJ),
  error = function(e) NULL)
tr_mask_shp <- file.path(TB_DATA_ROOT, "TR_mask", "TR_mask.shp")
tr_mask_sf <- if (file.exists(tr_mask_shp))
  sf::st_read(tr_mask_shp, quiet = TRUE) |> sf::st_transform(TB_CRS_PROJ) else NULL

e <- terra::ext(present_wmean); ext_pad <- 30000
xlim_p <- c(e$xmin - ext_pad, e$xmax + ext_pad)
ylim_p <- c(e$ymin - ext_pad, e$ymax + ext_pad)

## Threshold for "suitable" from present W_MEAN
thr_tbl <- read.table(file.path(TB_OUT_ENMTML, "Thresholds_Ensemble.txt"),
                      header = TRUE, sep = "\t", stringsAsFactors = FALSE)
thr_wmean <- thr_tbl$THR_VALUE[thr_tbl$Ensemble == "W_MEAN"]

cell_area_km2 <- prod(terra::res(present_wmean)) / 1e6

## ----------------------------------------------------------------------------
## 2) Copy / load all 18, summarise
## ----------------------------------------------------------------------------
tb_log_section("Process 18 futures")

summary_rows <- vector("list", nrow(scen_grid))
fut_stack    <- vector("list", nrow(scen_grid))
fut_delta    <- vector("list", nrow(scen_grid))

tb_tic()
for (i in seq_len(nrow(scen_grid))) {
  s <- scen_grid$scen[i]
  src <- scen_grid$src_path[i]
  if (!file.exists(src)) { tb_log(sprintf("skip %s", s), "WARN"); next }
  r <- terra::rast(src); names(r) <- "HS"

  ## Save copy
  dst <- file.path(TB_OUT_HS_FUTURE, paste0(s, ".tif"))
  terra::writeRaster(r, dst, overwrite = TRUE, datatype = "FLT4S",
                     gdal = c("COMPRESS=DEFLATE","PREDICTOR=2","TILED=YES"))

  ## Delta
  d <- r - present_wmean
  names(d) <- "delta"

  fut_stack[[i]] <- r
  fut_delta[[i]] <- d

  vals  <- terra::values(r, mat = FALSE); vals <- vals[!is.na(vals)]
  n_suit <- sum(vals >= thr_wmean, na.rm = TRUE)
  n_tot  <- length(vals)

  vals_d <- terra::values(d, mat = FALSE); vals_d <- vals_d[!is.na(vals_d)]

  summary_rows[[i]] <- data.frame(
    scenario    = s,
    period      = scen_grid$period[i],
    gcm         = scen_grid$gcm[i],
    ssp         = scen_grid$ssp[i],
    HS_mean     = mean(vals),
    HS_median   = median(vals),
    pct_suit    = 100 * n_suit / n_tot,
    area_suit_km2 = n_suit * cell_area_km2,
    delta_mean  = mean(vals_d),
    delta_median = median(vals_d),
    pct_gain_cells = 100 * mean(vals_d >  0.05, na.rm = TRUE),
    pct_loss_cells = 100 * mean(vals_d < -0.05, na.rm = TRUE)
  )
}
tb_toc("18 futures processed")

summary_df <- do.call(rbind, summary_rows)
tb_save_table(summary_df, "09_future_each_summary")

## ----------------------------------------------------------------------------
## FIGURES — 18-panel grids
## ----------------------------------------------------------------------------

## A consistent panel order: row = period × SSP (6 rows), col = GCM (3 cols)
make_label_df <- function(grid) {
  grid |>
    mutate(
      period_lbl = TB_PERIOD_LABELS[period],
      ssp_lbl    = TB_SSP_LABELS[ssp],
      gcm_lbl    = TB_GCM_LABELS[gcm],
      row_lbl    = sprintf("%s — %s", period_lbl, ssp_lbl)
    )
}
scen_grid <- make_label_df(scen_grid)

row_levels <- unique(scen_grid[order(scen_grid$period, scen_grid$ssp), "row_lbl"])
col_levels <- TB_GCMS

## ---- fig09a: 18 raw HS maps -------------------------------------------------
tb_log_section("fig09a 18 raw HS panel")

.tb_panel_hs <- function(r, title, fill_palette = "mako", fill_dir = 1,
                         limits = c(0, 1), fill_label = "HS",
                         show_legend = FALSE) {
  p <- ggplot()
  if (!is.null(world_sf)) p <- p + geom_sf(data = world_sf,
                                            fill = TB_FILL_LAND,
                                            color = TB_COLOR_LAND,
                                            linewidth = 0.2)
  p <- p +
    tidyterra::geom_spatraster(data = r, na.rm = TRUE) +
    scale_fill_viridis_c(option = fill_palette, direction = fill_dir,
                         na.value = "transparent", name = fill_label,
                         limits = limits)
  if (!is.null(tr_mask_sf)) p <- p +
    geom_sf(data = tr_mask_sf, fill = NA, color = TB_COLOR_FRAME,
            linewidth = 0.25)
  p <- p +
    coord_sf(xlim = xlim_p, ylim = ylim_p, datum = sf::st_crs(4326),
             expand = FALSE) +
    labs(title = title) +
    theme_trbear(base_size = 9) +
    theme(plot.title = element_text(size = 10, hjust = 0.5,
                                    margin = ggplot2::margin(b = 2)),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.margin = ggplot2::margin(1, 1, 1, 1))
  if (!show_legend) p <- p + theme(legend.position = "none")
  p
}

plots_a <- vector("list", nrow(scen_grid))
for (i in seq_len(nrow(scen_grid))) {
  r <- fut_stack[[i]]
  if (is.null(r)) {
    plots_a[[i]] <- patchwork::plot_spacer(); next
  }
  ttl <- sprintf("%s — %s | %s",
                 TB_PERIOD_LABELS[scen_grid$period[i]],
                 TB_SSP_LABELS[scen_grid$ssp[i]],
                 TB_GCM_LABELS[scen_grid$gcm[i]])
  plots_a[[i]] <- .tb_panel_hs(r, ttl, show_legend = (i == nrow(scen_grid)))
}

## Order plots into 6 rows × 3 cols by (row_lbl, gcm)
ord_idx <- with(scen_grid, order(match(row_lbl, row_levels),
                                  match(gcm, col_levels)))
plots_a <- plots_a[ord_idx]

p09a <- patchwork::wrap_plots(plots_a, ncol = 3, nrow = 6, guides = "collect") +
  patchwork::plot_annotation(
    title    = "Future habitat suitability",
    theme = theme(plot.title    = element_text(face = "bold", size = 18,
                                               color = TB_COLOR_FRAME)))
tb_save_fig(p09a, "fig09a_future_18panel", w = 14, h = 24, subdir = FIG_SUBDIR)

## ---- fig09b: 18 delta maps (future − present) -------------------------------
tb_log_section("fig09b 18 delta panel")

delta_limit <- 0.5  ## cap for visualization

.tb_panel_delta <- function(r, title, show_legend = FALSE) {
  p <- ggplot()
  if (!is.null(world_sf)) p <- p + geom_sf(data = world_sf,
                                            fill = TB_FILL_LAND,
                                            color = TB_COLOR_LAND,
                                            linewidth = 0.2)
  p <- p +
    tidyterra::geom_spatraster(data = r, na.rm = TRUE) +
    scale_fill_gradient2(low = "#D55E00", mid = "#F7F7F7", high = "#009E73",
                         midpoint = 0,
                         limits = c(-delta_limit, delta_limit),
                         oob = scales::squish,
                         na.value = "transparent",
                         name = "Δ HS\n(future − present)")
  if (!is.null(tr_mask_sf)) p <- p +
    geom_sf(data = tr_mask_sf, fill = NA, color = TB_COLOR_FRAME,
            linewidth = 0.25)
  p <- p +
    coord_sf(xlim = xlim_p, ylim = ylim_p, datum = sf::st_crs(4326),
             expand = FALSE) +
    labs(title = title) +
    theme_trbear(base_size = 9) +
    theme(plot.title = element_text(size = 10, hjust = 0.5,
                                    margin = ggplot2::margin(b = 2)),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.margin = ggplot2::margin(1, 1, 1, 1))
  if (!show_legend) p <- p + theme(legend.position = "none")
  p
}

plots_b <- vector("list", nrow(scen_grid))
for (i in seq_len(nrow(scen_grid))) {
  d <- fut_delta[[i]]
  if (is.null(d)) { plots_b[[i]] <- patchwork::plot_spacer(); next }
  ttl <- sprintf("%s — %s | %s",
                 TB_PERIOD_LABELS[scen_grid$period[i]],
                 TB_SSP_LABELS[scen_grid$ssp[i]],
                 TB_GCM_LABELS[scen_grid$gcm[i]])
  plots_b[[i]] <- .tb_panel_delta(d, ttl, show_legend = (i == nrow(scen_grid)))
}
plots_b <- plots_b[ord_idx]

p09b <- patchwork::wrap_plots(plots_b, ncol = 3, nrow = 6, guides = "collect") +
  patchwork::plot_annotation(
    title    = "Change in habitat suitability",
    theme = theme(plot.title    = element_text(face = "bold", size = 18,
                                               color = TB_COLOR_FRAME)))
tb_save_fig(p09b, "fig09b_future_delta_18panel", w = 14, h = 24, subdir = FIG_SUBDIR)

tb_log_session()
tb_log("09_postprocess_future_each DONE")
