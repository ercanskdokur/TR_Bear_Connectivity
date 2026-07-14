## ============================================================================
## 08_postprocess_present.R
## Project: TR_Bear_Connectivity
## Purpose: Present habitat suitability (HS) — final maps + summary tables.
##   - Copy/symlink ensemble rasters into derived/hs_present/
##   - Compute % suitable area, mean HS, median HS, IQR
##   - Build publication-grade maps: main W_MEAN, 8-algo panel, binary, MEAN/W_MEAN
## ============================================================================

suppressPackageStartupMessages({
  library(terra); library(sf); library(ggplot2); library(dplyr)
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
tb_log_init("08_postprocess_present")

SP    <- "Ursus_arctos"
ALGOS <- TB_ENM_ALGORITHMS
FIG_SUBDIR <- "08_present"

## ----------------------------------------------------------------------------
## 1) Inputs: ensemble rasters + per-algo + thresholds + presences
## ----------------------------------------------------------------------------
tb_log_section("Load inputs")

ens_wmean <- terra::rast(file.path(TB_OUT_ENMTML, "Ensemble", "W_MEAN", paste0(SP, ".tif")))
ens_mean  <- terra::rast(file.path(TB_OUT_ENMTML, "Ensemble", "MEAN",   paste0(SP, ".tif")))
ens_wmean_bin <- terra::rast(file.path(TB_OUT_ENMTML, "Ensemble", "W_MEAN", "MAX_TSS", paste0(SP, ".tif")))
ens_mean_bin  <- terra::rast(file.path(TB_OUT_ENMTML, "Ensemble", "MEAN",   "MAX_TSS", paste0(SP, ".tif")))

names(ens_wmean)    <- "HS"
names(ens_mean)     <- "HS"
names(ens_wmean_bin) <- "binary"
names(ens_mean_bin)  <- "binary"

thr_tbl <- read.table(file.path(TB_OUT_ENMTML, "Thresholds_Ensemble.txt"),
                      header = TRUE, sep = "\t", stringsAsFactors = FALSE)
thr_wmean <- thr_tbl$THR_VALUE[thr_tbl$Ensemble == "W_MEAN"]
thr_mean  <- thr_tbl$THR_VALUE[thr_tbl$Ensemble == "MEAN"]
tb_log(sprintf("thresholds: W_MEAN=%.4f, MEAN=%.4f", thr_wmean, thr_mean))

algo_stack <- terra::rast(lapply(ALGOS, function(a) {
  r <- terra::rast(file.path(TB_OUT_ENMTML, "Algorithm", a, paste0(SP, ".tif")))
  names(r) <- a; r
}))

occ <- read.table(file.path(TB_OUT_ENMTML, "Occurrences_Cleaned.txt"),
                  header = TRUE, sep = "\t", stringsAsFactors = FALSE)
occ_sf <- sf::st_as_sf(occ, coords = c("x","y"), crs = TB_CRS_PROJ)
tb_log(sprintf("presences cleaned: %d", nrow(occ_sf)))

## Country basemap (TR + neighbors) for map context
world_sf <- tryCatch(
  rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") |>
    sf::st_transform(TB_CRS_PROJ),
  error = function(e) { tb_log("world basemap failed; skipping", "WARN"); NULL }
)

tr_mask_shp <- file.path(TB_DATA_ROOT, "TR_mask", "TR_mask.shp")
tr_mask_sf <- if (file.exists(tr_mask_shp)) sf::st_read(tr_mask_shp, quiet = TRUE) |>
                                              sf::st_transform(TB_CRS_PROJ) else NULL

## Map extent — slight pad around bbox of ensemble
ext_pad <- 30000  # 30 km
e <- terra::ext(ens_wmean)
xlim_p <- c(e$xmin - ext_pad, e$xmax + ext_pad)
ylim_p <- c(e$ymin - ext_pad, e$ymax + ext_pad)

## ----------------------------------------------------------------------------
## 2) Copy ensembles to derived/hs_present/ for downstream
## ----------------------------------------------------------------------------
tb_log_section("Copy to derived/hs_present")

terra::writeRaster(ens_wmean,     file.path(TB_OUT_HS_PRESENT, "wmean.tif"),
                   overwrite = TRUE, datatype = "FLT4S",
                   gdal = c("COMPRESS=DEFLATE","PREDICTOR=2","TILED=YES"))
terra::writeRaster(ens_mean,      file.path(TB_OUT_HS_PRESENT, "mean.tif"),
                   overwrite = TRUE, datatype = "FLT4S",
                   gdal = c("COMPRESS=DEFLATE","PREDICTOR=2","TILED=YES"))
terra::writeRaster(ens_wmean_bin, file.path(TB_OUT_HS_BINARY,  "present_wmean.tif"),
                   overwrite = TRUE, datatype = "INT1U",
                   gdal = c("COMPRESS=DEFLATE","TILED=YES"))
terra::writeRaster(ens_mean_bin,  file.path(TB_OUT_HS_BINARY,  "present_mean.tif"),
                   overwrite = TRUE, datatype = "INT1U",
                   gdal = c("COMPRESS=DEFLATE","TILED=YES"))

## ----------------------------------------------------------------------------
## 3) Summary table
## ----------------------------------------------------------------------------
tb_log_section("Summary stats")

cell_area_km2 <- prod(terra::res(ens_wmean)) / 1e6  ## 1 km × 1 km = 1.0 km²

summary_row <- function(r_cont, r_bin, label) {
  vals  <- terra::values(r_cont, mat = FALSE)
  vals  <- vals[!is.na(vals)]
  bvals <- terra::values(r_bin, mat = FALSE)
  bvals <- bvals[!is.na(bvals)]
  n_suit <- sum(bvals == 1, na.rm = TRUE)
  n_tot  <- length(bvals)
  data.frame(
    ensemble       = label,
    area_total_km2 = n_tot * cell_area_km2,
    area_suit_km2  = n_suit * cell_area_km2,
    pct_suit       = 100 * n_suit / n_tot,
    HS_mean        = mean(vals),
    HS_median      = median(vals),
    HS_q25         = quantile(vals, 0.25, names = FALSE),
    HS_q75         = quantile(vals, 0.75, names = FALSE),
    HS_max         = max(vals)
  )
}
sum_df <- rbind(
  summary_row(ens_wmean, ens_wmean_bin, "W_MEAN"),
  summary_row(ens_mean,  ens_mean_bin,  "MEAN")
)
tb_save_table(sum_df, "08_present_summary")

## ----------------------------------------------------------------------------
## FIGURES
## ----------------------------------------------------------------------------

## Map helper
## Clip raster to TR mask before plotting so neighbouring countries (basemap)
## are visible — ENMTML output sometimes contains stray low-HS pixels outside
## the mask that hide the basemap fill.
.tb_clip_to_tr <- function(r) {
  if (is.null(tr_mask_sf)) return(r)
  terra::mask(r, terra::vect(tr_mask_sf))
}

.tb_map_hs <- function(r, title, subtitle,
                       fill_palette = "mako", fill_dir = 1,
                       fill_label   = "Habitat\nsuitability",
                       limits = c(0, 1),
                       overlay_pres = TRUE, show_legend = TRUE) {
  r_c <- .tb_clip_to_tr(r)
  p <- ggplot()
  if (!is.null(world_sf)) {
    p <- p + geom_sf(data = world_sf, fill = "#E8E8E8",
                     color = "#7C8A93", linewidth = 0.4)
  }
  p <- p +
    tidyterra::geom_spatraster(data = r_c, na.rm = TRUE) +
    scale_fill_viridis_c(option = fill_palette, direction = fill_dir,
                         na.value = "transparent", name = fill_label,
                         limits = limits)
  if (!is.null(tr_mask_sf)) {
    p <- p + geom_sf(data = tr_mask_sf, fill = NA,
                     color = TB_COLOR_FRAME, linewidth = 0.5)
  }
  if (overlay_pres) {
    p <- p + geom_sf(data = occ_sf, color = "#D55E00", fill = "white",
                     shape = 21, size = 0.7, stroke = 0.25, alpha = 0.8)
  }
  p <- p +
    coord_sf(xlim = xlim_p, ylim = ylim_p, datum = sf::st_crs(4326),
             expand = FALSE) +
    tb_map_decorations() +
    labs(title = title, subtitle = subtitle) +
    theme_trbear()
  if (!show_legend) p <- p + theme(legend.position = "none")
  p
}

## ---- fig08a: main W_MEAN with presences -------------------------------------
tb_log_section("fig08a present main")

p08a <- .tb_map_hs(
  ens_wmean,
  title    = "Brown bear habitat suitability across Türkiye for present-day",
  subtitle = NULL
)
tb_save_fig(p08a, "fig08a_present_hs_main", w = 13, h = 8.5, subdir = FIG_SUBDIR)

## ---- fig08b: 8-algorithm panel ----------------------------------------------
tb_log_section("fig08b 8-algo panel")

algo_plots <- lapply(ALGOS, function(a) {
  .tb_map_hs(algo_stack[[a]], title = a, subtitle = NULL,
             overlay_pres = FALSE, show_legend = (a == ALGOS[length(ALGOS)])) +
    theme(plot.title = element_text(size = 15, hjust = 0.5, face = "bold",
                                    margin = ggplot2::margin(b = 4)),
          axis.text = element_text(size = 8),
          plot.margin = ggplot2::margin(6, 6, 6, 6))
})
p08b <- patchwork::wrap_plots(algo_plots, ncol = 4, nrow = 2,
                              guides = "collect") +
  patchwork::plot_annotation(
    title = "Per-algorithm habitat suitability for present-day",
    theme = theme(plot.title = element_text(face = "bold", size = 22,
                                            color = TB_COLOR_FRAME)))
tb_save_fig(p08b, "fig08b_present_8algos_panel", w = 26, h = 13, subdir = FIG_SUBDIR)

## ---- fig08c: binary suitable/unsuitable -------------------------------------
tb_log_section("fig08c binary")

bin_factor <- .tb_clip_to_tr(ens_wmean_bin)
levels(bin_factor) <- data.frame(value = c(0, 1),
                                 class = c("Unsuitable", "Suitable"))

p08c <- ggplot()
if (!is.null(world_sf)) p08c <- p08c +
  geom_sf(data = world_sf, fill = "#E8E8E8",
          color = "#7C8A93", linewidth = 0.4)
p08c <- p08c +
  tidyterra::geom_spatraster(data = bin_factor, na.rm = TRUE) +
  scale_fill_manual(values = TB_PAL_BINARY, na.translate = FALSE,
                    name = "Class")
if (!is.null(tr_mask_sf)) p08c <- p08c +
  geom_sf(data = tr_mask_sf, fill = NA,
          color = TB_COLOR_FRAME, linewidth = 0.5)
p08c <- p08c +
  geom_sf(data = occ_sf, color = "#000000", fill = "white",
          shape = 21, size = 0.7, stroke = 0.25, alpha = 0.8) +
  coord_sf(xlim = xlim_p, ylim = ylim_p, datum = sf::st_crs(4326),
           expand = FALSE) +
  tb_map_decorations() +
  labs(title    = "Brown bear suitable habitats (binary) across Türkiye",
       subtitle = sprintf(
         "W_MEAN ≥ %.3f (MAX_TSS). Suitable area: %.0f km² (%.1f%%).",
         thr_wmean, sum_df$area_suit_km2[sum_df$ensemble == "W_MEAN"],
         sum_df$pct_suit[sum_df$ensemble == "W_MEAN"])) +
  theme_trbear()
tb_save_fig(p08c, "fig08c_present_binary", w = 13, h = 8.5, subdir = FIG_SUBDIR)

## ---- fig08d: MEAN vs W_MEAN side-by-side ------------------------------------
tb_log_section("fig08d mean vs wmean")

p_mean   <- .tb_map_hs(ens_mean,  title = "MEAN ensemble",
                       subtitle = sprintf("Threshold = %.3f", thr_mean),
                       overlay_pres = FALSE, show_legend = FALSE)
p_wmean  <- .tb_map_hs(ens_wmean, title = "W_MEAN (TSS-weighted)",
                       subtitle = sprintf("Threshold = %.3f", thr_wmean),
                       overlay_pres = FALSE, show_legend = TRUE)
p08d <- (p_mean | p_wmean) +
  patchwork::plot_annotation(
    title = "Ensemble comparison: simple MEAN vs TSS-weighted W_MEAN",
    subtitle = "W_MEAN gives more weight to algorithms with higher TSS (RDF, BRT, MXD top-ranked).",
    theme = theme(plot.title    = element_text(face = "bold", size = 16,
                                               color = TB_COLOR_FRAME),
                  plot.subtitle = element_text(size = 12, color = TB_COLOR_AXIS)))
tb_save_fig(p08d, "fig08d_mean_vs_wmean", w = 18, h = 8.5, subdir = FIG_SUBDIR)

tb_log_session()
tb_log("08_postprocess_present DONE")
