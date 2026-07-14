## ============================================================================
## 02_predictors_present_to_tif.R
## Project: TR_Bear_Connectivity — ENMTML pipeline
## Purpose: Build the present-day predictor TIFs ENMTML will consume.
##   (1) Build master target grid (TR Albers Equal-Area, 1 km).
##   (2) Reproject 30 source predictors (rTop + rBio + rHum) to target grid.
##   (3) Mask to Turkey landmass.
##   (4) Write one TIF per predictor into TB_PRED_ENMTML_PRESENT (ENMTML input dir).
##   (5) Save target grid + TR mask for downstream scripts.
##   (6) QC figures (per-layer + overview panel).
## Inputs:
##   data/Predictors_TIF/present/<layer>.tif       (pre-converted from RData)
## Outputs:
##   outputs/rasters/target_grid.tif
##   outputs/rasters/tr_landmask.tif
##   outputs/rasters/present_stack.tif
##   data/predictors_enmtml/present/<layer>.tif     ← ENMTML input
##   outputs/rds_files/02_target_grid_info.rds
##   outputs/figures/02_present/02_present_<layer>.png
##   outputs/figures/02_present/02_present_overview.png
## ============================================================================

suppressPackageStartupMessages({
  library(terra); library(sf); library(ggplot2); library(tidyterra)
  library(rnaturalearth); library(patchwork)
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
tb_log_init("02_predictors_present_to_tif")
tb_pkg_versions(c("terra","sf","tidyterra","rnaturalearth"))

OUT_FIG_DIR <- file.path(TB_OUT_FIGURES, "02_present")
dir.create(OUT_FIG_DIR, recursive = TRUE, showWarnings = FALSE)

## ============================================================================
## 1. LOAD PRESENT SOURCE TIFs
## ============================================================================
tb_log_section("1. LOAD PRESENT SOURCE TIFs")
tb_tic()

read_group <- function(nms, label) {
  files <- file.path(TB_PREDICTORS_TIF_PRESENT, sprintf("%s.tif", nms))
  miss  <- nms[!file.exists(files)]
  if (length(miss)) stop(sprintf("[%s] missing TIFs: %s", label, paste(miss, collapse=", ")))
  stk <- terra::rast(files)
  names(stk) <- nms
  tb_log(sprintf("%s: %d layers loaded from %s", label, terra::nlyr(stk), TB_PREDICTORS_TIF_PRESENT))
  stk
}

rTop <- read_group(TB_NAMES_TOP, "rTop")
rBio <- read_group(TB_NAMES_BIO, "rBio")
rHum <- read_group(TB_NAMES_HUM, "rHum")
tb_toc("load")

## ============================================================================
## 2. TURKEY LANDMASS POLYGON
## ============================================================================
tb_log_section("2. TR LANDMASS POLYGON")
tr_wgs <- tryCatch(
  rnaturalearth::ne_countries(country = "Turkey", scale = "large", returnclass = "sf"),
  error = function(e) {
    tb_log("ne_countries large failed; trying medium", "WARN")
    rnaturalearth::ne_countries(country = "Turkey", scale = "medium", returnclass = "sf")
  }
)
tr_proj <- sf::st_transform(tr_wgs, TB_CRS_PROJ)
tb_log(sprintf("TR bbox (m): %s", paste(round(sf::st_bbox(tr_proj)), collapse=",")))

## ============================================================================
## 3. TARGET GRID (1 km, TR Albers)
## ============================================================================
tb_log_section("3. TARGET GRID")

bbox <- sf::st_bbox(tr_proj)
BUFFER_M <- 50000
xmin <- floor((bbox["xmin"] - BUFFER_M) / TB_RES_M) * TB_RES_M
xmax <- ceiling((bbox["xmax"] + BUFFER_M) / TB_RES_M) * TB_RES_M
ymin <- floor((bbox["ymin"] - BUFFER_M) / TB_RES_M) * TB_RES_M
ymax <- ceiling((bbox["ymax"] + BUFFER_M) / TB_RES_M) * TB_RES_M

target <- terra::rast(
  xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
  resolution = TB_RES_M, crs = TB_CRS_PROJ
)
terra::values(target) <- 1L
names(target) <- "target_grid"
tb_log(sprintf("target grid: %d rows × %d cols (%d cells), res=%dm",
               nrow(target), ncol(target), terra::ncell(target), TB_RES_M))

tb_save_raster(target, "target_grid", datatype = "INT1U")

tb_save_rds(list(
  crs = TB_CRS_PROJ, res_m = TB_RES_M,
  extent = c(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
  ncell = terra::ncell(target), buffer_m = BUFFER_M,
  tr_polygon_proj = tr_proj
), "02_target_grid_info")

## TR landmask
tb_log("rasterizing TR landmask")
tr_mask <- terra::rasterize(terra::vect(tr_proj), target, field = 1)
tb_save_raster(tr_mask, "tr_landmask", datatype = "INT1U")

## ============================================================================
## 4. REPROJECT + MASK + WRITE TO ENMTML INPUT DIR
## ============================================================================
tb_log_section("4. REPROJECT + WRITE TIFs FOR ENMTML")
tb_tic()

reproj_one <- function(stk_name, stk, method_default = "bilinear") {
  out_layers <- list()
  for (nm in names(stk)) {
    method <- if (nm == "Aspect") "near" else method_default
    r_src <- stk[[nm]]
    if (terra::crs(r_src) == "") {
      tb_log(sprintf("%s/%s: src CRS empty; assuming WGS84", stk_name, nm), "WARN")
      terra::crs(r_src) <- TB_CRS_WGS
    }
    r_out <- terra::project(r_src, target, method = method, threads = TRUE)
    r_out <- terra::mask(r_out, tr_mask)
    names(r_out) <- nm

    ## ENMTML input directory
    fn_enm <- file.path(TB_PRED_ENMTML_PRESENT, sprintf("%s.tif", nm))
    terra::writeRaster(r_out, fn_enm, overwrite = TRUE, datatype = "FLT4S",
                       gdal = c("COMPRESS=DEFLATE","PREDICTOR=2","TILED=YES"))
    mm <- terra::minmax(r_out)
    tb_log(sprintf("  %s/%s (%s): [%.3f, %.3f] -> %s",
                   stk_name, nm, method, mm[1], mm[2], basename(fn_enm)))
    out_layers[[nm]] <- r_out
  }
  out_layers
}

present_layers <- c(
  reproj_one("rTop", rTop, "bilinear"),
  reproj_one("rBio", rBio, "bilinear"),
  reproj_one("rHum", rHum, "bilinear")
)
tb_toc("reprojection")

## also save stack for QC / quick inspection (outside ENMTML input dir)
present_stack <- terra::rast(present_layers)
terra::writeRaster(present_stack,
                   file.path(TB_OUT_RASTERS, "present_stack.tif"),
                   overwrite = TRUE, datatype = "FLT4S",
                   gdal = c("COMPRESS=DEFLATE","PREDICTOR=2","TILED=YES"))

## ============================================================================
## 5. QC FIGURES
## ============================================================================
tb_log_section("5. QC FIGURES")
tb_tic()

bbox_proj <- sf::st_bbox(tr_proj)
basemap   <- tb_basemap_world(TB_CRS_PROJ, scale_res = "medium")

draw_layer <- function(r, title, sub = NULL) {
  ggplot() +
    geom_sf(data = basemap, fill = TB_FILL_LAND, color = TB_COLOR_LAND, linewidth = 0.3) +
    geom_spatraster(data = r, maxcell = 6e5) +
    scale_fill_viridis_c(option = "viridis", na.value = "transparent", name = names(r)) +
    geom_sf(data = tr_proj, fill = NA, color = TB_COLOR_FRAME, linewidth = 0.4) +
    coord_sf(
      xlim = c(bbox_proj$xmin - 80000, bbox_proj$xmax + 80000),
      ylim = c(bbox_proj$ymin - 80000, bbox_proj$ymax + 80000),
      expand = FALSE, datum = sf::st_crs(4326)
    ) +
    theme_trbear(base_size = 10) + tb_map_decorations() +
    labs(title = title, subtitle = sub %||% "",
         caption = "TR Albers Equal-Area | 1 km")
}

for (nm in names(present_stack)) {
  p <- draw_layer(present_stack[[nm]], title = nm,
                  sub = sprintf("Present-day  |  ENMTML input"))
  tb_save_fig(p, sprintf("02_present_%s", nm), w = 12, h = 7, dpi = 300, subdir = "02_present")
}

## overview small-multiples
mk_thumb <- function(r) {
  ggplot() +
    geom_spatraster(data = r, maxcell = 1.5e5) +
    scale_fill_viridis_c(option = "viridis", na.value = "transparent", guide = "none") +
    geom_sf(data = tr_proj, fill = NA, color = TB_COLOR_FRAME, linewidth = 0.25) +
    coord_sf(xlim = c(bbox_proj$xmin, bbox_proj$xmax),
             ylim = c(bbox_proj$ymin, bbox_proj$ymax),
             expand = FALSE) +
    theme_void() +
    theme(plot.title = element_text(size = 8, face = "bold",
                                    color = TB_COLOR_FRAME, hjust = 0.5)) +
    labs(title = names(r))
}
thumbs <- lapply(names(present_stack), function(nm) mk_thumb(present_stack[[nm]]))
overview <- patchwork::wrap_plots(thumbs, ncol = 6) +
  patchwork::plot_annotation(
    title    = "Present-day predictors (1 km)",
    subtitle = "19 bioclimatic, 4 topographic, 7 human-related predictors",
    theme = theme(plot.title    = element_text(face = "bold", size = 16,
                                                color = TB_COLOR_FRAME),
                  plot.subtitle = element_text(size = 12,
                                                color = TB_COLOR_FRAME))
  )
tb_save_fig(overview, "02_present_overview", w = 18, h = 13, dpi = 220, subdir = "02_present")
tb_toc("figures")

tb_log_session()
tb_log("02_predictors_present_to_tif DONE")
