## ============================================================================
## 26_compound_risk.R
## Project: TR_Bear_Connectivity
## Purpose: Identify spatial hotspots where THREE pressures co-occur:
##   (1) Top-5% UNICOR corridor (present, global threshold from 16/15)
##   (2) Paved roads within a 1-km buffer
##   (3) Conflict SDM (W_MEAN, present) above its 75th percentile
##
##   These are priority intervention sites: places where bears must move,
##   where they encounter roads, and where bear-human conflict risk is high.
##
##   Approach:
##     - Build binary layers for each pressure
##     - Multiply / AND them to get a compound-risk binary
##     - Extract connected components (terra::patches)
##     - Save top-N (by area) as gpkg + table; map them on a TR-scale figure
##
## Inputs:
##   present KDE       : <UNICOR_DIR>/present/results/*.kdepaths
##   roads             : <TB_OUT_VECTORS>/roads_paved.gpkg
##   conflict SDM      : <TB_OUT_ENMTML_CONFLICT>/Ensemble/W_MEAN/Ursus_arctos_conflict.tif
##
## Outputs (tables/):
##   26_compound_hotspots.csv          (one row per hotspot polygon: id, area_km2, x, y, ...)
##   26_compound_risk_summary.csv       (total area, % of corridor, etc.)
## Outputs (vectors/):
##   26_compound_risk_hotspots.gpkg
##   26_compound_risk_mask.tif          (binary raster, 1 = compound risk)
## Outputs (figures/26_compound/):
##   fig26a_overview.png                 (TR-wide overview, hotspots labelled)
##   fig26b_top10_inset.png              (zoomed inset of top-10 hotspots)
## ============================================================================

suppressPackageStartupMessages({
  library(terra); library(sf); library(ggplot2); library(dplyr); library(tidyterra)
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
source("00_paths.R"); source("00_helpers.R")
tb_log_init("26_compound_risk")

ROAD_BUFFER_M     <- 1000
CONFLICT_QTILE    <- 0.75
TOP_N_HOTSPOTS    <- 10
MIN_HOTSPOT_KM2   <- 1        # ignore micro-clusters

## ---- Inputs ----------------------------------------------------------------
tb_log_section("Inputs")

## Present top-5% corridor
res_dir <- file.path(TB_OUT_UNICOR_DIR, "present", "results")
kde_src <- list.files(res_dir, pattern = "\\.kdepaths$",
                       full.names = TRUE, ignore.case = TRUE)[1]
stopifnot(file.exists(kde_src))
proxy <- tempfile(fileext = ".asc"); file.copy(kde_src, proxy, overwrite = TRUE)
kde <- terra::rast(proxy)
if (is.na(terra::crs(kde)) || terra::crs(kde) == "") terra::crs(kde) <- TB_CRS_PROJ
kde <- terra::ifel(kde < 0, 0, kde)

## Global threshold (from 16_corridor_thresholds.csv if available, else compute)
thr_csv <- file.path(TB_OUT_TABLES, "16_corridor_thresholds.csv")
if (file.exists(thr_csv)) {
  thr_df <- read.csv(thr_csv)
  top5_thr <- thr_df$threshold[thr_df$top_pct == "top5"]
  if (!length(top5_thr)) top5_thr <- quantile(terra::values(kde),
                                                0.95, na.rm = TRUE)
} else {
  top5_thr <- quantile(terra::values(kde), 0.95, na.rm = TRUE)
}
tb_log(sprintf("top-5%% threshold = %.4g", top5_thr))
corr_bin <- terra::ifel(kde >= top5_thr, 1, 0)

## Conflict SDM
conf_tif <- file.path(TB_OUT_ENMTML_CONFLICT, "Ensemble", "W_MEAN",
                       "Ursus_arctos_conflict.tif")
stopifnot(file.exists(conf_tif))
conf <- terra::rast(conf_tif)
if (is.na(terra::crs(conf)) || terra::crs(conf) == "") terra::crs(conf) <- TB_CRS_PROJ
if (!terra::compareGeom(corr_bin, conf, stopOnError = FALSE))
  conf <- terra::resample(conf, corr_bin, method = "bilinear")
conf_thr <- quantile(terra::values(conf), CONFLICT_QTILE, na.rm = TRUE)
tb_log(sprintf("conflict SDM Q%.0f threshold = %.4g",
                CONFLICT_QTILE * 100, conf_thr))
conf_bin <- terra::ifel(conf >= conf_thr, 1, 0)

## Roads buffer raster
roads_gpkg <- file.path(TB_OUT_VECTORS, "roads_paved.gpkg")
stopifnot(file.exists(roads_gpkg))
roads <- sf::st_read(roads_gpkg, quiet = TRUE) |>
  sf::st_transform(terra::crs(corr_bin))
road_buf <- sf::st_buffer(roads, dist = ROAD_BUFFER_M) |>
  sf::st_union() |> sf::st_make_valid()
road_bin <- terra::rasterize(terra::vect(sf::st_sf(geometry = road_buf)),
                              corr_bin, field = 1, background = 0)

## TR mask
tr_shp <- file.path(TB_DATA_ROOT, "TR_mask", "TR_mask.shp")
tr_sf  <- sf::st_read(tr_shp, quiet = TRUE) |>
  sf::st_transform(terra::crs(corr_bin))
.clip <- function(r) terra::mask(r, terra::vect(tr_sf))

## ---- Compound binary -------------------------------------------------------
tb_log_section("Compound risk")

compound <- corr_bin * road_bin * conf_bin
compound <- .clip(compound)
names(compound) <- "compound"

n_compound_cells <- sum(terra::values(compound, mat = FALSE) == 1, na.rm = TRUE)
n_corridor_cells <- sum(terra::values(.clip(corr_bin), mat = FALSE) == 1, na.rm = TRUE)
cell_km2 <- prod(terra::res(corr_bin)) / 1e6
compound_km2 <- n_compound_cells * cell_km2
corridor_km2 <- n_corridor_cells * cell_km2
pct_corridor <- 100 * compound_km2 / max(corridor_km2, 1)
tb_log(sprintf("compound risk area = %.1f km² (%.2f%% of corridor)",
                compound_km2, pct_corridor))

terra::writeRaster(compound,
                   file.path(TB_OUT_RASTERS, "26_compound_risk_mask.tif"),
                   overwrite = TRUE)

## ---- Connected components → polygons ---------------------------------------
tb_log_section("Hotspot polygons")
clumps <- terra::patches(compound, directions = 8, zeroAsNA = TRUE)
hot_poly <- terra::as.polygons(clumps, dissolve = TRUE, values = TRUE)
hot_sf <- sf::st_as_sf(hot_poly) |>
  dplyr::mutate(area_km2 = as.numeric(sf::st_area(.data$geometry)) / 1e6) |>
  dplyr::filter(.data$area_km2 >= MIN_HOTSPOT_KM2) |>
  dplyr::arrange(dplyr::desc(.data$area_km2)) |>
  dplyr::mutate(rank = dplyr::row_number())

if (nrow(hot_sf) == 0) {
  tb_log("No hotspots ≥ MIN_HOTSPOT_KM2; relaxing minimum to 0.25 km²", "WARN")
  hot_sf <- sf::st_as_sf(hot_poly) |>
    dplyr::mutate(area_km2 = as.numeric(sf::st_area(.data$geometry)) / 1e6) |>
    dplyr::filter(.data$area_km2 >= 0.25) |>
    dplyr::arrange(dplyr::desc(.data$area_km2)) |>
    dplyr::mutate(rank = dplyr::row_number())
}

## Centroids for labeling
centroids <- sf::st_centroid(hot_sf$geometry) |> sf::st_coordinates()
hot_sf$x <- centroids[, 1]; hot_sf$y <- centroids[, 2]

## Reproject to WGS84 for human-readable centroid lon/lat
ll <- hot_sf |> sf::st_transform(TB_CRS_WGS)
ll_xy <- sf::st_coordinates(sf::st_centroid(ll))
hot_sf$lon <- ll_xy[, 1]; hot_sf$lat <- ll_xy[, 2]

sf::st_write(hot_sf, file.path(TB_OUT_VECTORS, "26_compound_risk_hotspots.gpkg"),
             delete_dsn = TRUE, quiet = TRUE)

hot_tbl <- sf::st_drop_geometry(hot_sf)
tb_save_table(hot_tbl, "26_compound_hotspots")

summary_df <- data.frame(
  corridor_km2 = corridor_km2,
  compound_km2 = compound_km2,
  pct_corridor_under_compound_risk = pct_corridor,
  conflict_quantile_threshold = CONFLICT_QTILE,
  road_buffer_m   = ROAD_BUFFER_M,
  n_hotspots = nrow(hot_sf))
tb_save_table(summary_df, "26_compound_risk_summary")

## ---- Fig26a — overview -----------------------------------------------------
tb_log_section("Fig26a overview")
world_sf <- tryCatch(
  rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") |>
    sf::st_transform(TB_CRS_PROJ),
  error = function(e) NULL)

bear_tif <- file.path(TB_OUT_ENMTML, "Ensemble", "W_MEAN", "Ursus_arctos.tif")
bear_r <- if (file.exists(bear_tif)) {
  br <- terra::rast(bear_tif); .clip(br)
} else NULL

p_a <- ggplot()
if (!is.null(world_sf)) p_a <- p_a + geom_sf(data = world_sf, fill = "#E8E8E8",
                                              color = "#7C8A93", linewidth = 0.3)
if (!is.null(bear_r))
  p_a <- p_a + tidyterra::geom_spatraster(data = bear_r, na.rm = TRUE) +
    scale_fill_viridis_c(option = "viridis", limits = c(0, 1),
                          name = "Bear HS", na.value = NA, alpha = 0.7)
tr_bb_a <- sf::st_bbox(tr_sf)
pad_a   <- 30000
p_a <- p_a +
  geom_sf(data = tr_sf, fill = NA, color = TB_COLOR_FRAME, linewidth = 0.4) +
  geom_sf(data = hot_sf, fill = "#9E2A2B", color = "#9E2A2B",
          alpha = 0.85, linewidth = 0.3) +
  ggrepel::geom_label_repel(
    data = hot_sf[hot_sf$rank <= TOP_N_HOTSPOTS, ],
    aes(x = x, y = y, label = sprintf("#%d", rank)),
    size = 3.4, fontface = "bold", label.size = 0.2,
    fill = alpha("white", 0.85),
    box.padding = 0.4, max.overlaps = Inf, segment.size = 0.3) +
  coord_sf(xlim = c(tr_bb_a["xmin"] - pad_a, tr_bb_a["xmax"] + pad_a),
           ylim = c(tr_bb_a["ymin"] - pad_a, tr_bb_a["ymax"] + pad_a),
           datum = sf::st_crs(4326), expand = FALSE) +
  labs(title    = sprintf("Compound-risk hotspots: Corridors for suitability ∩ paved-road buffer (%.0f km) ∩ conflict risky model (Q%.0f)",
                          ROAD_BUFFER_M / 1000, CONFLICT_QTILE * 100),
       subtitle = sprintf("Red = compound risk cells, n = %d hotspots ≥ %.0f km². Top-%d labelled. Compound area = %s km² (%.1f%% of present top-5 corridor).",
                          nrow(hot_sf), MIN_HOTSPOT_KM2, TOP_N_HOTSPOTS,
                          formatC(compound_km2, format = "d", big.mark = ","),
                          pct_corridor)) +
  theme_trbear(base_size = 11)
tb_save_fig(p_a, "fig26a_overview", w = 14, h = 8, subdir = "26_compound")

## ---- Fig26b — top-10 inset table figure ------------------------------------
tb_log_section("Fig26b top-10 bars")
top10 <- hot_sf[hot_sf$rank <= TOP_N_HOTSPOTS, ]
top10$lbl <- factor(sprintf("#%d  (%.2f, %.2f)",
                             top10$rank, top10$lon, top10$lat),
                    levels = sprintf("#%d  (%.2f, %.2f)",
                                      top10$rank, top10$lon, top10$lat))
p_b <- ggplot(sf::st_drop_geometry(top10), aes(area_km2, lbl)) +
  geom_col(fill = "#9E2A2B", color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f km²", area_km2)),
            hjust = -0.1, size = 3.4) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(title = sprintf("Top-%d compound-risk hotspots", TOP_N_HOTSPOTS),
       x = "Hotspot area (km²)", y = "Rank (lon, lat)") +
  theme_trbear_bar(base_size = 11) +
  theme(axis.text.y = element_text(face = "bold"))
tb_save_fig(p_b, "fig26b_top10_bars", w = 10, h = 6, subdir = "26_compound")

tb_log_session()
tb_log("26_compound_risk DONE")
