## ============================================================================
## 19_roads_overlay.R
## Project: TR_Bear_Connectivity
## Purpose: Road–corridor conflict overlay.
##   Uses OSM gis_osm_roads_free_1.shp, retained classes:
##     motorway, trunk, primary, secondary  (paved high-traffic)
##   Reports per scenario:
##     - Road km inside top-5% corridor cells
##     - % of corridor cells crossed by a paved road
##     - Pinch-point count: corridor cells crossed by ≥ 1 road segment
##
## Outputs (tables/):
##   19_roads_overlay.csv        scenario → corridor cells, intersected cells,
##                                          road km, pinch %
##
## Vectors (vectors/):
##   roads_paved.gpkg            filtered + projected paved roads
##   road_pinchpoints_present.gpkg  cells where road crosses top-5% corridor
##
## Figures (figures/19_roads_overlay/):
##   fig19a_roads_overview.png   paved roads + corridor (present)
##   fig19b_pinch_density.png    road-km density inside top-5% corridor cells
##   fig19c_pinchpoint_map.png   pinch-point cells highlighted, present
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
tb_log_init("19_roads_overlay")

FIG_SUBDIR <- "19_roads_overlay"

PAVED_CLASSES <- c("motorway", "trunk", "primary", "secondary")

scenarios <- c("present",
               sprintf("%s_%s", rep(TB_PERIODS, each = length(TB_SSPS)),
                                rep(TB_SSPS,    times = length(TB_PERIODS))))

## ----------------------------------------------------------------------------
## 1) Read + filter OSM roads
## ----------------------------------------------------------------------------
tb_log_section("Load OSM roads")

if (!file.exists(TB_ROADS_SHP)) {
  tb_log(sprintf("Roads shapefile missing: %s", TB_ROADS_SHP), "ERROR")
  tb_log_session(); quit(status = 1)
}

## Bounding-box filter — read TR mask first, transform to roads CRS, then read
## roads with bbox to keep memory under control on 1.2 GB shapefile.
tr_mask_shp <- file.path(TB_DATA_ROOT, "TR_mask", "TR_mask.shp")
tr_mask_sf  <- if (file.exists(tr_mask_shp))
  sf::st_read(tr_mask_shp, quiet = TRUE) |> sf::st_transform(TB_CRS_PROJ) else NULL

## Detect roads CRS
roads_layer <- sf::st_layers(TB_ROADS_SHP)
roads_crs   <- roads_layer$crs[[1]]
tb_log(sprintf("OSM roads CRS: %s", sf::st_crs(roads_crs)$input))

if (!is.null(tr_mask_sf)) {
  bb <- sf::st_bbox(sf::st_transform(tr_mask_sf, roads_crs))
  wkt_filter <- sf::st_as_text(sf::st_as_sfc(bb))
} else {
  wkt_filter <- NULL
}

## SQL filter to drop everything but the 4 paved classes at read time.
sql_filter <- sprintf(
  "SELECT fclass, name FROM \"%s\" WHERE fclass IN (%s)",
  tools::file_path_sans_ext(basename(TB_ROADS_SHP)),
  paste(sprintf("'%s'", PAVED_CLASSES), collapse = ", "))

roads <- tryCatch(
  sf::st_read(TB_ROADS_SHP, query = sql_filter,
              wkt_filter = wkt_filter, quiet = TRUE),
  error = function(e) {
    tb_log(sprintf("SQL filter failed (%s); falling back to full read.",
                    conditionMessage(e)), "WARN")
    r <- sf::st_read(TB_ROADS_SHP, quiet = TRUE,
                     wkt_filter = wkt_filter)
    r[r$fclass %in% PAVED_CLASSES, ]
  })
tb_log(sprintf("paved roads read: %d features", nrow(roads)))

roads <- sf::st_transform(roads, TB_CRS_PROJ)
if (!is.null(tr_mask_sf)) {
  roads <- sf::st_intersection(roads, sf::st_union(tr_mask_sf))
  tb_log(sprintf("after TR clip: %d features", nrow(roads)))
}

tot_road_km <- as.numeric(sum(sf::st_length(roads))) / 1000
tb_log(sprintf("total paved road length within TR = %.0f km", tot_road_km))

tb_save_vector(roads[, "fclass"], "roads_paved", subdir = NULL)

## ----------------------------------------------------------------------------
## 2) Rasterize: per-cell road length (km) on KDE grid
## ----------------------------------------------------------------------------
tb_log_section("Rasterize road density")

.find_one <- function(dir, pattern) {
  hits <- list.files(dir, pattern = pattern, full.names = TRUE,
                     ignore.case = TRUE)
  if (length(hits)) hits[1] else NA_character_
}
.read_aaigrid <- function(src) {
  if (is.na(src) || !file.exists(src)) return(NULL)
  proxy <- tempfile(fileext = ".asc"); file.copy(src, proxy, overwrite = TRUE)
  r <- terra::rast(proxy)
  if (is.na(terra::crs(r)) || terra::crs(r) == "") terra::crs(r) <- TB_CRS_PROJ
  r <- terra::ifel(r < 0, 0, r); names(r) <- "value"; r
}

kde_raw <- lapply(scenarios, function(s) {
  res_dir <- file.path(TB_OUT_UNICOR_DIR, s, "results")
  .read_aaigrid(.find_one(res_dir, "\\.kdepaths$"))
})
names(kde_raw) <- scenarios
ref_kde <- kde_raw[["present"]]
for (s in scenarios) {
  if (s == "present" || is.null(kde_raw[[s]])) next
  if (!terra::compareGeom(ref_kde, kde_raw[[s]], stopOnError = FALSE)) {
    kde_raw[[s]] <- terra::resample(kde_raw[[s]], ref_kde, method = "bilinear")
  }
}

## Global top-5% threshold
v_all <- unlist(lapply(kde_raw, function(r) {
  if (is.null(r)) NULL else terra::values(r, mat = FALSE)
}))
v_all <- v_all[!is.na(v_all) & v_all > 0]
top5_thr <- quantile(v_all, 0.95, names = FALSE)
tb_log(sprintf("top-5%% global threshold = %.4g", top5_thr))

cell_km2 <- prod(terra::res(ref_kde)) / 1e6

## Rasterise road length per cell (km). terra::rasterize for SpatVector lines
## with fun="length" was added in 1.7-3 — but the safe approach is to use
## terra::extract by cell. We use rasterize + length="line" if available; else
## a chunked sf-based approach.
roads_v <- terra::vect(roads)

road_len_km_r <- tryCatch({
  ## Newer terra: fun = "length"
  rr <- terra::rasterize(roads_v, ref_kde, fun = "length",
                         background = 0, touches = FALSE)
  rr / 1000        # metres → km per cell
}, error = function(e) {
  tb_log(sprintf("fun='length' rasterize failed (%s); using extract-by-cell fallback.",
                  conditionMessage(e)), "WARN")
  cells_with_road <- terra::cells(ref_kde, roads_v)[, "cell"]
  ## Approximate by counting touched cells × cell-side
  rr <- ref_kde * 0; names(rr) <- "road_km"
  tab <- table(cells_with_road)
  rr[as.integer(names(tab))] <- as.numeric(tab) * (terra::res(ref_kde)[1] / 1000)
  rr
})
names(road_len_km_r) <- "road_km"

## Binary "road present" mask
road_mask <- road_len_km_r > 0

## ----------------------------------------------------------------------------
## 3) Overlay stats per scenario
## ----------------------------------------------------------------------------
tb_log_section("Overlay stats")

overlay_rows <- list()
for (s in scenarios) {
  k <- kde_raw[[s]]
  if (is.null(k)) next
  corr <- k >= top5_thr
  corr_total   <- sum(terra::values(corr), na.rm = TRUE)
  corr_w_road  <- sum(terra::values(corr & road_mask), na.rm = TRUE)
  road_km_in_corr <- sum(terra::values(road_len_km_r * corr), na.rm = TRUE)

  overlay_rows[[s]] <- data.frame(
    scenario             = s,
    corridor_total_km2   = corr_total * cell_km2,
    corridor_with_road   = corr_w_road,
    pct_corr_with_road   = 100 * corr_w_road / max(corr_total, 1),
    road_km_in_corridor  = road_km_in_corr,
    road_km_total        = tot_road_km,
    pct_roads_in_corr    = 100 * road_km_in_corr / max(tot_road_km, 1)
  )
  tb_log(sprintf("[%s] Corr=%.0f km² | road-cells=%.0f%% | road-km in corr=%.0f",
                 s, overlay_rows[[s]]$corridor_total_km2,
                 overlay_rows[[s]]$pct_corr_with_road,
                 overlay_rows[[s]]$road_km_in_corridor))
}
overlay_df <- do.call(rbind, overlay_rows)
tb_save_table(overlay_df, "19_roads_overlay")

## ----------------------------------------------------------------------------
## 4) Pinch-point cells (present) — save as polygon vector
## ----------------------------------------------------------------------------
tb_log_section("Pinch points")

corr_pres <- (kde_raw[["present"]] >= top5_thr) & road_mask
corr_pres <- terra::mask(corr_pres, kde_raw[["present"]])
pinch_v <- terra::as.polygons(corr_pres, dissolve = FALSE, na.rm = TRUE)
pinch_v <- pinch_v[terra::values(pinch_v)[, 1] == 1, ]
if (nrow(pinch_v)) {
  pinch_sf <- sf::st_as_sf(pinch_v)
  tb_save_vector(pinch_sf, "road_pinchpoints_present", subdir = NULL)
  tb_log(sprintf("pinch-point cells: %d", nrow(pinch_sf)))
}

## ----------------------------------------------------------------------------
## 5) FIGURES
## ----------------------------------------------------------------------------
tb_log_section("Figures")

world_sf <- tryCatch(
  rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") |>
    sf::st_transform(TB_CRS_PROJ),
  error = function(e) NULL)

e <- terra::ext(ref_kde); pad <- 30000
xl <- c(e$xmin - pad, e$xmax + pad)
yl <- c(e$ymin - pad, e$ymax + pad)

.clip <- function(r) if (is.null(tr_mask_sf)) r else terra::mask(r, terra::vect(tr_mask_sf))

## ---- fig19a: roads + corridor (present) -------------------------------------
corr_pres_fac <- (kde_raw[["present"]] >= top5_thr) |> terra::as.int() |> terra::as.factor()
levels(corr_pres_fac) <- data.frame(id = c(0, 1),
                                     class = c("Matrix", "Top-5% corridor"))
names(corr_pres_fac) <- "class"
corr_pres_fac <- .clip(corr_pres_fac)

p19a <- ggplot()
if (!is.null(world_sf)) p19a <- p19a +
  geom_sf(data = world_sf, fill = "#E8E8E8",
          color = "#7C8A93", linewidth = 0.4)
p19a <- p19a +
  tidyterra::geom_spatraster(data = corr_pres_fac, na.rm = TRUE) +
  scale_fill_manual(values = c("Matrix"          = "#F2F2F2",
                               "Top-5% corridor" = "#9E2A2B"),
                    na.translate = FALSE, name = "Corridor")
if (!is.null(tr_mask_sf)) p19a <- p19a +
  geom_sf(data = tr_mask_sf, fill = NA,
          color = TB_COLOR_FRAME, linewidth = 0.5)
p19a <- p19a +
  geom_sf(data = roads, color = "#1F3A93", linewidth = 0.25, alpha = 0.85) +
  coord_sf(xlim = xl, ylim = yl, datum = sf::st_crs(4326), expand = FALSE) +
  tb_map_decorations() +
  labs(title    = "Paved roads over present top-5% bear corridor",
       subtitle = sprintf(
         "OSM motorway/trunk/primary/secondary (%.0f km total). Red = corridor backbone.",
         tot_road_km)) +
  theme_trbear()
tb_save_fig(p19a, "fig19a_roads_overview", w = 14, h = 9, subdir = FIG_SUBDIR)

## ---- fig19b: per-cell road density INSIDE present corridor ------------------
corr_pres_bin <- kde_raw[["present"]] >= top5_thr
road_in_corr  <- road_len_km_r * corr_pres_bin
road_in_corr  <- .clip(road_in_corr)
road_in_corr  <- terra::ifel(road_in_corr <= 0, NA, road_in_corr)

p19b <- ggplot()
if (!is.null(world_sf)) p19b <- p19b +
  geom_sf(data = world_sf, fill = "#E8E8E8",
          color = "#7C8A93", linewidth = 0.4)
p19b <- p19b +
  tidyterra::geom_spatraster(data = road_in_corr, na.rm = TRUE) +
  scale_fill_viridis_c(option = "rocket", direction = -1,
                       trans = "sqrt",
                       na.value = "transparent",
                       name = "Road km\nper cell")
if (!is.null(tr_mask_sf)) p19b <- p19b +
  geom_sf(data = tr_mask_sf, fill = NA,
          color = TB_COLOR_FRAME, linewidth = 0.5)
p19b <- p19b +
  coord_sf(xlim = xl, ylim = yl, datum = sf::st_crs(4326), expand = FALSE) +
  tb_map_decorations() +
  labs(title    = "Paved-road density inside present top-5% corridor",
       subtitle = sprintf(
         "Cells coloured by road km within the cell (1 km² grid). Total: %.0f road km in corridor (%.1f%% of TR paved network).",
         overlay_df$road_km_in_corridor[overlay_df$scenario == "present"],
         overlay_df$pct_roads_in_corr[overlay_df$scenario == "present"])) +
  theme_trbear()
tb_save_fig(p19b, "fig19b_pinch_density", w = 14, h = 9, subdir = FIG_SUBDIR)

## ---- fig19c: pinch-point cells (binary highlight) ---------------------------
pinch_r <- (kde_raw[["present"]] >= top5_thr) & road_mask
pinch_r <- terra::as.factor(terra::as.int(pinch_r))
levels(pinch_r) <- data.frame(id = c(0, 1),
                               class = c("Other", "Road × corridor"))
names(pinch_r) <- "class"
pinch_r <- .clip(pinch_r)

p19c <- ggplot()
if (!is.null(world_sf)) p19c <- p19c +
  geom_sf(data = world_sf, fill = "#E8E8E8",
          color = "#7C8A93", linewidth = 0.4)
p19c <- p19c +
  tidyterra::geom_spatraster(data = pinch_r, na.rm = TRUE) +
  scale_fill_manual(values = c("Other"           = "#F2F2F2",
                               "Road × corridor" = "#9E2A2B"),
                    na.translate = FALSE, name = "Pinch")
if (!is.null(tr_mask_sf)) p19c <- p19c +
  geom_sf(data = tr_mask_sf, fill = NA,
          color = TB_COLOR_FRAME, linewidth = 0.5)
p19c <- p19c +
  coord_sf(xlim = xl, ylim = yl, datum = sf::st_crs(4326), expand = FALSE) +
  tb_map_decorations() +
  labs(title    = "Road-corridor pinch-point cells (present)",
       subtitle = sprintf(
         "Red = 1 km² cells where a paved road crosses the top-5%% corridor (%.0f%% of corridor cells affected).",
         overlay_df$pct_corr_with_road[overlay_df$scenario == "present"])) +
  theme_trbear()
tb_save_fig(p19c, "fig19c_pinchpoint_map", w = 14, h = 9, subdir = FIG_SUBDIR)

tb_log_session()
tb_log("19_roads_overlay DONE")
