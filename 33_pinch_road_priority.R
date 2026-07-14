## ============================================================================
## 33_pinch_road_priority.R
## Project: TR_Bear_Connectivity
## Purpose: Turn the corridor x paved-road intersections into a RANKED, actionable
##   list of mitigation priority sites (candidate wildlife crossings). A pinch
##   point matters most where (i) a lot of modelled movement funnels through it
##   (high corridor KDE intensity), (ii) human-bear conflict risk is high, and
##   (iii) it lies outside protected areas. We combine these into a priority
##   score and rank discrete road-crossing segments.
##
## Method (present scenario):
##   - Top-5% corridor (UNICOR KDE), paved roads (motorway/trunk/primary/
##     secondary), conflict-risk surface (conflict ENMTML W_MEAN ensemble).
##   - Pinch cells = corridor ∩ road. Cluster contiguous pinch cells into
##     segments (terra::patches). Per segment: mean KDE intensity, mean conflict
##     risk, area, centroid (lon/lat), % inside PA, nearest province.
##   - priority_score = z(KDE_intensity) + z(conflict_risk) + (unprotected ? +0.5)
##     ranked descending; top sites reported as crossing-structure candidates.
##
## Outputs (tables/):
##   33_pinch_priority.csv        ranked pinch segments (the deliverable table)
## Figures (figures/33_pinch_priority/):
##   fig33a_priority_map.png      top-N priority pinch segments on corridor+roads
##   fig33b_priority_bars.png     priority score decomposition for top-N
## ============================================================================

suppressPackageStartupMessages({
  library(terra); library(sf); library(dplyr); library(tidyr); library(ggplot2)
  library(patchwork); library(rnaturalearth); library(tidyterra); library(ggrepel)
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
tb_log_init("33_pinch_road_priority")

FIG_SUBDIR    <- "33_pinch_priority"
TOP_PCT       <- 0.05
TOP_N         <- 20
PAVED_CLASSES <- c("motorway", "trunk", "primary", "secondary")

## ----------------------------------------------------------------------------
## 1) Corridor (continuous intensity + top-5% mask)
## ----------------------------------------------------------------------------
tb_log_section("Corridor + roads + conflict")
.read_aaigrid <- function(src) {
  proxy <- tempfile(fileext = ".asc"); file.copy(src, proxy, overwrite = TRUE)
  r <- terra::rast(proxy)
  if (is.na(terra::crs(r)) || terra::crs(r) == "") terra::crs(r) <- TB_CRS_PROJ
  r <- terra::ifel(r < 0, 0, r); names(r) <- "kde"; r
}
kde_src <- list.files(file.path(TB_OUT_UNICOR_DIR, "present", "results"),
                      pattern = "\\.kdepaths$", full.names = TRUE)[1]
kde <- .read_aaigrid(kde_src)
vv  <- terra::values(kde)[, 1]; pos <- which(!is.na(vv) & vv > 0)
thr <- stats::quantile(vv[pos], 1 - TOP_PCT, names = FALSE)
corr_mask <- terra::ifel(kde >= thr, 1, NA); names(corr_mask) <- "corr"

## roads
roads_gpkg <- file.path(TB_OUT_VECTORS, "roads_paved.gpkg")
if (file.exists(roads_gpkg)) {
  roads <- sf::st_read(roads_gpkg, quiet = TRUE) |> sf::st_transform(TB_CRS_PROJ)
} else {
  tb_log("roads_paved.gpkg missing — reading from source shp", "WARN")
  cls <- paste(sprintf("'%s'", PAVED_CLASSES), collapse = ",")
  sql <- sprintf("SELECT fclass, name FROM \"%s\" WHERE fclass IN (%s)",
                 tools::file_path_sans_ext(basename(TB_ROADS_SHP)), cls)
  roads <- tryCatch(sf::st_read(TB_ROADS_SHP, query = sql, quiet = TRUE),
                    error = function(e) {
                      r <- sf::st_read(TB_ROADS_SHP, quiet = TRUE)
                      r[r$fclass %in% PAVED_CLASSES, ]
                    }) |> sf::st_transform(TB_CRS_PROJ)
}
road_mask <- terra::rasterize(terra::vect(roads), kde, field = 1, background = 0,
                              touches = TRUE)
names(road_mask) <- "road"

## conflict-risk surface
conf_tif <- file.path(TB_OUT_ENMTML_CONFLICT, "Ensemble", "W_MEAN",
                      "Ursus_arctos_conflict.tif")
conf <- if (file.exists(conf_tif)) {
  cr <- terra::rast(conf_tif)
  if (!terra::compareGeom(kde, cr, stopOnError = FALSE))
    cr <- terra::resample(cr, kde, method = "bilinear")
  names(cr) <- "conflict"; cr
} else { tb_log("conflict raster missing", "WARN"); NULL }

## PA mask
pa_gpkg <- file.path(TB_OUT_VECTORS, "pa_combined.gpkg")
pa_mask <- if (file.exists(pa_gpkg)) {
  pa <- sf::st_read(pa_gpkg, quiet = TRUE) |> sf::st_transform(TB_CRS_PROJ)
  m <- terra::rasterize(terra::vect(pa), kde, field = 1, background = 0); names(m) <- "pa"; m
} else { tb_log("pa_combined.gpkg missing — PA flag NA", "WARN"); NULL }

## ----------------------------------------------------------------------------
## 2) Pinch cells -> segments
## ----------------------------------------------------------------------------
tb_log_section("Pinch segments")
pinch <- terra::ifel((corr_mask == 1) & (road_mask == 1), 1, NA); names(pinch) <- "pinch"
seg <- terra::patches(pinch, directions = 8, zeroAsNA = TRUE); names(seg) <- "seg"
nseg <- length(stats::na.omit(unique(terra::values(seg)[, 1])))
tb_log(sprintf("pinch cells=%d | segments=%d",
               sum(!is.na(terra::values(pinch)[, 1])), nseg))

cell_km2 <- prod(terra::res(kde)) / 1e6
zonal_mean <- function(x) {
  z <- terra::zonal(x, seg, fun = "mean", na.rm = TRUE); names(z) <- c("seg", "v"); z
}
kz <- zonal_mean(kde); names(kz)[2] <- "kde_mean"
fz <- if (!is.null(conf)) { z <- zonal_mean(conf); names(z)[2] <- "conflict_mean"; z } else NULL
pz <- if (!is.null(pa_mask)) { z <- zonal_mean(pa_mask); names(z)[2] <- "pa_frac"; z } else NULL
nz <- terra::zonal(terra::cellSize(seg, unit = "km"), seg, fun = "sum", na.rm = TRUE)
names(nz) <- c("seg", "area_km2")
## centroids
seg_pts <- terra::as.polygons(seg) |> terra::centroids()
seg_xy  <- terra::crds(seg_pts)
seg_ids <- terra::values(seg_pts)[, 1]
cent <- data.frame(seg = seg_ids, x = seg_xy[, 1], y = seg_xy[, 2])

dat <- cent |>
  dplyr::left_join(kz, by = "seg") |>
  dplyr::left_join(nz, by = "seg")
if (!is.null(fz)) dat <- dplyr::left_join(dat, fz, by = "seg") else dat$conflict_mean <- NA
if (!is.null(pz)) dat <- dplyr::left_join(dat, pz, by = "seg") else dat$pa_frac <- NA

## lon/lat + nearest province (best effort)
ll <- sf::st_as_sf(dat, coords = c("x","y"), crs = TB_CRS_PROJ) |> sf::st_transform(TB_CRS_WGS)
cc <- sf::st_coordinates(ll); dat$lon <- cc[,1]; dat$lat <- cc[,2]
prov <- tryCatch({
  st <- rnaturalearth::ne_states(country = "Turkey", returnclass = "sf") |>
    sf::st_transform(TB_CRS_WGS)
  idx <- sf::st_nearest_feature(ll, st)
  nm <- if ("name" %in% names(st)) st$name else st[[1]]
  nm[idx]
}, error = function(e) rep(NA_character_, nrow(dat)))
dat$province <- prov

## ----------------------------------------------------------------------------
## 3) Priority score
## ----------------------------------------------------------------------------
tb_log_section("Priority scoring")
z <- function(v) { v <- as.numeric(v); s <- stats::sd(v, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(v))); (v - mean(v, na.rm = TRUE)) / s }
dat$unprotected <- ifelse(is.na(dat$pa_frac), TRUE, dat$pa_frac < 0.5)
dat$priority_score <- z(dat$kde_mean) +
  (if (!all(is.na(dat$conflict_mean))) z(dat$conflict_mean) else 0) +
  ifelse(dat$unprotected, 0.5, 0)
dat <- dat |> dplyr::arrange(dplyr::desc(priority_score)) |>
  dplyr::mutate(rank = dplyr::row_number())
tb_save_table(dat, "33_pinch_priority")

top <- dat |> dplyr::slice_head(n = min(TOP_N, nrow(dat)))

## ----------------------------------------------------------------------------
## 4) FIGURES
## ----------------------------------------------------------------------------
tb_log_section("Figures")
tr_shp <- file.path(TB_DATA_ROOT, "TR_mask", "TR_mask.shp")
tr_sf  <- sf::st_read(tr_shp, quiet = TRUE) |> sf::st_transform(TB_CRS_PROJ)
world_sf <- tryCatch(rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") |>
                       sf::st_transform(TB_CRS_PROJ), error = function(e) NULL)
corr_show <- terra::mask(corr_mask, terra::vect(tr_sf))
top_sf <- sf::st_as_sf(top, coords = c("x","y"), crs = TB_CRS_PROJ)
e <- terra::ext(kde); pad <- 30000
xl <- c(e$xmin - pad, e$xmax + pad); yl <- c(e$ymin - pad, e$ymax + pad)

p33a <- ggplot()
if (!is.null(world_sf)) p33a <- p33a +
  geom_sf(data = world_sf, fill = "#F2F2F2", color = "#CFCFCF", linewidth = 0.3)
p33a <- p33a +
  tidyterra::geom_spatraster(data = corr_show, na.rm = TRUE) +
  scale_fill_gradient(low = "#BBD4E6", high = "#BBD4E6", na.value = "transparent",
                      guide = "none") +
  geom_sf(data = roads, color = "#9E9E9E", linewidth = 0.2, alpha = 0.6) +
  geom_sf(data = tr_sf, fill = NA, color = TB_COLOR_FRAME, linewidth = 0.5) +
  geom_sf(data = top_sf, aes(size = priority_score, fill = priority_score),
          shape = 21, color = "black", stroke = 0.5) +
  scale_fill_viridis_c(option = "inferno", direction = -1, name = "Priority") +
  scale_size_continuous(range = c(3, 9), guide = "none") +
  ggrepel::geom_text_repel(data = top,
                           aes(x = x, y = y, label = rank),
                           size = 3, color = "black", max.overlaps = 30,
                           box.padding = 0.3) +
  coord_sf(xlim = xl, ylim = yl, datum = sf::st_crs(4326), expand = FALSE) +
  tb_map_decorations() +
  labs(title = "Priority road-crossing mitigation sites",
       subtitle = sprintf("Top-%d pinch segments where the top-5%% corridor meets a paved road, ranked by flow intensity × conflict risk × exposure.", nrow(top))) +
  theme_trbear()
tb_save_fig(p33a, "fig33a_priority_map", w = 15, h = 9, subdir = FIG_SUBDIR)

## ---- fig33b: priority decomposition bars ----------------------------------
top_b <- top |>
  dplyr::mutate(z_kde = z(dat$kde_mean)[rank],
                z_conf = if (!all(is.na(dat$conflict_mean))) z(dat$conflict_mean)[rank] else 0,
                lbl = factor(sprintf("#%d %s", rank, ifelse(is.na(province), "", province)),
                             levels = rev(sprintf("#%d %s", rank, ifelse(is.na(province), "", province))))) |>
  dplyr::select(lbl, `Flow (KDE)` = z_kde, `Conflict risk` = z_conf) |>
  tidyr::pivot_longer(-lbl, names_to = "component", values_to = "z")
p33b <- ggplot(top_b, aes(z, lbl, fill = component)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("Flow (KDE)" = "#0072B2", "Conflict risk" = "#9E2A2B"),
                    name = "Component (z-score)") +
  labs(title = sprintf("Priority-score components — top %d crossing sites", nrow(top)),
       subtitle = "Stacked standardized contributions; an unprotected bonus (+0.5) is added to the total score.",
       x = "Standardized contribution (z)", y = NULL) +
  theme_trbear_bar(base_size = 12)
tb_save_fig(p33b, "fig33b_priority_bars", w = 13, h = 8, subdir = FIG_SUBDIR)

tb_log_session()
tb_log("33_pinch_road_priority DONE")
