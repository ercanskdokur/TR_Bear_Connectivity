## ============================================================================
## 31_corridor_centroid_shift.R
## Project: TR_Bear_Connectivity
## Purpose: Does the connectivity BACKBONE shift northward and upslope under
##   climate change, mirroring the habitat range shift reported by
##   Sıkdokur et al. (2025)? We quantify the displacement of the corridor
##   network's centre of mass and the change in its mean elevation/latitude.
##
## Method:
##   - Load each scenario's UNICOR KDE corridor surface (.kdepaths), align to
##     the present grid (bilinear resample), define each scenario's own
##     top-5% corridor backbone.
##   - Connectivity-weighted centroid = KDE-weighted mean of backbone cell
##     coordinates. Shift vector (present -> future): distance (km) + bearing.
##   - Mean corridor elevation (Elevation held constant across scenarios, so a
##     rise means corridors migrate to higher ground) and mean / distribution
##     of corridor latitude per scenario.
##
## Outputs (tables/):
##   31_corridor_centroid_shift.csv  scenario -> centroid lon/lat, shift_km,
##                                    bearing, mean_elev_m, mean_lat
## Figures (figures/31_centroid/):
##   fig31a_centroid_shift_map.png   centroids + present->future arrows on TR
##   fig31b_corridor_elev_lat.png    mean corridor elevation & latitude per scen
## ============================================================================

suppressPackageStartupMessages({
  library(terra); library(sf); library(dplyr); library(tidyr)
  library(ggplot2); library(patchwork); library(rnaturalearth); library(tidyterra)
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
tb_log_init("31_corridor_centroid_shift")

FIG_SUBDIR <- "31_centroid"
TOP_PCT    <- 0.05

scenarios <- c("present",
               sprintf("%s_%s", rep(TB_PERIODS, each = length(TB_SSPS)),
                                rep(TB_SSPS,    times = length(TB_PERIODS))))
scen_label <- c(
  "present"          = "Present",
  "2041_2070_ssp126" = "2070s SSP126",  "2041_2070_ssp370" = "2070s SSP370",
  "2041_2070_ssp585" = "2070s SSP585",  "2071_2100_ssp126" = "2100s SSP126",
  "2071_2100_ssp370" = "2100s SSP370",  "2071_2100_ssp585" = "2100s SSP585")

## ----------------------------------------------------------------------------
## 1) UNICOR KDE loader (Arc/Info ASCII grid with .kdepaths suffix)
## ----------------------------------------------------------------------------
tb_log_section("Load UNICOR KDE rasters")
.read_aaigrid <- function(src) {
  if (is.na(src) || !file.exists(src)) return(NULL)
  proxy <- tempfile(fileext = ".asc"); file.copy(src, proxy, overwrite = TRUE)
  r <- terra::rast(proxy)
  if (is.na(terra::crs(r)) || terra::crs(r) == "") terra::crs(r) <- TB_CRS_PROJ
  r <- terra::ifel(r < 0, 0, r); names(r) <- "value"; r
}
.find_kde <- function(s) {
  rd <- file.path(TB_OUT_UNICOR_DIR, s, "results")
  if (!dir.exists(rd)) return(NA_character_)
  h <- list.files(rd, pattern = "\\.kdepaths$", full.names = TRUE)
  if (length(h)) h[1] else NA_character_
}
kde <- lapply(scenarios, function(s) .read_aaigrid(.find_kde(s)))
names(kde) <- scenarios
ref <- kde[["present"]]
if (is.null(ref)) { tb_log("present KDE missing", "ERROR"); tb_log_session(); quit(status = 1) }
for (s in scenarios) {
  if (s == "present" || is.null(kde[[s]])) next
  if (!terra::compareGeom(ref, kde[[s]], stopOnError = FALSE))
    kde[[s]] <- terra::resample(kde[[s]], ref, method = "bilinear")
}

## ----------------------------------------------------------------------------
## 2) Elevation on the corridor grid
## ----------------------------------------------------------------------------
tb_log_section("Elevation raster")
get_elev <- function() {
  ps_stack <- file.path(TB_OUT_RASTERS, "present_stack.tif")
  if (file.exists(ps_stack)) {
    st <- terra::rast(ps_stack)
    if ("Elevation" %in% names(st)) {
      e <- st[["Elevation"]]
      if (terra::compareGeom(ref, e, stopOnError = FALSE)) return(e)
      return(terra::resample(e, ref, method = "bilinear"))
    }
  }
  ef <- file.path(TB_PREDICTORS_TIF_PRESENT, "Elevation.tif")
  if (file.exists(ef)) {
    e <- terra::rast(ef)
    e <- terra::project(e, ref, method = "bilinear")
    return(e)
  }
  tb_log("Elevation raster not found — elevation columns will be NA", "WARN")
  NULL
}
elev <- get_elev()
if (!is.null(elev)) names(elev) <- "elev"

## ----------------------------------------------------------------------------
## 3) Per-scenario centroid, shift, elevation, latitude
## ----------------------------------------------------------------------------
tb_log_section("Centroid + shift metrics")

backbone_stats <- function(r) {
  v <- terra::values(r)[, 1]
  pos <- which(!is.na(v) & v > 0)
  thr <- stats::quantile(v[pos], 1 - TOP_PCT, names = FALSE)
  cells <- which(!is.na(v) & v >= thr)
  xy <- terra::xyFromCell(r, cells)
  w  <- v[cells]
  cx <- sum(xy[,1] * w) / sum(w); cy <- sum(xy[,2] * w) / sum(w)
  ev <- if (!is.null(elev)) terra::extract(elev, xy)[, 1] else rep(NA_real_, length(cells))
  list(cells = cells, xy = xy, w = w, cx = cx, cy = cy,
       mean_elev = stats::weighted.mean(ev, w, na.rm = TRUE),
       elev = ev)
}

to_lonlat <- function(x, y) {
  p <- sf::st_as_sf(data.frame(x = x, y = y), coords = c("x","y"), crs = TB_CRS_PROJ) |>
    sf::st_transform(TB_CRS_WGS)
  sf::st_coordinates(p)
}

bb <- lapply(scenarios, function(s) if (is.null(kde[[s]])) NULL else backbone_stats(kde[[s]]))
names(bb) <- scenarios

pres_c <- c(bb[["present"]]$cx, bb[["present"]]$cy)
rows <- list(); lat_box <- list()
for (s in scenarios) {
  if (is.null(bb[[s]])) next
  cll <- to_lonlat(bb[[s]]$cx, bb[[s]]$cy)
  dx <- bb[[s]]$cx - pres_c[1]; dy <- bb[[s]]$cy - pres_c[2]
  ## weighted-mean latitude of backbone cells
  ll_all <- to_lonlat(bb[[s]]$xy[,1], bb[[s]]$xy[,2])
  rows[[s]] <- data.frame(
    scenario = s, label = scen_label[s],
    centroid_lon = cll[1], centroid_lat = cll[2],
    centroid_x = bb[[s]]$cx, centroid_y = bb[[s]]$cy,
    shift_km = sqrt(dx^2 + dy^2) / 1000,
    bearing_deg = (atan2(dx, dy) * 180 / pi) %% 360,
    northward_km = dy / 1000,
    mean_elev_m = bb[[s]]$mean_elev,
    mean_lat = stats::weighted.mean(ll_all[,2], bb[[s]]$w))
  ## sample for latitude boxplot (weighted by taking cells; subsample for size)
  idx <- if (length(bb[[s]]$w) > 4000) sample(seq_along(bb[[s]]$w), 4000) else seq_along(bb[[s]]$w)
  lat_box[[s]] <- data.frame(scenario = s, label = scen_label[s],
                             lat = ll_all[idx, 2], elev = bb[[s]]$elev[idx])
}
shift_df <- dplyr::bind_rows(rows)
shift_df$label <- factor(shift_df$label, levels = scen_label[scenarios])
tb_save_table(shift_df, "31_corridor_centroid_shift")
box_df <- dplyr::bind_rows(lat_box)
box_df$label <- factor(box_df$label, levels = scen_label[scenarios])

PAL_SCEN <- c(
  "Present" = "#000000", "2070s SSP126" = "#56B4E9", "2070s SSP370" = "#E69F00",
  "2070s SSP585" = "#D55E00", "2100s SSP126" = "#0072B2",
  "2100s SSP370" = "#CC79A7", "2100s SSP585" = "#9E2A2B")

## ----------------------------------------------------------------------------
## 4) FIGURES
## ----------------------------------------------------------------------------
tb_log_section("Figures")
tr_shp <- file.path(TB_DATA_ROOT, "TR_mask", "TR_mask.shp")
tr_sf  <- sf::st_read(tr_shp, quiet = TRUE) |> sf::st_transform(TB_CRS_PROJ)
world_sf <- tryCatch(rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") |>
                       sf::st_transform(TB_CRS_PROJ), error = function(e) NULL)

## present backbone for context
pres_bb_cells <- bb[["present"]]$cells
pres_back <- terra::rast(ref); terra::values(pres_back) <- NA
pres_back[pres_bb_cells] <- 1
pres_back <- terra::mask(pres_back, terra::vect(tr_sf)); names(pres_back) <- "corr"

cent_sf <- sf::st_as_sf(shift_df, coords = c("centroid_x","centroid_y"), crs = TB_CRS_PROJ)
fut <- shift_df |> dplyr::filter(scenario != "present")
arrows_df <- data.frame(x = pres_c[1], y = pres_c[2],
                        xend = fut$centroid_x, yend = fut$centroid_y,
                        label = fut$label)
e <- terra::ext(ref); pad <- 30000
xl <- c(e$xmin - pad, e$xmax + pad); yl <- c(e$ymin - pad, e$ymax + pad)

p31a <- ggplot()
if (!is.null(world_sf)) p31a <- p31a +
  geom_sf(data = world_sf, fill = "#F2F2F2", color = "#CFCFCF", linewidth = 0.3)
p31a <- p31a +
  tidyterra::geom_spatraster(data = pres_back, na.rm = TRUE) +
  scale_fill_gradient(low = "#BBD4E6", high = "#BBD4E6", na.value = "transparent",
                      guide = "none") +
  geom_sf(data = tr_sf, fill = NA, color = TB_COLOR_FRAME, linewidth = 0.5) +
  geom_segment(data = arrows_df, aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
               linewidth = 0.6, color = "#37474F", alpha = 0.8) +
  geom_sf(data = cent_sf, aes(color = label), size = 4.5) +
  scale_color_manual(values = PAL_SCEN, name = "Scenario centroid") +
  coord_sf(xlim = xl, ylim = yl, datum = sf::st_crs(4326), expand = FALSE) +
  tb_map_decorations() +
  labs(title = "Climate-driven displacement of the connectivity backbone",
       subtitle = sprintf("Light band = present top-%.0f%% corridor. Points = KDE-weighted centroid of each scenario's backbone; arrows = present → future shift.",
                          100 * TOP_PCT)) +
  theme_trbear()
tb_save_fig(p31a, "fig31a_centroid_shift_map", w = 15, h = 9, subdir = FIG_SUBDIR)

## ---- fig31b: corridor elevation + latitude per scenario --------------------
pe <- ggplot(shift_df, aes(label, mean_elev_m, fill = label)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.0f", mean_elev_m)), vjust = -0.5, size = 3) +
  scale_fill_manual(values = PAL_SCEN, guide = "none") +
  labs(title = "A. Mean elevation of corridor backbone", x = NULL, y = "Elevation (m)") +
  theme_trbear_bar(base_size = 11) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
pl <- ggplot(box_df, aes(label, lat, fill = label)) +
  geom_boxplot(outlier.size = 0.3, linewidth = 0.3, alpha = 0.85) +
  scale_fill_manual(values = PAL_SCEN, guide = "none") +
  labs(title = "B. Latitude of corridor cells (northward shift)", x = NULL,
       y = "Latitude (°N)") +
  theme_trbear_bar(base_size = 11) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
p31b <- (pe / pl) +
  patchwork::plot_annotation(
    title = "Are corridors migrating upslope and northward?",
    theme = theme(plot.title = element_text(face = "bold", size = 15, color = TB_COLOR_FRAME)))
tb_save_fig(p31b, "fig31b_corridor_elev_lat", w = 13, h = 10, subdir = FIG_SUBDIR)

tb_log_session()
tb_log("31_corridor_centroid_shift DONE")
