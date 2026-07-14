## ============================================================================
## 18_pa_overlay.R
## Project: TR_Bear_Connectivity
## Purpose: Protected-area gap analysis for habitat suitability and connectivity
##   under present + 6 future scenarios. Reads 12 layers from PAs.gpkg:
##     hassas_sukutle, millipark, MUHAZAORM, OZELCEVREKORUMA,
##     REKR_KENTORMANI, REKR_MESIREALAN,
##     sulak_MahOnHaSuAl, sulak_Ramsar, sulak_UlnHaSuAl,
##     tabiat_koruma_alani, tabiat_parki, YABANHAYATIGELSAH
##
##   For each scenario it reports:
##     - % of suitable habitat inside PA (any of the 12 layers)
##     - % of top-5% corridor cells inside PA
##     - Gap area (suitable but outside PA, km²)
##     - Per-PA-type breakdown (% of suitable habitat per layer)
##
## Outputs (tables/):
##   18_pa_overlay_summary.csv   scenario × {suitable, corridor} → in/out, km²
##   18_pa_by_layer.csv          PA layer → scenario → habitat km² inside
##
## Figures (figures/18_pa_overlay/):
##   fig18a_pa_overview.png      PA polygons + present binary HS
##   fig18b_coverage_bars.png    Coverage % bars (habitat, corridor) per scenario
##   fig18c_gap_map.png          Present suitable habitat OUTSIDE PA (gap)
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
tb_log_init("18_pa_overlay")

FIG_SUBDIR <- "18_pa_overlay"

PA_LAYERS <- c(
  "hassas_sukutle", "millipark", "MUHAZAORM", "OZELCEVREKORUMA",
  "REKR_KENTORMANI", "REKR_MESIREALAN",
  "sulak_MahOnHaSuAl", "sulak_Ramsar", "sulak_UlnHaSuAl",
  "tabiat_koruma_alani", "tabiat_parki", "YABANHAYATIGELSAH"
)

scenarios <- c("present",
               sprintf("%s_%s", rep(TB_PERIODS, each = length(TB_SSPS)),
                                rep(TB_SSPS,    times = length(TB_PERIODS))))

.bin_path <- function(s) {
  if (s == "present") file.path(TB_OUT_HS_BINARY, "present_wmean.tif")
  else                file.path(TB_OUT_HS_BINARY, sprintf("future_%s.tif", s))
}

## ----------------------------------------------------------------------------
## 1) Read PA layers from File Geodatabase
## ----------------------------------------------------------------------------
tb_log_section("Load PA layers")

avail <- tryCatch(sf::st_layers(TB_PA_GDB)$name, error = function(e) character())
tb_log(sprintf("layers in PAs.gpkg: %d found", length(avail)))

read_pa_layer <- function(layer) {
  if (!(layer %in% avail)) {
    tb_log(sprintf("[%s] not present in PAs.gpkg — skipping", layer), "WARN")
    return(NULL)
  }
  v <- tryCatch(sf::st_read(TB_PA_GDB, layer = layer, quiet = TRUE),
                error = function(e) {
                  tb_log(sprintf("[%s] read failed: %s", layer,
                                  conditionMessage(e)), "WARN")
                  NULL})
  if (is.null(v) || !nrow(v)) return(NULL)
  v <- sf::st_make_valid(v)
  v <- v[sf::st_geometry_type(v) %in%
           c("POLYGON", "MULTIPOLYGON",
             "GEOMETRYCOLLECTION", "CURVEPOLYGON", "MULTISURFACE"), ]
  if (!nrow(v)) return(NULL)
  v <- sf::st_transform(v, TB_CRS_PROJ)
  v$pa_layer <- layer
  ## Keep only essentials for memory
  v[, "pa_layer"]
}

pa_list <- lapply(PA_LAYERS, read_pa_layer)
names(pa_list) <- PA_LAYERS
pa_list <- pa_list[!sapply(pa_list, is.null)]
tb_log(sprintf("loaded %d / %d PA layers", length(pa_list), length(PA_LAYERS)))

pa_all <- do.call(rbind, pa_list)
pa_all <- sf::st_make_valid(pa_all)
tb_log(sprintf("combined PA features: %d (total km² = %.0f)",
               nrow(pa_all),
               as.numeric(sum(sf::st_area(pa_all))) / 1e6))

## Save merged PA polygons for downstream + QC
tb_save_vector(pa_all, "pa_combined", subdir = NULL)

## ----------------------------------------------------------------------------
## 2) Rasterize PA mask once on present grid
## ----------------------------------------------------------------------------
tb_log_section("Rasterize PA mask")

ref_path <- .bin_path("present")
ref <- terra::rast(ref_path)
cell_km2 <- prod(terra::res(ref)) / 1e6

pa_vect    <- terra::vect(pa_all)
pa_mask    <- terra::rasterize(pa_vect, ref, field = 1, background = 0)
names(pa_mask) <- "pa"

## Per-layer rasters (for breakdown)
pa_layer_rasters <- lapply(pa_list, function(v) {
  if (is.null(v)) return(NULL)
  r <- terra::rasterize(terra::vect(v), ref, field = 1, background = 0)
  names(r) <- "pa"; r
})

## ----------------------------------------------------------------------------
## 3) Load aligned KDE rasters (for top-5% corridor mask)
## ----------------------------------------------------------------------------
tb_log_section("Load KDE rasters")

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
  if (s == "present") next
  if (is.null(kde_raw[[s]])) next
  if (!terra::compareGeom(ref_kde, kde_raw[[s]], stopOnError = FALSE)) {
    kde_raw[[s]] <- terra::resample(kde_raw[[s]], ref_kde, method = "bilinear")
  }
}

## Global top-5% threshold (same logic as 15_unicor_post)
v_all <- unlist(lapply(kde_raw, function(r) {
  if (is.null(r)) NULL else terra::values(r, mat = FALSE)
}))
v_all <- v_all[!is.na(v_all) & v_all > 0]
top5_thr <- quantile(v_all, 0.95, names = FALSE)
tb_log(sprintf("top-5%% global threshold = %.4g", top5_thr))

## Resample PA mask onto KDE grid if needed
if (!terra::compareGeom(ref_kde, pa_mask, stopOnError = FALSE)) {
  pa_mask_kde <- terra::resample(pa_mask, ref_kde, method = "near")
} else {
  pa_mask_kde <- pa_mask
}

## ----------------------------------------------------------------------------
## 4) Overlay statistics
## ----------------------------------------------------------------------------
tb_log_section("Overlay stats")

overlay_rows <- list()
for (s in scenarios) {
  bp <- .bin_path(s)
  if (!file.exists(bp)) { tb_log(sprintf("missing %s", bp), "WARN"); next }
  hs <- terra::rast(bp)
  hs <- terra::resample(hs, ref, method = "near")

  hs_suit <- hs == 1
  pa_in_hs <- pa_mask == 1
  inside_hab <- hs_suit & pa_in_hs

  hs_total <- sum(terra::values(hs_suit), na.rm = TRUE)
  hs_in_pa <- sum(terra::values(inside_hab), na.rm = TRUE)

  ## Corridor stats
  k <- kde_raw[[s]]
  if (!is.null(k)) {
    corr <- k >= top5_thr
    corr_total <- sum(terra::values(corr), na.rm = TRUE)
    corr_in_pa <- sum(terra::values(corr & (pa_mask_kde == 1)), na.rm = TRUE)
  } else {
    corr_total <- NA_integer_; corr_in_pa <- NA_integer_
  }

  overlay_rows[[s]] <- data.frame(
    scenario          = s,
    hs_total_km2      = hs_total * cell_km2,
    hs_in_pa_km2      = hs_in_pa * cell_km2,
    hs_pct_in_pa      = 100 * hs_in_pa / max(hs_total, 1),
    hs_gap_km2        = (hs_total - hs_in_pa) * cell_km2,
    corr_total_km2    = corr_total * cell_km2,
    corr_in_pa_km2    = corr_in_pa * cell_km2,
    corr_pct_in_pa    = 100 * corr_in_pa / max(corr_total, 1)
  )
  tb_log(sprintf("[%s] HS %.0f km² (%.1f%% in PA) | Corr %.0f km² (%.1f%% in PA)",
                 s,
                 overlay_rows[[s]]$hs_total_km2,
                 overlay_rows[[s]]$hs_pct_in_pa,
                 overlay_rows[[s]]$corr_total_km2 %||% NA,
                 overlay_rows[[s]]$corr_pct_in_pa %||% NA))
}
overlay_df <- do.call(rbind, overlay_rows)
tb_save_table(overlay_df, "18_pa_overlay_summary")

## ---- Per-layer breakdown (km² of suitable habitat inside each PA layer) -----
tb_log_section("Per-layer breakdown")

layer_rows <- list()
for (s in scenarios) {
  bp <- .bin_path(s)
  if (!file.exists(bp)) next
  hs <- terra::rast(bp); hs <- terra::resample(hs, ref, method = "near")
  for (lay in names(pa_layer_rasters)) {
    pl <- pa_layer_rasters[[lay]]
    if (is.null(pl)) next
    n_in <- sum(terra::values((hs == 1) & (pl == 1)), na.rm = TRUE)
    layer_rows[[paste(s, lay, sep = "_")]] <- data.frame(
      scenario = s, pa_layer = lay,
      hs_in_layer_km2 = n_in * cell_km2)
  }
}
layer_df <- do.call(rbind, layer_rows)
tb_save_table(layer_df, "18_pa_by_layer")

## ----------------------------------------------------------------------------
## 5) FIGURES
## ----------------------------------------------------------------------------
tb_log_section("Figures")

world_sf <- tryCatch(
  rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") |>
    sf::st_transform(TB_CRS_PROJ),
  error = function(e) NULL)
tr_mask_shp <- file.path(TB_DATA_ROOT, "TR_mask", "TR_mask.shp")
tr_mask_sf  <- if (file.exists(tr_mask_shp))
  sf::st_read(tr_mask_shp, quiet = TRUE) |> sf::st_transform(TB_CRS_PROJ) else NULL

e <- terra::ext(ref); pad <- 30000
xl <- c(e$xmin - pad, e$xmax + pad)
yl <- c(e$ymin - pad, e$ymax + pad)

.clip <- function(r) if (is.null(tr_mask_sf)) r else terra::mask(r, terra::vect(tr_mask_sf))

## ---- fig18a: PA polygons over present binary HS -----------------------------
pres_bin <- terra::rast(.bin_path("present")) |> .clip()
pres_bin_fac <- terra::as.factor(terra::as.int(pres_bin))
levels(pres_bin_fac) <- data.frame(id = c(0, 1),
                                    class = c("Unsuitable", "Suitable"))
names(pres_bin_fac) <- "class"

p18a <- ggplot()
if (!is.null(world_sf)) p18a <- p18a +
  geom_sf(data = world_sf, fill = "#E8E8E8",
          color = "#7C8A93", linewidth = 0.4)
p18a <- p18a +
  tidyterra::geom_spatraster(data = pres_bin_fac, na.rm = TRUE) +
  scale_fill_manual(values = TB_PAL_BINARY, na.translate = FALSE,
                    name = "Habitat") +
  geom_sf(data = pa_all, fill = "#1B5E20", color = "#003300",
          alpha = 0.45, linewidth = 0.15)
if (!is.null(tr_mask_sf)) p18a <- p18a +
  geom_sf(data = tr_mask_sf, fill = NA,
          color = TB_COLOR_FRAME, linewidth = 0.5)
p18a <- p18a +
  coord_sf(xlim = xl, ylim = yl, datum = sf::st_crs(4326), expand = FALSE) +
  tb_map_decorations() +
  labs(title    = "Protected areas over present-day bear habitat suitability",
       subtitle = sprintf("Combined area = %s km².",
                          formatC(as.numeric(sum(sf::st_area(pa_all))) / 1e6,
                                   format = "d", big.mark = ","))) +
  theme_trbear()
tb_save_fig(p18a, "fig18a_pa_overview", w = 14, h = 9, subdir = FIG_SUBDIR)

## ---- fig18b: coverage bars --------------------------------------------------
scen_label <- c(
  "present"          = "Present",
  "2041_2070_ssp126" = "2070s SSP126",
  "2041_2070_ssp370" = "2070s SSP370",
  "2041_2070_ssp585" = "2070s SSP585",
  "2071_2100_ssp126" = "2100s SSP126",
  "2071_2100_ssp370" = "2100s SSP370",
  "2071_2100_ssp585" = "2100s SSP585")

bar_df <- overlay_df |>
  mutate(scen_lbl = factor(scen_label[scenario], levels = scen_label[scenarios])) |>
  select(scen_lbl, hs_pct_in_pa, corr_pct_in_pa) |>
  tidyr::pivot_longer(-scen_lbl,
                      names_to = "metric", values_to = "pct") |>
  mutate(metric = factor(metric,
                         levels = c("hs_pct_in_pa", "corr_pct_in_pa"),
                         labels = c("Suitable habitat",
                                    "Top-5% corridor")))

p18b <- ggplot(bar_df, aes(scen_lbl, pct, fill = metric)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7,
           color = "white", linewidth = 0.2) +
  geom_text(aes(label = sprintf("%.1f", pct)),
            position = position_dodge(width = 0.8),
            vjust = -0.3, size = 3.2) +
  scale_fill_manual(values = c("Suitable habitat" = "#009E73",
                               "Top-5% corridor"  = "#9E2A2B"),
                    name = "Layer") +
  labs(title    = "Protected-area coverage of habitat and corridors",
       subtitle = "Per-scenario % of suitable habitat / top-5% corridor cells inside the PA layers.",
       x = NULL, y = "% inside PA") +
  theme_trbear_bar(base_size = 12)
tb_save_fig(p18b, "fig18b_coverage_bars", w = 14, h = 7, subdir = FIG_SUBDIR)

## ---- fig18c: present gap map (suitable & outside PA) ------------------------
gap_r <- (pres_bin == 1) & (pa_mask == 0)
gap_r <- terra::as.factor(terra::as.int(gap_r))
levels(gap_r) <- data.frame(id = c(0, 1),
                             class = c("In-PA / unsuitable",
                                       "Gap: suitable & unprotected"))
names(gap_r) <- "class"
gap_r <- .clip(gap_r)

p18c <- ggplot()
if (!is.null(world_sf)) p18c <- p18c +
  geom_sf(data = world_sf, fill = "#E8E8E8",
          color = "#7C8A93", linewidth = 0.4)
p18c <- p18c +
  tidyterra::geom_spatraster(data = gap_r, na.rm = TRUE) +
  scale_fill_manual(values = c("In-PA / unsuitable"          = "#F2F2F2",
                               "Gap: suitable & unprotected" = "#9E2A2B"),
                    na.translate = FALSE, name = NULL)
if (!is.null(tr_mask_sf)) p18c <- p18c +
  geom_sf(data = tr_mask_sf, fill = NA,
          color = TB_COLOR_FRAME, linewidth = 0.5)
p18c <- p18c +
  coord_sf(xlim = xl, ylim = yl, datum = sf::st_crs(4326), expand = FALSE) +
  tb_map_decorations() +
  labs(title    = "Conservation gap — suitable bear habitat outside protected areas (present)",
       subtitle = sprintf(
         "Red = suitable & unprotected (%.0f km², %.1f%% of total suitable habitat).",
         overlay_df$hs_gap_km2[overlay_df$scenario == "present"],
         100 - overlay_df$hs_pct_in_pa[overlay_df$scenario == "present"])) +
  theme_trbear()
tb_save_fig(p18c, "fig18c_gap_map", w = 14, h = 9, subdir = FIG_SUBDIR)

tb_log_session()
tb_log("18_pa_overlay DONE")
