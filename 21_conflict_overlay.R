## ============================================================================
## 21_conflict_overlay.R
## Project: TR_Bear_Connectivity
## Purpose: Descriptive cross-tabulation of 478 bear–human conflict points with
##   three spatial covariates (NO SDM; complement to the conflict ENMTML run):
##     (a) Biogeographic region (Euro-Siberian / Irano-Turanian / Mediterranean)
##     (b) Protected-area membership (any of the 12 PA layers in PAs.gpkg)
##     (c) Distance to nearest paved road (motorway/trunk/primary/secondary)
##   ... each crossed with Activity type (7 levels: Beekeeping, Livestock, ...).
##
## Inputs:
##   data/points/conflict_table_df.xlsx (478 × 17 with Biogeographic_Region, Activity)
##   data/points/ConflictPoints.txt     (fallback for coordinates if needed)
##   outputs/vectors/pa_combined.gpkg   (from 18_pa_overlay)
##   outputs/vectors/roads_paved.gpkg   (from 19_roads_overlay)
##
## Outputs (tables/):
##   21_conflict_clean.csv           cleaned points with overlay columns
##   21_xtab_bioregion_activity.csv  bioregion × activity counts + col %
##   21_xtab_pa_activity.csv         PA in/out × activity counts + col %
##   21_road_dist_by_activity.csv    distance summary stats per activity
##
## Figures (figures/21_conflict/):
##   fig21a_bioregion_activity_heatmap.png
##   fig21b_pa_activity_bars.png
##   fig21c_road_dist_box.png
## ============================================================================

suppressPackageStartupMessages({
  for (.p in c("readxl")) {
    if (!requireNamespace(.p, quietly = TRUE))
      install.packages(.p, repos = "https://cloud.r-project.org")
  }
  library(terra); library(sf); library(ggplot2); library(dplyr); library(tidyr)
  library(readxl); library(patchwork)
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
tb_log_init("21_conflict_overlay")

FIG_SUBDIR <- "21_conflict"

## ----------------------------------------------------------------------------
## 1) Load + clean conflict points
## ----------------------------------------------------------------------------
tb_log_section("Load conflict points")

cdf <- readxl::read_excel(TB_CONFLICT_XLSX)
tb_log(sprintf("xlsx rows = %d, cols = %d", nrow(cdf), ncol(cdf)))

## Locate coordinate columns (Long/Lat or x/y) — flexible
lon_col <- intersect(c("Long","long","Longitude","x","X"), names(cdf))[1]
lat_col <- intersect(c("Lat","lat","Latitude","y","Y"),   names(cdf))[1]
if (is.na(lon_col) || is.na(lat_col)) {
  tb_log("xlsx coordinate columns not found; falling back to ConflictPoints.txt", "WARN")
  tx <- read.table(TB_CONFLICT_TXT, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE)
  ## Coordinates in ConflictPoints.txt are (lon, lat) per peek (x=lon, y=lat)
  cdf$Long <- tx$x; cdf$Lat <- tx$y
  lon_col <- "Long"; lat_col <- "Lat"
}
tb_log(sprintf("lon=%s, lat=%s", lon_col, lat_col))

## Normalise Biogeographic_Region spelling
br_col <- intersect(c("Biogeographic_Region", "Biogeographic_region",
                      "BiogeographicRegion"), names(cdf))[1]
if (is.na(br_col)) stop("Biogeographic_Region column not found")
br_clean <- function(x) {
  x <- gsub("\\s+", "", x); x <- gsub("-", "", x)
  out <- toupper(x)
  ifelse(grepl("EUROSIB", out),  "Euro-Siberian",
  ifelse(grepl("IRANOTUR", out), "Irano-Turanian",
  ifelse(grepl("MEDITERR", out), "Mediterranean", x)))
}
cdf$bioregion <- br_clean(cdf[[br_col]])

## Activity column — keep as-is, drop NAs
act_col <- intersect(c("Activity", "activity"), names(cdf))[1]
if (is.na(act_col)) stop("Activity column not found")
cdf$activity <- as.character(cdf[[act_col]])
cdf <- cdf[!is.na(cdf$activity) & nchar(cdf$activity) > 0, ]

tb_log("bioregion counts:"); print(table(cdf$bioregion, useNA = "ifany"))
tb_log("activity counts:");  print(table(cdf$activity,  useNA = "ifany"))

## Convert to sf, project to TB_CRS_PROJ
pts_wgs <- sf::st_as_sf(cdf, coords = c(lon_col, lat_col), crs = TB_CRS_WGS,
                       remove = FALSE)
pts <- sf::st_transform(pts_wgs, TB_CRS_PROJ)

## ----------------------------------------------------------------------------
## 2) PA overlay (any-of-12)
## ----------------------------------------------------------------------------
tb_log_section("PA overlay")

pa_combined_gpkg <- file.path(TB_OUT_VECTORS, "pa_combined.gpkg")
if (!file.exists(pa_combined_gpkg)) {
  tb_log("pa_combined.gpkg missing — run 18_pa_overlay first", "ERROR")
  tb_log_session(); quit(status = 1)
}
pa <- sf::st_read(pa_combined_gpkg, quiet = TRUE) |> sf::st_transform(TB_CRS_PROJ)
pa_union <- sf::st_union(sf::st_make_valid(pa))

pts$in_pa <- as.logical(lengths(sf::st_intersects(pts, pa_union)) > 0)

tb_log(sprintf("conflicts inside PA: %d / %d (%.1f%%)",
               sum(pts$in_pa), nrow(pts),
               100 * mean(pts$in_pa)))

## ----------------------------------------------------------------------------
## 3) Road distance overlay
## ----------------------------------------------------------------------------
tb_log_section("Road distance overlay")

roads_gpkg <- file.path(TB_OUT_VECTORS, "roads_paved.gpkg")
if (!file.exists(roads_gpkg)) {
  tb_log("roads_paved.gpkg missing — run 19_roads_overlay first", "ERROR")
  tb_log_session(); quit(status = 1)
}
roads <- sf::st_read(roads_gpkg, quiet = TRUE) |> sf::st_transform(TB_CRS_PROJ)
roads_union <- sf::st_union(roads)

## st_distance to MULTILINESTRING is per-point shortest distance — fast for ~478
d_road_m <- as.numeric(sf::st_distance(pts, roads_union))
pts$dist_road_m <- d_road_m

tb_log(sprintf("road distance: median = %.0f m | mean = %.0f m | min = %.0f | max = %.0f",
               median(d_road_m), mean(d_road_m), min(d_road_m), max(d_road_m)))

## ----------------------------------------------------------------------------
## 4) Save clean joined table
## ----------------------------------------------------------------------------
clean_df <- sf::st_drop_geometry(pts) |>
  dplyr::select(any_of(c(lon_col, lat_col,
                         "bioregion", "activity",
                         "in_pa", "dist_road_m")))
tb_save_table(clean_df, "21_conflict_clean")
tb_save_vector(pts[, c("bioregion","activity","in_pa","dist_road_m")],
               "conflict_overlay", subdir = NULL)

## ----------------------------------------------------------------------------
## 5) Cross-tabs
## ----------------------------------------------------------------------------
tb_log_section("Cross-tabs")

## (a) bioregion × activity
xt_br <- clean_df |>
  dplyr::count(bioregion, activity) |>
  dplyr::group_by(activity) |>
  dplyr::mutate(pct_within_activity = 100 * n / sum(n)) |>
  dplyr::ungroup()
tb_save_table(xt_br, "21_xtab_bioregion_activity")

## (b) PA in/out × activity
xt_pa <- clean_df |>
  dplyr::mutate(pa_status = ifelse(in_pa, "Inside PA", "Outside PA")) |>
  dplyr::count(pa_status, activity) |>
  dplyr::group_by(activity) |>
  dplyr::mutate(pct_within_activity = 100 * n / sum(n)) |>
  dplyr::ungroup()
tb_save_table(xt_pa, "21_xtab_pa_activity")

## (c) road distance summary per activity
xt_rd <- clean_df |>
  dplyr::group_by(activity) |>
  dplyr::summarise(
    n         = dplyr::n(),
    mean_m    = mean(dist_road_m,   na.rm = TRUE),
    median_m  = median(dist_road_m, na.rm = TRUE),
    q25_m     = quantile(dist_road_m, 0.25, na.rm = TRUE),
    q75_m     = quantile(dist_road_m, 0.75, na.rm = TRUE),
    p_within_500m  = mean(dist_road_m <= 500,  na.rm = TRUE) * 100,
    p_within_1000m = mean(dist_road_m <= 1000, na.rm = TRUE) * 100,
    .groups   = "drop") |>
  dplyr::arrange(median_m)
tb_save_table(xt_rd, "21_road_dist_by_activity")

## ----------------------------------------------------------------------------
## 6) FIGURES
## ----------------------------------------------------------------------------
tb_log_section("Figures")

## ---- fig21a: bioregion × activity heatmap -----------------------------------
br_levels <- c("Euro-Siberian", "Irano-Turanian", "Mediterranean")
xt_br_wide <- clean_df |>
  dplyr::count(bioregion, activity) |>
  tidyr::complete(bioregion = br_levels, activity, fill = list(n = 0))

p21a <- ggplot(xt_br_wide, aes(activity, bioregion, fill = n)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = n), color = "black", size = 4) +
  scale_fill_gradient(low = "#FCE4EC", high = "#9E2A2B", name = "Count") +
  scale_x_discrete(guide = guide_axis(angle = 25)) +
  labs(title    = "Conflict events: biogeographic regions x activity reason",
       subtitle = sprintf("n = %d", nrow(clean_df)),
       x = NULL, y = NULL) +
  theme_trbear_bar(base_size = 12) +
  theme(panel.grid.major = element_blank(),
        axis.text.x = element_text(face = "bold"))
tb_save_fig(p21a, "fig21a_bioregion_activity_heatmap",
            w = 14, h = 6, subdir = FIG_SUBDIR)

## ---- fig21b: PA in/out × activity (stacked bars) ----------------------------
## Stack order: the LARGER-proportion class (Outside PA, red) is at the BOTTOM
## of the bar; the smaller proportion (Inside PA, green) is at the top.  This
## matches the requested visual emphasis (most-frequent class anchored to base).
xt_pa_p <- xt_pa
xt_pa_p$pa_status <- factor(xt_pa_p$pa_status,
                            levels = c("Inside PA", "Outside PA"))
p21b <- ggplot(xt_pa_p, aes(activity, n, fill = pa_status)) +
  geom_col(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%d (%.0f%%)", n, pct_within_activity)),
            position = position_stack(vjust = 0.5),
            color = "white", size = 3.4) +
  scale_fill_manual(values = c("Outside PA" = "#9E2A2B",
                               "Inside PA"  = "#009E73"),
                    breaks = c("Outside PA", "Inside PA"),
                    name   = NULL) +
  scale_x_discrete(guide = guide_axis(angle = 25)) +
  labs(title    = "Conflict events by activity, split by protected areas",
       subtitle = sprintf(
         "%d of %d (%.1f%%) conflict points fall inside any of the 12 PA layers.",
         sum(clean_df$in_pa), nrow(clean_df), 100 * mean(clean_df$in_pa)),
       x = NULL, y = "Number of conflict incidents") +
  theme_trbear_bar(base_size = 12) +
  theme(axis.text.x = element_text(face = "bold"))
tb_save_fig(p21b, "fig21b_pa_activity_bars",
            w = 14, h = 7, subdir = FIG_SUBDIR)

## ---- fig21c: road distance boxplot per activity -----------------------------
clean_df$activity_f <- factor(
  clean_df$activity,
  levels = xt_rd$activity)   # ordered by median dist (asc)

## Linear km axis, capped at the 99th percentile for readability; outliers
## above the cap shown as a top-of-panel rug.  Avoids log10(0) = -Inf artefacts
## from points that sit directly on a road.
clean_df$dist_km     <- clean_df$dist_road_m / 1000
y_cap                <- as.numeric(quantile(clean_df$dist_km, 0.99, na.rm = TRUE))
clean_df$dist_km_cap <- pmin(clean_df$dist_km, y_cap)
clean_df$capped      <- clean_df$dist_km > y_cap

## More vivid palette for the box-plot (Tableau-10 ordered to match the
## ascending-median activity ordering).
PAL_BOX <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2",
              "#59A14F", "#EDC948", "#B07AA1")
names(PAL_BOX) <- levels(clean_df$activity_f)
p21c <- ggplot(clean_df, aes(activity_f, dist_km_cap)) +
  geom_boxplot(aes(fill = activity_f),
               width = 0.65, alpha = 0.9, outlier.size = 0.8,
               outlier.alpha = 0.5) +
  scale_fill_manual(values = PAL_BOX, guide = "none") +
  scale_x_discrete(guide = guide_axis(angle = 25)) +
  scale_y_continuous(labels = scales::label_comma(),
                     breaks = scales::pretty_breaks(6),
                     expand = expansion(mult = c(0.02, 0.08))) +
  geom_rug(data = subset(clean_df, capped),
           aes(x = activity_f, y = y_cap),
           sides = "r", color = "#9E2A2B", linewidth = 0.4, alpha = 0.7) +
  labs(title    = "Distance from conflict events to nearest paved road",
       subtitle = NULL,
       x = NULL, y = "Distance to nearest paved road (km)") +
  theme_trbear_bar(base_size = 12) +
  theme(axis.text.x = element_text(face = "bold"))
tb_save_fig(p21c, "fig21c_road_dist_box",
            w = 14, h = 7, subdir = FIG_SUBDIR)

tb_log_session()
tb_log("21_conflict_overlay DONE")
