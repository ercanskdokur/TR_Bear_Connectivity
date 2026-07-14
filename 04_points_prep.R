## ============================================================================
## 04_points_prep.R
## Project: TR_Bear_Connectivity — ENMTML pipeline
## Purpose: Prepare presence + conflict points for ENMTML.
##   (1) Reproject presence to TR Albers EA.
##   (2) Drop points falling on NA in predictor mask.
##   (3) Write ENMTML occurrence TXT (tab-separated, columns: species, x, y).
##   (4) Same for conflict points (separate ENMTML run later).
##   (5) Normalize conflict typology and save full conflict GPKG with metadata.
##   (6) QC figures.
## ENMTML expects: tab-delimited TXT with columns named exactly as passed to sp/x/y.
##   We use: species, x, y. Coordinates in the SAME CRS as predictors (AEA, meters).
## Inputs:
##   data/points/PresencePoints.txt
##   data/points/conflict_table_df.xlsx
##   outputs/rasters/target_grid.tif
##   data/predictors_enmtml/present/Bio01.tif (anchor)
## Outputs:
##   data/points/occ_enmtml.txt              ← presence for SDM
##   data/points/occ_conflict_enmtml.txt     ← conflict for separate SDM
##   outputs/vectors/presence_clean.gpkg
##   outputs/vectors/conflict_points.gpkg
##   outputs/tables/04_points_summary.csv
##   outputs/tables/04_conflict_by_activity.csv
##   outputs/figures/04_points/04_presence_map.png
##   outputs/figures/04_points/04_conflict_by_activity.png
## ============================================================================

suppressPackageStartupMessages({
  library(terra); library(sf); library(dplyr); library(ggplot2); library(readxl)
  library(rnaturalearth)
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
tb_log_init("04_points_prep")
tb_pkg_versions(c("terra","sf","readxl"))

## ============================================================================
## 1. PRESENCE — read, project, NA-drop, write ENMTML TXT
## ============================================================================
tb_log_section("1. PRESENCE")
tb_tic()

pres <- read.table(TB_PRESENCE_TXT, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, quote = "\"")
tb_log(sprintf("raw presence n = %d", nrow(pres)))

## column name detection (handle x/y or lon/lat or Long/Lat)
xcol <- intersect(c("x","X","lon","Long","Longitude","longitude"), names(pres))[1]
ycol <- intersect(c("y","Y","lat","Lat","Latitude","latitude"), names(pres))[1]
if (is.na(xcol) || is.na(ycol))
  stop("could not detect coordinate columns in PresencePoints.txt — found cols: ",
       paste(names(pres), collapse=","))
tb_log(sprintf("detected coord cols: %s, %s", xcol, ycol))

pres_sf <- sf::st_as_sf(pres, coords = c(xcol, ycol), crs = TB_CRS_WGS, remove = FALSE)
pres_sf <- sf::st_transform(pres_sf, TB_CRS_PROJ)

## drop points on NA in anchor predictor
anchor <- terra::rast(file.path(TB_PRED_ENMTML_PRESENT, "Bio01.tif"))
vals <- terra::extract(anchor, terra::vect(pres_sf))[, 2]
n0 <- nrow(pres_sf)
pres_sf <- pres_sf[!is.na(vals), , drop = FALSE]
tb_log(sprintf("dropped %d presence points on NA cells; final n = %d",
               n0 - nrow(pres_sf), nrow(pres_sf)))

tb_save_vector(pres_sf, "presence_clean")

## --- write ENMTML occurrence TXT ---
coords <- sf::st_coordinates(pres_sf)
occ_df <- data.frame(species = "Ursus_arctos",
                     x       = coords[, "X"],
                     y       = coords[, "Y"],
                     stringsAsFactors = FALSE)
write.table(occ_df, TB_OCC_ENMTML, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)
tb_log(sprintf("wrote ENMTML occurrence TXT: %s (%d rows)",
               TB_OCC_ENMTML, nrow(occ_df)))
tb_toc("presence")

## ============================================================================
## 2. CONFLICT — read xlsx, normalize, project, write ENMTML TXT
## ============================================================================
tb_log_section("2. CONFLICT POINTS + TYPOLOGY")
tb_tic()

conf <- readxl::read_excel(TB_CONFLICT_XLSX)
tb_log(sprintf("raw conflict n = %d", nrow(conf)))

## normalize Biogeographic_Region spelling
br <- conf$Biogeographic_Region
br <- gsub("^Eurosiberian$",  "Euro-Siberian",  br)
br <- gsub("^IranoTuranian$", "Irano-Turanian", br)
conf$Biogeographic_Region <- br

nc0 <- nrow(conf)
conf <- conf[!is.na(conf$Long) & !is.na(conf$Lat), ]
if (nrow(conf) < nc0)
  tb_log(sprintf("dropped %d conflict rows missing coords", nc0 - nrow(conf)))

conf_sf <- sf::st_as_sf(conf, coords = c("Long","Lat"), crs = TB_CRS_WGS, remove = FALSE)
conf_sf <- sf::st_transform(conf_sf, TB_CRS_PROJ)

vals_c <- terra::extract(anchor, terra::vect(conf_sf))[, 2]
nc_before <- nrow(conf_sf)
conf_sf <- conf_sf[!is.na(vals_c), , drop = FALSE]
tb_log(sprintf("dropped %d conflict points on NA cells; final n = %d",
               nc_before - nrow(conf_sf), nrow(conf_sf)))

tb_save_vector(conf_sf, "conflict_points")

## ENMTML conflict TXT
cc <- sf::st_coordinates(conf_sf)
occ_conf <- data.frame(species = "Ursus_arctos_conflict",
                       x = cc[, "X"], y = cc[, "Y"],
                       stringsAsFactors = FALSE)
write.table(occ_conf, TB_OCC_CONFLICT_ENMTML, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)
tb_log(sprintf("wrote ENMTML conflict TXT: %s (%d rows)",
               TB_OCC_CONFLICT_ENMTML, nrow(occ_conf)))

## activity table
act_tab <- conf_sf |> sf::st_drop_geometry() |>
  count(Activity, name = "n") |> arrange(-n)
tb_save_table(act_tab, "04_conflict_by_activity")
tb_toc("conflict")

## ============================================================================
## 3. SUMMARY
## ============================================================================
tb_log_section("3. SUMMARY")
summary_df <- data.frame(
  metric = c("presence_raw","presence_after_NA_drop",
             "conflict_raw","conflict_after_NA_drop",
             "occ_enmtml_file","occ_conflict_enmtml_file"),
  value  = c(nrow(pres), nrow(pres_sf),
             nc0, nrow(conf_sf),
             TB_OCC_ENMTML, TB_OCC_CONFLICT_ENMTML),
  stringsAsFactors = FALSE
)
tb_save_table(summary_df, "04_points_summary")

## ============================================================================
## 4. FIGURES
## ============================================================================
tb_log_section("4. FIGURES")
tb_tic()
tr_proj <- readRDS(file.path(TB_OUT_RDS, "02_target_grid_info.rds"))$tr_polygon_proj
bbox <- sf::st_bbox(tr_proj)
basemap <- tb_basemap_world(TB_CRS_PROJ, "medium")

p1 <- ggplot() +
  geom_sf(data = basemap, fill = TB_FILL_LAND, color = TB_COLOR_LAND, linewidth = 0.3) +
  geom_sf(data = tr_proj, fill = "white", color = TB_COLOR_FRAME, linewidth = 0.4) +
  geom_sf(data = pres_sf, color = "#0072B2", size = 1.0, alpha = 0.7) +
  coord_sf(xlim = c(bbox$xmin - 60000, bbox$xmax + 60000),
           ylim = c(bbox$ymin - 60000, bbox$ymax + 60000),
           expand = FALSE, datum = sf::st_crs(4326)) +
  theme_trbear() + tb_map_decorations() +
  labs(title    = "Brown bear presence points",
       subtitle = sprintf("n=%d", nrow(pres_sf)))
tb_save_fig(p1, "04_presence_map", w = 14, h = 8, subdir = "04_points")

conf_df <- sf::st_drop_geometry(conf_sf)
conf_df$Activity <- factor(conf_df$Activity,
                           levels = names(sort(table(conf_df$Activity), decreasing = TRUE)))
p2 <- ggplot() +
  geom_sf(data = basemap, fill = TB_FILL_LAND, color = TB_COLOR_LAND, linewidth = 0.3) +
  geom_sf(data = tr_proj, fill = "white", color = TB_COLOR_FRAME, linewidth = 0.4) +
  geom_sf(data = conf_sf, aes(color = Activity), size = 1.4, alpha = 0.85) +
  scale_color_manual(values = TB_PAL_ACTIVITY, name = "Activity") +
  coord_sf(xlim = c(bbox$xmin - 60000, bbox$xmax + 60000),
           ylim = c(bbox$ymin - 60000, bbox$ymax + 60000),
           expand = FALSE, datum = sf::st_crs(4326)) +
  theme_trbear() + tb_map_decorations() +
  labs(title    = "Human–bear conflict points by activity",
       subtitle = sprintf("n=%d events (2017–2025)", nrow(conf_sf)))
tb_save_fig(p2, "04_conflict_by_activity", w = 14, h = 8, subdir = "04_points")
tb_toc("figures")

tb_save_rds(list(
  presence_sf = pres_sf, conflict_sf = conf_sf,
  occ_enmtml = TB_OCC_ENMTML, occ_conflict_enmtml = TB_OCC_CONFLICT_ENMTML,
  summary = summary_df
), "04_points_prep")

tb_log_session()
tb_log("04_points_prep DONE")
