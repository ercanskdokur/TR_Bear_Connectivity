## ============================================================================
## 05_accessible_area.R
## Project: TR_Bear_Connectivity — ENMTML pipeline
## Purpose: Produce the TR boundary mask shapefile ENMTML uses as
##   sp_accessible_area = c(method='MASK', filepath=TB_TR_MASK_SHP).
##   Restricts model calibration + projection to Turkey landmass (excludes sea).
##
##   The shapefile must be a polygon (sf POLYGON / MULTIPOLYGON) and ENMTML
##   reads it via sf/raster. We project it to TR Albers EA to match predictors.
## Inputs:
##   outputs/rds_files/02_target_grid_info.rds  (contains tr_polygon_proj)
## Outputs:
##   data/TR_mask/TR_mask.shp                   ← TB_TR_MASK_SHP
##   outputs/figures/05_mask/05_TR_mask.png
## ============================================================================

suppressPackageStartupMessages({
  library(sf); library(terra); library(ggplot2); library(rnaturalearth)
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
tb_log_init("05_accessible_area")
tb_pkg_versions(c("sf","terra","rnaturalearth"))

## ============================================================================
## 1. LOAD TR POLYGON
## ============================================================================
tb_log_section("1. LOAD TR POLYGON")
rds_path <- file.path(TB_OUT_RDS, "02_target_grid_info.rds")
if (file.exists(rds_path)) {
  tr_proj <- readRDS(rds_path)$tr_polygon_proj
  tb_log("loaded TR polygon from 02_target_grid_info.rds")
} else {
  tb_log("02 RDS not found; fetching from rnaturalearth", "WARN")
  tr_wgs <- tryCatch(
    rnaturalearth::ne_countries(country = "Turkey", scale = "large", returnclass = "sf"),
    error = function(e)
      rnaturalearth::ne_countries(country = "Turkey", scale = "medium", returnclass = "sf")
  )
  tr_proj <- sf::st_transform(tr_wgs, TB_CRS_PROJ)
}

tb_log(sprintf("TR polygon: n_features=%d, bbox(m)=%s",
               nrow(tr_proj), paste(round(sf::st_bbox(tr_proj)), collapse=",")))

## ============================================================================
## 2. CLEAN + SIMPLIFY (optional)
## ============================================================================
tb_log_section("2. CLEAN")
## ensure valid geometries
tr_proj <- sf::st_make_valid(tr_proj)
## cast to MULTIPOLYGON for ENMTML compatibility
tr_proj <- sf::st_cast(tr_proj, "MULTIPOLYGON", warn = FALSE)
## keep only essential attributes (avoid encoding warnings on Turkish characters)
tr_proj <- tr_proj["geometry"]
tr_proj$id <- 1L
tr_proj$name <- "Turkey"

## ============================================================================
## 3. WRITE SHAPEFILE
## ============================================================================
tb_log_section("3. WRITE TR_mask.shp")
dir.create(dirname(TB_TR_MASK_SHP), recursive = TRUE, showWarnings = FALSE)
sf::st_write(tr_proj, TB_TR_MASK_SHP, delete_dsn = TRUE, quiet = TRUE)
tb_log(sprintf("wrote: %s", TB_TR_MASK_SHP))

## sanity: re-read to verify
tr_back <- sf::st_read(TB_TR_MASK_SHP, quiet = TRUE)
tb_log(sprintf("verify read: %d features, CRS=%s",
               nrow(tr_back), sf::st_crs(tr_back)$proj4string))

## ============================================================================
## 4. FIGURE
## ============================================================================
tb_log_section("4. FIGURE")
basemap <- tb_basemap_world(TB_CRS_PROJ, "medium")
bbox <- sf::st_bbox(tr_proj)
p <- ggplot() +
  geom_sf(data = basemap, fill = TB_FILL_LAND, color = TB_COLOR_LAND, linewidth = 0.3) +
  geom_sf(data = tr_proj, fill = "#56B4E9", color = TB_COLOR_FRAME, linewidth = 0.5, alpha = 0.5) +
  coord_sf(xlim = c(bbox$xmin - 80000, bbox$xmax + 80000),
           ylim = c(bbox$ymin - 80000, bbox$ymax + 80000),
           expand = FALSE, datum = sf::st_crs(4326)) +
  theme_trbear() + tb_map_decorations() +
  labs(title = "TR mask (ENMTML accessible area)",
       subtitle = sprintf("MASK method  |  %s", basename(TB_TR_MASK_SHP)),
       caption = "Restricts model calibration + projection to Turkey landmass")
tb_save_fig(p, "05_TR_mask", w = 14, h = 8, subdir = "05_mask")

tb_log_session()
tb_log("05_accessible_area DONE")
