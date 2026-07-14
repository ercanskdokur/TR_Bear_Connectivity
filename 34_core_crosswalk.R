## ============================================================================
## 34_core_crosswalk.R
## Project: TR_Bear_Connectivity
## Purpose: Assign clean, consecutive, manuscript-facing labels to the 93 present
##   source patches. terra::patches() labels every connected component in the
##   full binary raster (including sub-threshold specks), so the retained patches
##   keep arbitrary ids (236, 31, 1098, ...). Here we map each retained patch to
##   a consecutive identity C01..C93, ranked by present dPC at the focal
##   dispersal distance (d = 100 km): C01 = most connectivity-important core.
##
## This crosswalk is the single source of truth read by scripts 27/28/30 so that
## the same core is named identically across every figure and table.
##
## Output (tables/):
##   34_core_crosswalk.csv   patch_id, core_id, core_rank, area_km2,
##                            dPC_d100, dIIC_d100, x, y, lon, lat
## ============================================================================

suppressPackageStartupMessages({ library(sf); library(dplyr) })

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
tb_log_init("34_core_crosswalk")

FOCAL_D <- 100
dpc <- read.csv(file.path(TB_OUT_TABLES, "17_conefor_dpc_dii.csv"))
cw <- dpc |>
  dplyr::filter(scenario == "present", d_km == FOCAL_D) |>
  dplyr::arrange(dplyr::desc(dPC)) |>
  dplyr::mutate(core_rank = dplyr::row_number(),
                core_id   = sprintf("C%02d", core_rank)) |>
  dplyr::select(patch_id, core_id, core_rank,
                area_km2, dPC_d100 = dPC, dIIC_d100 = dIIC, x, y)

ll <- sf::st_as_sf(cw, coords = c("x","y"), crs = TB_CRS_PROJ) |>
  sf::st_transform(TB_CRS_WGS) |> sf::st_coordinates()
cw$lon <- ll[, 1]; cw$lat <- ll[, 2]
cw$x <- dpc$x[match(cw$patch_id, dpc$patch_id)]   # keep projected coords too
cw$y <- dpc$y[match(cw$patch_id, dpc$patch_id)]

tb_save_table(cw, "34_core_crosswalk")
tb_log(sprintf("crosswalk: %d cores | C01=patch %d (%.0f km², dPC=%.1f%%)",
               nrow(cw), cw$patch_id[1], cw$area_km2[1], cw$dPC_d100[1]))
tb_log_session()
tb_log("34_core_crosswalk DONE")
