#!/usr/bin/env Rscript
# ============================================================================
# E_source_independence.R
# TR_Bear_Connectivity: Presence / HBC-conflict source independence check (eval item 8)
#
#
# This script quantifies that coincidence at several distance thresholds and
# reports it explicitly as a spatial proxy, not a database-level provenance
#
# Outputs:
#   tables/E_source_independence_summary.csv   threshold_km -> n_coincident, pct_presence, pct_conflict
#   tables/E_source_independence_pairs.csv     nearest conflict point per presence point (distance_km)
#   figures/E_source_independence/fig_E_coincidence_map.png
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(sf); library(ggplot2)
})

.tb_find_paths_R <- function() {
  a <- commandArgs(trailingOnly = FALSE); f <- a[grepl("--file=", a)]
  if (length(f)) { d <- dirname(normalizePath(sub("--file=", "", f[1]), mustWork = FALSE))
    if (file.exists(file.path(d, "00_paths.R"))) return(d) }
  env_dir <- Sys.getenv("TB_SCRIPTS", unset = "")
  if (nzchar(env_dir) && file.exists(file.path(env_dir, "00_paths.R"))) return(env_dir)
  if (file.exists("00_paths.R")) return(getwd()); stop("Cannot find 00_paths.R")
}
setwd(.tb_find_paths_R()); source("00_paths.R"); source("00_helpers.R")
tb_log_init("E_source_independence")

FIG_SUBDIR <- "E_source_independence"
THRESH_KM  <- c(0.05, 1, 5, 10)   # ~50 m (same GPS fix), 1 km (thinning cell), 5 km, 10 km

## ---- load ------------------------------------------------------------------
pres <- read_delim(TB_PRESENCE_TXT, delim = "\t", show_col_types = FALSE)
conf <- read_delim(TB_CONFLICT_TXT, delim = "\t", show_col_types = FALSE)
tb_log(sprintf("presence n=%d, conflict n=%d (as used in the main pipeline)", nrow(pres), nrow(conf)))

pres_sf <- st_as_sf(pres, coords = c("x", "y"), crs = 4326) |> st_transform(TB_CRS_PROJ)
conf_sf <- st_as_sf(conf, coords = c("x", "y"), crs = 4326) |> st_transform(TB_CRS_PROJ)

## ---- nearest-neighbour distance, presence -> conflict ---------------------
nn <- st_nearest_feature(pres_sf, conf_sf)
d_m <- as.numeric(st_distance(pres_sf, conf_sf[nn, ], by_element = TRUE))
pairs <- pres |>
  mutate(nearest_conflict_idx = nn, distance_km = d_m / 1000) |>
  arrange(distance_km)
tb_save_table(pairs, "E_source_independence_pairs")

## ---- summary at multiple thresholds ----------------------------------------
summary_rows <- lapply(THRESH_KM, function(th) {
  n_coincident_from_presence <- sum(pairs$distance_km <= th)
  ## symmetric check: how many conflict points have a presence point within th
  nn2 <- st_nearest_feature(conf_sf, pres_sf)
  d2_m <- as.numeric(st_distance(conf_sf, pres_sf[nn2, ], by_element = TRUE))
  n_coincident_from_conflict <- sum(d2_m / 1000 <= th)
  data.frame(
    threshold_km = th,
    n_presence_near_conflict = n_coincident_from_presence,
    pct_presence_near_conflict = 100 * n_coincident_from_presence / nrow(pres),
    n_conflict_near_presence = n_coincident_from_conflict,
    pct_conflict_near_presence = 100 * n_coincident_from_conflict / nrow(conf))
}) |> bind_rows()
tb_save_table(summary_rows, "E_source_independence_summary")
tb_log("Coincidence summary:")
for (i in seq_len(nrow(summary_rows))) {
  r <- summary_rows[i, ]
  tb_log(sprintf("  <= %.2f km: %d/%d presence pts (%.1f%%), %d/%d conflict pts (%.1f%%)",
    r$threshold_km, r$n_presence_near_conflict, nrow(pres), r$pct_presence_near_conflict,
    r$n_conflict_near_presence, nrow(conf), r$pct_conflict_near_presence))
}

## ---- exact-coordinate duplicates (strongest signal of shared origin) ------
exact_dupe <- sum(pairs$distance_km < 0.001)  # < 1 m: effectively identical GPS fix
tb_log(sprintf("Exact-coordinate duplicates (< 1 m): %d of %d presence points", exact_dupe, nrow(pres)))

## ---- figure: where do the two datasets coincide? ---------------------------
tb_log_section("Figure")
tr_mask_sf <- sf::st_read(file.path(TB_DATA_ROOT, "TR_mask", "TR_mask.shp"), quiet = TRUE) |>
  sf::st_transform(TB_CRS_PROJ)
pres_plot <- pres_sf |> mutate(near_conflict_1km = pairs$distance_km <= 1)

p <- ggplot() +
  geom_sf(data = tr_mask_sf, fill = "#F2F2F2", color = "#BDBDBD", linewidth = 0.3) +
  geom_sf(data = conf_sf, color = "#D55E00", size = 1.1, alpha = 0.55, shape = 17) +
  geom_sf(data = pres_plot, aes(color = near_conflict_1km), size = 1.1, alpha = 0.7) +
  scale_color_manual(values = c(`FALSE` = "#0072B2", `TRUE` = "#CC79A7"),
                      labels = c(`FALSE` = "Presence, > 1 km from any conflict record",
                                 `TRUE`  = "Presence, \u2264 1 km from a conflict record"),
                      name = NULL) +
  labs(title = "Spatial coincidence of presence and HBC-conflict records",
       subtitle = sprintf("Triangles = conflict events (n=%d). %.1f%% of presence points (n=%d) fall within 1 km of a conflict record.",
                           nrow(conf), summary_rows$pct_presence_near_conflict[summary_rows$threshold_km == 1], nrow(pres))) +
  theme_trbear(base_size = 12) +
  theme(legend.position = "bottom")
tb_save_fig(p, "fig_E_coincidence_map", w = 9, h = 8, subdir = FIG_SUBDIR)

tb_log_session(); tb_log("E_source_independence DONE")
