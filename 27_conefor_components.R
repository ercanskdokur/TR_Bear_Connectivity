## ============================================================================
## 27_conefor_components.R
## Project: TR_Bear_Connectivity
## Purpose: Decompose per-patch connectivity importance into its intra-patch,
##   flux and connector fractions (Saura & Rubio 2010) across dispersal scales,
##   separating cores that matter for their own habitat area from those that
##   matter as topological stepping stones (high connector fraction).
##
##   For each of the six dispersal distances used throughout (50, 100, 150, 200,
##   300, 400 km) we compute dPC and dIIC and their three fractions, then plot
##   the top-6 cores (ranked by mean dPC across distances) as one panel per core
##   with four lines each — total, intra, flux, connector — and likewise for
##   dIIC. A companion map shows the top-6 cores at d = 100 km coloured by their
##   dominant fraction.
##
##   Decomposition (Saura & Rubio 2010, Ecography):
##     dPC_k        = dPCintra_k + dPCflux_k + dPCconnector_k
##     dPCintra_k   = 100 · a_k² / (A_L² · PC)                              %
##     dPCflux_k    = 100 · 2 · a_k · Σ_{i≠k} a_i · P*_ki(full) / (A_L² · PC)
##     dPCconnector_k = dPC_k − dPCintra_k − dPCflux_k                      %
##   (and analogous IIC formulas with P* replaced by 1/(1+nl))
##
## Inputs:
##   tables/35_costdist_dpc_dii.csv   (patch coords + area, present, all d)
##   tables/35_costdist_indices.csv   (PC, IIC headline values for normalisation)
##
## Outputs (tables/):
##   27_conefor_components.csv        (one row per patch × d_km)
##   27_top6_dPC.csv                  (wide, 6 patches × 6 distances × 4 cols)
##   27_top6_dIIC.csv
## Outputs (figures/27_components/):
##   fig27a_dpc_top6_distance.png     (dPC total/intra/flux/connector, 6 panels)
##   fig27b_diic_top6_distance.png    (dIIC total/intra/flux/connector, 6 panels)
##   fig27c_components_map.png         (TR map, top-6 patches at d = 100 km
##                                       coloured by dominant component)
## ============================================================================

suppressPackageStartupMessages({
  library(igraph); library(terra); library(sf); library(dplyr); library(tidyr)
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
tb_log_init("27_conefor_components")

FIG_SUBDIR     <- "27_components"
DIST_KM        <- c(50, 100, 150, 200, 300, 400)
PROB_AT_THRESH <- 0.5
TOP_N          <- 6   # top 6 cores by mean dPC

## ---- Inputs ----------------------------------------------------------------
tb_log_section("Inputs")

dpc_df <- read.csv(file.path(TB_OUT_TABLES, "35_costdist_dpc_dii.csv"))  # cost-distance
idx_df <- read.csv(file.path(TB_OUT_TABLES, "35_costdist_indices.csv"))
## cost-distance matrix + threshold calibration from script 35
.cm   <- readRDS(file.path(TB_OUT_RDS, "35_cost_matrices.rds"))
R_EFF <- read.csv(file.path(TB_OUT_TABLES, "35_costdist_calibration.csv"))$R_eff[1]

## Canonical core labels (C01..C93) from 34_core_crosswalk.R; fall back to raw id
.cw_file <- file.path(TB_OUT_TABLES, "34_core_crosswalk.csv")
core_cw <- if (file.exists(.cw_file)) {
  read.csv(.cw_file)[, c("patch_id", "core_id")]
} else {
  tb_log("34_core_crosswalk.csv missing — using raw patch ids", "WARN")
  data.frame(patch_id = unique(dpc_df$patch_id),
             core_id  = paste0("id", unique(dpc_df$patch_id)))
}
core_of <- function(pid) {
  out <- core_cw$core_id[match(pid, core_cw$patch_id)]
  ifelse(is.na(out), paste0("id", pid), out)
}

## present-scenario patches are identical across distances (same partition):
ps_present <- dpc_df[dpc_df$scenario == "present" & dpc_df$d_km == 100,
                      c("patch_id","x","y","area_km2")]
stopifnot(nrow(ps_present) > 1)
coords <- as.matrix(ps_present[, c("x","y")])
area   <- ps_present$area_km2
n      <- nrow(ps_present)

tr_shp <- file.path(TB_DATA_ROOT, "TR_mask", "TR_mask.shp")
tr_sf  <- sf::st_read(tr_shp, quiet = TRUE) |> sf::st_transform(TB_CRS_PROJ)
AL_km2 <- as.numeric(sum(sf::st_area(tr_sf))) / 1e6
tb_log(sprintf("A_L = %.0f km² | n_patches = %d", AL_km2, n))

## cost-distance (cost-km) matrix for the present cores, reordered to ps_present
.mo      <- match(ps_present$patch_id, .cm[["present"]]$tab$patch_id)
D_km_mat <- .cm[["present"]]$D[.mo, .mo]

## ---- One-distance decomposition --------------------------------------------
.compute_components_at_d <- function(d_km) {
  d_eff <- d_km * R_EFF                       # effective (cost) dispersal threshold
  alpha <- -log(PROB_AT_THRESH) / d_eff
  W <- -log(pmax(exp(-alpha * D_km_mat), .Machine$double.eps))
  diag(W) <- 0
  g_pc <- igraph::graph_from_adjacency_matrix(W, mode = "undirected",
                                                weighted = TRUE)
  Pstar_full <- exp(-igraph::distances(g_pc, weights = igraph::E(g_pc)$weight))
  PC_full <- sum(outer(area, area) * Pstar_full) / AL_km2^2

  Adj   <- (D_km_mat <= d_eff) & (D_km_mat > 0)
  g_iic <- igraph::graph_from_adjacency_matrix(Adj, mode = "undirected")
  nl_full <- igraph::distances(g_iic)
  IIC_kernel_full <- 1 / (1 + nl_full)
  IIC_full <- sum(outer(area, area) * IIC_kernel_full, na.rm = TRUE) / AL_km2^2

  denom_PC  <- AL_km2^2 * PC_full
  denom_IIC <- AL_km2^2 * IIC_full

  dPC_total <- dPCintra <- dPCflux <- numeric(n)
  dIIC_total <- dIICintra <- dIICflux <- numeric(n)

  for (k in seq_len(n)) {
    a_k <- area[k]; oth <- seq_len(n)[-k]
    g_pc_k  <- igraph::delete_vertices(g_pc, k)
    L_k     <- igraph::distances(g_pc_k, weights = igraph::E(g_pc_k)$weight)
    Pstar_k <- exp(-L_k)
    PC_minus_k <- sum(outer(area[oth], area[oth]) * Pstar_k) / AL_km2^2
    dPC_total[k] <- 100 * (PC_full - PC_minus_k) / PC_full
    dPCintra[k]  <- 100 * (a_k^2)                                 / denom_PC
    dPCflux[k]   <- 100 * 2 * a_k * sum(area[oth] * Pstar_full[k, oth]) / denom_PC

    g_iic_k <- igraph::delete_vertices(g_iic, k)
    nl_k    <- igraph::distances(g_iic_k)
    IIC_minus_k <- sum(outer(area[oth], area[oth]) / (1 + nl_k), na.rm = TRUE) /
                     AL_km2^2
    dIIC_total[k] <- 100 * (IIC_full - IIC_minus_k) / IIC_full
    dIICintra[k]  <- 100 * (a_k^2)                                / denom_IIC
    dIICflux[k]   <- 100 * 2 * a_k * sum(area[oth] * IIC_kernel_full[k, oth]) /
                       denom_IIC
  }
  dPCconnector  <- dPC_total  - dPCintra  - dPCflux
  dIICconnector <- dIIC_total - dIICintra - dIICflux

  data.frame(
    d_km          = d_km,
    patch_id      = ps_present$patch_id,
    x             = ps_present$x, y = ps_present$y,
    area_km2      = ps_present$area_km2,
    dPC_total     = dPC_total,
    dPCintra      = dPCintra,
    dPCflux       = dPCflux,
    dPCconnector  = dPCconnector,
    dIIC_total    = dIIC_total,
    dIICintra     = dIICintra,
    dIICflux      = dIICflux,
    dIICconnector = dIICconnector)
}

tb_log_section("Compute components at every dispersal distance")
all_components <- do.call(rbind, lapply(DIST_KM, function(d) {
  tb_log(sprintf("d = %d km", d))
  .compute_components_at_d(d)
}))
tb_save_table(all_components, "27_conefor_components")

## ---- Rank top patches by MEAN dPC across distances -------------------------
patch_mean <- all_components |>
  dplyr::group_by(patch_id) |>
  dplyr::summarise(
    area_km2     = mean(area_km2),
    x            = mean(x), y = mean(y),
    mean_dPC     = mean(dPC_total),
    mean_dIIC    = mean(dIIC_total),
    .groups      = "drop")

top_pc  <- patch_mean |> dplyr::arrange(dplyr::desc(mean_dPC))  |>
  dplyr::slice_head(n = TOP_N) |> dplyr::mutate(rank = dplyr::row_number(),
                                                core_id = core_of(patch_id))
## ---- fig27b cores selected AT THE FOCAL DISTANCE (d = 100), not by mean ------
## Mean-across-distance ranking buries the focal-distance connector stepping
## stones (e.g. C53 has ~0 dIIC at d=100 but a connector blip at d=200 that
## inflates its mean). We instead show the 3 dominant cores plus the 3 genuine
## connector stepping stones (small, dIIC >> dPC) at d = 100, matching the text
## and Fig. 4 (C40/C39/C37).
d100c     <- all_components |> dplyr::filter(d_km == 100)
.dom_iic  <- d100c |> dplyr::arrange(dplyr::desc(dIIC_total)) |>
  dplyr::slice_head(n = 3) |> dplyr::pull(patch_id)
.step_iic <- d100c |> dplyr::filter(area_km2 < 500, dIIC_total > dPC_total) |>
  dplyr::arrange(dplyr::desc(dIIC_total)) |>
  dplyr::slice_head(n = 3) |> dplyr::pull(patch_id)
.sel_iic  <- c(.dom_iic, .step_iic)
top_iic <- patch_mean |> dplyr::filter(patch_id %in% .sel_iic) |>
  dplyr::arrange(match(patch_id, .sel_iic)) |>
  dplyr::mutate(rank = dplyr::row_number(), core_id = core_of(patch_id))
tb_save_table(top_pc,  "27_top6_dPC")
tb_save_table(top_iic, "27_top6_dIIC")

## ---- WGS84 centroid lon/lat for panel titles -------------------------------
.wgs_centroids <- function(df) {
  s <- sf::st_as_sf(df, coords = c("x","y"), crs = TB_CRS_PROJ) |>
    sf::st_transform(TB_CRS_WGS)
  cc <- sf::st_coordinates(s)
  df$lon <- cc[, 1]; df$lat <- cc[, 2]
  df
}
top_pc  <- .wgs_centroids(top_pc)
top_iic <- .wgs_centroids(top_iic)

## ---- Long-format for plotting ----------------------------------------------
.make_long <- function(top_df, prefix) {
  base <- all_components |>
    dplyr::filter(patch_id %in% top_df$patch_id) |>
    dplyr::left_join(top_df |> dplyr::select(patch_id, rank, core_id, lon, lat),
                     by = "patch_id")
  ## Build panel-label factor levels in ascending rank order so panels appear
  ## as Core #1 .. Core #6 (left-to-right, top-to-bottom), not in patch-id
  ## order.  (Earlier version used unique() on data order, which produced
  ## #2, #5, #1 in the first row — that was the bug.)
  panel_order <- top_df |>
    dplyr::arrange(rank) |>
    dplyr::transmute(lbl = sprintf("%s (%.0f km²)\n(%.2f, %.2f)",
                                    core_id, area_km2, lon, lat)) |>
    dplyr::pull(lbl)
  long <- base |>
    dplyr::transmute(
      patch_id, rank, core_id, lon, lat, area_km2, d_km,
      total     = .data[[paste0(prefix, "_total")]],
      intra     = .data[[paste0(prefix, "intra")]],
      flux      = .data[[paste0(prefix, "flux")]],
      connector = .data[[paste0(prefix, "connector")]]) |>
    tidyr::pivot_longer(c("total","intra","flux","connector"),
                        names_to = "component", values_to = "value") |>
    dplyr::mutate(
      component = factor(component,
                          levels = c("total", "intra", "flux", "connector"),
                          labels = c(paste0(prefix, " (total)"),
                                      paste0(prefix, "intra"),
                                      paste0(prefix, "flux"),
                                      paste0(prefix, "connector"))),
      panel_lbl = factor(sprintf("%s (%.0f km²)\n(%.2f, %.2f)",
                                    core_id, area_km2, lon, lat),
                          levels = panel_order))
  long
}

long_pc  <- .make_long(top_pc,  "dPC")
long_iic <- .make_long(top_iic, "dIIC")

## palette (4 lines per panel: total / intra / flux / connector)
PAL_COMP4 <- c(
  "dPC (total)"   = "#000000",
  "dPCintra"      = "#0072B2",
  "dPCflux"       = "#E69F00",
  "dPCconnector"  = "#009E73",
  "dIIC (total)"  = "#000000",
  "dIICintra"     = "#0072B2",
  "dIICflux"      = "#E69F00",
  "dIICconnector" = "#009E73")

.make_component_panel <- function(long_df, ytitle, fig_title) {
  ggplot(long_df, aes(d_km, value, color = component, group = component)) +
    geom_line(aes(linewidth = component == levels(component)[1])) +
    geom_point(size = 2.4) +
    scale_color_manual(values = PAL_COMP4[levels(long_df$component)],
                       name = NULL) +
    scale_linewidth_manual(values = c("TRUE" = 1.2, "FALSE" = 0.7),
                            guide = "none") +
    facet_wrap(~ panel_lbl, ncol = 3, scales = "free_y") +
    scale_x_continuous(breaks = DIST_KM,
                       labels = paste0(DIST_KM, " km")) +
    labs(title = fig_title,
         x = "Dispersal distance threshold (km)",
         y = ytitle) +
    theme_trbear_bar(base_size = 11) +
    theme(axis.text.x  = element_text(angle = 25, hjust = 1, size = 9),
          strip.text   = element_text(face = "bold", size = 10),
          panel.spacing = unit(0.8, "lines"),
          legend.position = "top")
}

tb_log_section("Fig27a / Fig27b — dPC / dIIC component decomposition across distances")
fig27a <- .make_component_panel(
  long_pc, "dPC (% of PC)",
  sprintf("Mean values of dPC and its three fractions (intra, flux, connector) — top-%d core habitats × %d dispersal scenarios",
           TOP_N, length(DIST_KM)))
tb_save_fig(fig27a, "fig27a_dpc_top6_distance", w = 16, h = 9,
            subdir = FIG_SUBDIR)

fig27b <- .make_component_panel(
  long_iic, "dIIC (% of IIC)",
  sprintf("Mean values of dIIC and its three fractions (intra, flux, connector) — top-%d core habitats × %d dispersal scenarios",
          TOP_N, length(DIST_KM)))
tb_save_fig(fig27b, "fig27b_diic_top6_distance", w = 16, h = 9,
            subdir = FIG_SUBDIR)

## ---- Fig27c: TR map of top-6 patches at d=100 km, coloured by dominant -----
tb_log_section("Fig27c map of top-6 patches (d=100 km)")

at_100 <- all_components |>
  dplyr::filter(d_km == 100,
                patch_id %in% top_pc$patch_id) |>
  dplyr::left_join(top_pc |> dplyr::select(patch_id, rank, core_id, lon, lat),
                   by = "patch_id")
dom_label <- function(intra, flux, conn) {
  vals <- c(intra = intra, flux = flux, connector = conn)
  names(which.max(vals))
}
at_100$dominant <- mapply(dom_label, at_100$dPCintra, at_100$dPCflux,
                           at_100$dPCconnector)
top_sf <- sf::st_as_sf(at_100, coords = c("x","y"), crs = TB_CRS_PROJ)

bear_tif <- file.path(TB_OUT_HS_BINARY, "present_wmean.tif")
bear_r <- if (file.exists(bear_tif)) {
  br <- terra::rast(bear_tif)
  br <- terra::mask(br, terra::vect(tr_sf))
  b2 <- terra::as.factor(terra::as.int(br))
  levels(b2) <- data.frame(id = c(0, 1), class = c("Unsuitable", "Suitable"))
  names(b2) <- "class"; b2
} else NULL

world_sf <- tryCatch(
  rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") |>
    sf::st_transform(TB_CRS_PROJ), error = function(e) NULL)
tr_bb <- sf::st_bbox(tr_sf); pad <- 30000
xl <- c(tr_bb["xmin"] - pad, tr_bb["xmax"] + pad)
yl <- c(tr_bb["ymin"] - pad, tr_bb["ymax"] + pad)

PAL_DOM <- c("intra" = "#0072B2", "flux" = "#E69F00", "connector" = "#009E73")
p27c <- ggplot()
if (!is.null(world_sf)) p27c <- p27c + geom_sf(data = world_sf, fill = "#E8E8E8",
                                                color = "#7C8A93", linewidth = 0.3)
if (!is.null(bear_r))
  p27c <- p27c + tidyterra::geom_spatraster(data = bear_r, na.rm = TRUE) +
    scale_fill_manual(values = TB_PAL_BINARY, na.translate = FALSE,
                       name = "Habitat")
p27c <- p27c +
  geom_sf(data = tr_sf, fill = NA, color = TB_COLOR_FRAME, linewidth = 0.4) +
  geom_sf(data = top_sf, aes(color = dominant, size = dPC_total),
          alpha = 0.95, stroke = 0.8) +
  scale_color_manual(values = PAL_DOM, name = "Dominant\ncomponent (dPC)",
                      labels = c("intra"     = "Intra (within-patch)",
                                  "flux"      = "Flux (source/sink)",
                                  "connector" = "Connector (stepping-stone)")) +
  scale_size_continuous(name = "dPC at d=100 km (%)", range = c(4, 11)) +
  geom_sf_text(data = top_sf, aes(label = core_id),
               size = 2.9, color = "white", fontface = "bold") +
  coord_sf(xlim = xl, ylim = yl, datum = sf::st_crs(4326), expand = FALSE) +
  tb_map_decorations() +
  labs(title    = "Top-6 priority core habitats (present, d = 100 km)",
       subtitle = "Cores labelled C01..C93 by present dPC importance (d=100 km). Colour = which dPC fraction dominates at d=100 km.") +
  theme_trbear(base_size = 11)
tb_save_fig(p27c, "fig27c_components_map", w = 14, h = 8, subdir = FIG_SUBDIR)

tb_log_session()
tb_log("27_conefor_components DONE")
