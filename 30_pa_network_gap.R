## ============================================================================
## 30_pa_network_gap.R
## Project: TR_Bear_Connectivity
## Purpose: Is Türkiye's protected-area (PA) system a CONNECTED network, or a
##   set of isolated islands? And which UNPROTECTED cores are the critical glue
##   that should become new PAs? This directly answers the explicit call of
##   Sıkdokur et al. (2025) to "design new PAs that enhance habitat connectivity".
##
## Method (present scenario, focal d = 100 km):
##   - Re-extract present cores (terra::patches on present_wmean binary; ids
##     match 17_conefor_dpc_dii.csv).
##   - Classify each core by % of its area inside the combined PA layer
##     (pa_combined.gpkg from script 18; rebuilt from PAs.gpkg if absent):
##         Protected   : >= 50% inside PA
##         Partial     : 5–50% inside PA
##         Unprotected : < 5% inside PA
##   - PA-subnetwork: restrict the graph to Protected cores, count connected
##     components at d=100, and compute PC of the PA-only network relative to
##     the full network (how much national connectivity the PA system captures).
##   - Conservation priority: rank UNPROTECTED cores by full-network dPC and by
##     the connector fraction (dPCconnector / dIICconnector from script 27) —
##     unprotected stepping-stones are the highest-value new-PA candidates.
##
## Outputs (tables/):
##   30_core_protection.csv       patch_id, area, dPC, dIIC, pct_in_pa, status
##   30_pa_network_summary.csv    components, PA-network PC share, etc.
##   30_unprotected_priority.csv  ranked unprotected cores (new-PA candidates)
## Figures (figures/30_pa_gap/):
##   fig30a_pa_gap_network.png    map: cores by status (size=dPC) + d=100 links
##   fig30b_unprotected_priority.png  bars of top unprotected cores (dPC/dIIC)
## ============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("igraph", quietly = TRUE))
    install.packages("igraph", repos = "https://cloud.r-project.org")
  library(terra); library(sf); library(dplyr); library(tidyr)
  library(ggplot2); library(igraph); library(patchwork)
  library(rnaturalearth); library(tidyterra)
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
tb_log_init("30_pa_network_gap")

FIG_SUBDIR <- "30_pa_gap"
FOCAL_D    <- 100
PROB_AT_THRESH <- 0.5
PROT_HI <- 0.50   # >= -> Protected
PROT_LO <- 0.05   # < -> Unprotected; between -> Partial

PA_LAYERS <- c(
  "hassas_sukutle", "millipark", "MUHAZAORM", "OZELCEVREKORUMA",
  "REKR_KENTORMANI", "REKR_MESIREALAN",
  "sulak_MahOnHaSuAl", "sulak_Ramsar", "sulak_UlnHaSuAl",
  "tabiat_koruma_alani", "tabiat_parki", "YABANHAYATIGELSAH")

## ----------------------------------------------------------------------------
## 1) Present cores + dPC/dIIC + connector fractions
## ----------------------------------------------------------------------------
tb_log_section("Inputs")
dpc <- read.csv(file.path(TB_OUT_TABLES, "35_costdist_dpc_dii.csv"))  # cost-distance
ps  <- dpc |> dplyr::filter(scenario == "present", d_km == FOCAL_D) |>
  dplyr::select(patch_id, x, y, area_km2, dPC, dIIC)
n <- nrow(ps)
## cost-distance matrix + calibration from script 35
.cm   <- readRDS(file.path(TB_OUT_RDS, "35_cost_matrices.rds"))
R_EFF <- read.csv(file.path(TB_OUT_TABLES, "35_costdist_calibration.csv"))$R_eff[1]

## Canonical core labels (C01..C93) from 34_core_crosswalk.R; fall back to raw id
.cw_file <- file.path(TB_OUT_TABLES, "34_core_crosswalk.csv")
core_cw <- if (file.exists(.cw_file)) read.csv(.cw_file)[, c("patch_id","core_id")] else
  data.frame(patch_id = ps$patch_id, core_id = paste0("id", ps$patch_id))
ps$core_id <- core_cw$core_id[match(ps$patch_id, core_cw$patch_id)]
ps$core_id <- ifelse(is.na(ps$core_id), paste0("id", ps$patch_id), ps$core_id)

comp_f <- file.path(TB_OUT_TABLES, "27_conefor_components.csv")
if (file.exists(comp_f)) {
  comp <- read.csv(comp_f) |> dplyr::filter(d_km == FOCAL_D) |>
    dplyr::select(patch_id, dPCconnector, dIICconnector, dPCflux, dIICflux)
  ps <- dplyr::left_join(ps, comp, by = "patch_id")
}

bin <- terra::rast(file.path(TB_OUT_HS_BINARY, "present_wmean.tif"))
pat <- terra::patches(bin, directions = 8, zeroAsNA = TRUE); names(pat) <- "pid"
cell_km2 <- prod(terra::res(bin)) / 1e6

## ----------------------------------------------------------------------------
## 2) PA mask
## ----------------------------------------------------------------------------
tb_log_section("Protected-area mask")
pa_gpkg <- file.path(TB_OUT_VECTORS, "pa_combined.gpkg")
if (file.exists(pa_gpkg)) {
  pa_all <- sf::st_read(pa_gpkg, quiet = TRUE) |> sf::st_transform(TB_CRS_PROJ)
} else {
  tb_log("pa_combined.gpkg missing — rebuilding from PAs.gpkg", "WARN")
  avail <- tryCatch(sf::st_layers(TB_PA_GDB)$name, error = function(e) character())
  rd <- function(l) {
    if (!(l %in% avail)) return(NULL)
    v <- tryCatch(sf::st_read(TB_PA_GDB, layer = l, quiet = TRUE), error = function(e) NULL)
    if (is.null(v) || !nrow(v)) return(NULL)
    v <- sf::st_make_valid(v)
    v <- v[sf::st_geometry_type(v) %in%
             c("POLYGON","MULTIPOLYGON","GEOMETRYCOLLECTION","CURVEPOLYGON","MULTISURFACE"), ]
    if (!nrow(v)) return(NULL)
    v <- sf::st_transform(v, TB_CRS_PROJ); v$pa_layer <- l; v[, "pa_layer"]
  }
  pl <- lapply(PA_LAYERS, rd); pl <- pl[!sapply(pl, is.null)]
  pa_all <- sf::st_make_valid(do.call(rbind, pl))
}
pa_mask <- terra::rasterize(terra::vect(pa_all), bin, field = 1, background = 0)
names(pa_mask) <- "pa"

## ----------------------------------------------------------------------------
## 3) % of each core inside PA
## ----------------------------------------------------------------------------
tb_log_section("Core protection level")
## cells per patch, and protected cells per patch, via cross-tabulation
ct_all <- terra::freq(pat)                         # value=pid, count
ct_all <- as.data.frame(ct_all) |>
  dplyr::transmute(patch_id = value, cells_tot = count)
prot_cells <- terra::zonal(pa_mask, pat, fun = "sum", na.rm = TRUE)
names(prot_cells) <- c("patch_id", "cells_pa")
prot <- ct_all |>
  dplyr::left_join(prot_cells, by = "patch_id") |>
  dplyr::mutate(cells_pa = ifelse(is.na(cells_pa), 0, cells_pa),
                pct_in_pa = 100 * cells_pa / cells_tot)

ps <- ps |> dplyr::left_join(prot[, c("patch_id","pct_in_pa")], by = "patch_id")
ps$status <- dplyr::case_when(
  ps$pct_in_pa / 100 >= PROT_HI ~ "Protected",
  ps$pct_in_pa / 100 <  PROT_LO ~ "Unprotected",
  TRUE                          ~ "Partial")
tb_save_table(ps, "30_core_protection")
tb_log(sprintf("Protected=%d | Partial=%d | Unprotected=%d cores",
               sum(ps$status == "Protected"), sum(ps$status == "Partial"),
               sum(ps$status == "Unprotected")))

## ----------------------------------------------------------------------------
## 4) Graph at d=100 + PA-subnetwork connectivity
## ----------------------------------------------------------------------------
tb_log_section("PA-subnetwork connectivity")
coords <- as.matrix(ps[, c("x","y")]); area <- ps$area_km2
.mo  <- match(ps$patch_id, .cm[["present"]]$tab$patch_id)
D_km <- .cm[["present"]]$D[.mo, .mo]            # cost-km (resistance-weighted)
D_EFF <- FOCAL_D * R_EFF                         # effective dispersal threshold
Adj  <- (D_km <= D_EFF) & (D_km > 0)
g_full <- igraph::graph_from_adjacency_matrix(Adj, mode = "undirected")

tr_shp <- file.path(TB_DATA_ROOT, "TR_mask", "TR_mask.shp")
tr_sf  <- sf::st_read(tr_shp, quiet = TRUE) |> sf::st_transform(TB_CRS_PROJ)
AL_km2 <- as.numeric(sum(sf::st_area(tr_sf))) / 1e6

PC_of <- function(idx) {
  if (length(idx) < 2) return(if (length(idx) == 1) area[idx]^2 / AL_km2^2 else 0)
  Dk <- D_km[idx, idx, drop = FALSE]
  alpha <- -log(PROB_AT_THRESH) / D_EFF
  W <- -log(pmax(exp(-alpha * Dk), .Machine$double.eps)); diag(W) <- 0
  gg <- igraph::graph_from_adjacency_matrix(W, mode = "undirected", weighted = TRUE)
  Ps <- exp(-igraph::distances(gg, weights = igraph::E(gg)$weight))
  sum(outer(area[idx], area[idx]) * Ps) / AL_km2^2
}

prot_idx <- which(ps$status == "Protected")
g_pa <- igraph::induced_subgraph(g_full, prot_idx)
comp_pa <- igraph::components(g_pa)
pc_full <- PC_of(seq_len(n))
pc_pa   <- PC_of(prot_idx)

summary_df <- data.frame(
  n_cores            = n,
  n_protected        = length(prot_idx),
  n_partial          = sum(ps$status == "Partial"),
  n_unprotected      = sum(ps$status == "Unprotected"),
  pa_area_km2        = sum(area[prot_idx]),
  pa_components_d100  = comp_pa$no,
  pa_largest_comp_n   = max(comp_pa$csize),
  full_components_d100 = igraph::components(g_full)$no,
  PC_full            = pc_full,
  PC_pa_only         = pc_pa,
  PC_pa_share_pct    = 100 * pc_pa / pc_full)
tb_save_table(summary_df, "30_pa_network_summary")
tb_log(sprintf("PA-only network: %d components, captures %.1f%% of full PC",
               comp_pa$no, 100 * pc_pa / pc_full))

## ----------------------------------------------------------------------------
## 5) Unprotected priority cores
## ----------------------------------------------------------------------------
unprot <- ps |>
  dplyr::filter(status %in% c("Unprotected", "Partial")) |>
  dplyr::arrange(dplyr::desc(dPC)) |>
  dplyr::mutate(priority_rank = dplyr::row_number())
tb_save_table(unprot, "30_unprotected_priority")

## ----------------------------------------------------------------------------
## 6) FIGURES
## ----------------------------------------------------------------------------
tb_log_section("Figures")
PAL_STATUS <- c("Protected" = "#009E73", "Partial" = "#E69F00", "Unprotected" = "#9E2A2B")

## edges as line segments
edge_idx <- which(Adj & upper.tri(Adj), arr.ind = TRUE)
edges_df <- data.frame(
  x = coords[edge_idx[,1],1], y = coords[edge_idx[,1],2],
  xend = coords[edge_idx[,2],1], yend = coords[edge_idx[,2],2],
  s1 = ps$status[edge_idx[,1]], s2 = ps$status[edge_idx[,2]])
edges_df$etype <- dplyr::case_when(
  edges_df$s1 == "Protected" & edges_df$s2 == "Protected" ~ "PA–PA",
  edges_df$s1 == "Protected" | edges_df$s2 == "Protected" ~ "PA–gap",
  TRUE ~ "gap–gap")
PAL_EDGE <- c("PA–PA" = "#009E73", "PA–gap" = "#E69F00", "gap–gap" = "#BDBDBD")

nodes_sf <- sf::st_as_sf(ps, coords = c("x","y"), crs = TB_CRS_PROJ)
world_sf <- tryCatch(rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") |>
                       sf::st_transform(TB_CRS_PROJ), error = function(e) NULL)
e <- terra::ext(bin); pad <- 30000
xl <- c(e$xmin - pad, e$xmax + pad); yl <- c(e$ymin - pad, e$ymax + pad)

p30a <- ggplot()
if (!is.null(world_sf)) p30a <- p30a +
  geom_sf(data = world_sf, fill = "#F2F2F2", color = "#CFCFCF", linewidth = 0.3)
p30a <- p30a +
  geom_sf(data = pa_all, fill = "#009E73", color = NA, alpha = 0.18) +
  geom_sf(data = tr_sf, fill = NA, color = TB_COLOR_FRAME, linewidth = 0.5) +
  geom_segment(data = edges_df,
               aes(x = x, y = y, xend = xend, yend = yend, color = etype),
               linewidth = 0.35, alpha = 0.6) +
  scale_color_manual(values = PAL_EDGE, name = sprintf("Link (d=%d km)", FOCAL_D)) +
  geom_sf(data = nodes_sf, aes(size = dPC, fill = status),
          shape = 21, color = "black", stroke = 0.4, alpha = 0.9) +
  geom_sf_text(data = nodes_sf[order(-nodes_sf$dPC), ][1:6, ],
               aes(label = core_id), size = 2.6, color = "black",
               fontface = "bold", nudge_y = 28000) +
  scale_fill_manual(values = PAL_STATUS, name = "Core protection") +
  scale_size_continuous(name = "dPC (%)", range = c(2, 11)) +
  coord_sf(xlim = xl, ylim = yl, datum = sf::st_crs(4326), expand = FALSE) +
  tb_map_decorations() +
  labs(title = "Is the protected-area network connected?",
       subtitle = sprintf("Present cores by protection status (fill) and dPC importance (size); links at d = %d km. PA system = %d disconnected components capturing %.0f%% of national PC.",
                          FOCAL_D, comp_pa$no, 100 * pc_pa / pc_full)) +
  theme_trbear()
tb_save_fig(p30a, "fig30a_pa_gap_network", w = 15, h = 9, subdir = FIG_SUBDIR)

## ---- fig30b: unprotected priority bars ------------------------------------
## Top unprotected cores by dPC PLUS the genuine connector stepping stones
## (small cores with dIIC >> dPC) that a pure dPC ranking buries below the cut —
## these ARE priority new-PA candidates per the stepping-stone (dIIC) argument.
topU_dpc  <- unprot |> dplyr::slice_head(n = 10)
topU_step <- unprot |>
  dplyr::filter(area_km2 < 500, dIIC > dPC, !patch_id %in% topU_dpc$patch_id) |>
  dplyr::arrange(dplyr::desc(dIIC)) |>
  dplyr::slice_head(n = 3)
topU <- dplyr::bind_rows(topU_dpc, topU_step)
bar_long <- topU |>
  dplyr::transmute(lbl = factor(sprintf("%s (%.0f km², %.0f%% PA)",
                                        core_id, area_km2, pct_in_pa),
                                levels = rev(sprintf("%s (%.0f km², %.0f%% PA)",
                                                     core_id, area_km2, pct_in_pa))),
                   dPC, dIIC) |>
  tidyr::pivot_longer(c(dPC, dIIC), names_to = "index", values_to = "value")
p30b <- ggplot(bar_long, aes(value, lbl, fill = index)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  scale_fill_manual(values = c("dPC" = "#0072B2", "dIIC" = "#9E2A2B"), name = "Index") +
  labs(title = "Highest-priority UNPROTECTED cores (new-PA candidates)",
       subtitle = "Cores with < 50% area inside any PA: top 10 by full-network dPC, plus the three small connector stepping stones (dIIC ≫ dPC) that a dPC ranking buries below the cut.",
       x = "Patch importance (%)", y = NULL) +
  theme_trbear_bar(base_size = 12)
tb_save_fig(p30b, "fig30b_unprotected_priority", w = 13, h = 8, subdir = FIG_SUBDIR)

tb_log_session()
tb_log("30_pa_network_gap DONE")
