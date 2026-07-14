## ============================================================================
## 35_costdist_conefor.R
## Project: TR_Bear_Connectivity
## Purpose: Recompute the graph-theoretic connectivity indices (PC, IIC) and
##   per-patch importance (dPC, dIIC) using LEAST-COST (resistance-weighted,
##   "effective") distances between source patches INSTEAD of straight-line
##   Euclidean distances.
##
##   This corrects a limitation whereby the original graph
##   indices (script 17) were Euclidean and therefore only partly reflected the
##   resistance surface that defines the corridors. Here the SAME resistance
##   surface that UNICOR uses for the corridors is used to weight the graph, so
##   corridors and graph indices share one cost model.
##
## Engine: terra::costDist (least-cost accumulated cost over friction = R).
##   For each source patch i: mark its cells as target (cost 0), run costDist
##   over the scenario resistance surface, then take the MINIMUM accumulated
##   cost to every other patch j with terra::zonal -> patch-to-patch effective
##   distance matrix D_cost (cost-km).
##
## Threshold calibration (so cost-distance indices stay on a scale comparable
##   to the Euclidean version and to biological dispersal):
##     R_eff = median( D_cost_km / D_euclid_km ) over all present patch pairs
##     a geographic dispersal threshold d (km) maps to a cost threshold d*R_eff
##     alpha = -ln(0.5) / (d * R_eff);  p_ij = exp(-alpha * D_cost_km)
##   R_eff is fixed at its PRESENT value for all scenarios, so future increases
##   in matrix resistance correctly depress future connectivity.
##
## Validation: present terra::costDist matrix is correlated against the UNICOR
##   .cdmatrix.csv (independent least-cost engine) and r is logged + saved.
##
## Outputs (tables/):
##   35_costdist_indices.csv      scenario x d_km -> PC, IIC, n_patches, hab_km2
##   35_costdist_dpc_dii.csv      scenario x d_km x patch_id -> dPC, dIIC, area
##   35_costdist_calibration.csv  R_eff, validation correlation vs UNICOR
## Figures (figures/35_costdist/):
##   fig35a_pc_vs_distance.png    PC vs distance, lines per scenario (cost-dist)
##   fig35b_iic_vs_distance.png   IIC vs distance, lines per scenario
##   fig35c_euclid_vs_cost.png    present dPC: Euclidean (script 17) vs cost-dist
## ============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("igraph", quietly = TRUE))
    install.packages("igraph", repos = "https://cloud.r-project.org")
  library(terra); library(sf); library(ggplot2); library(dplyr); library(tidyr)
  library(igraph); library(ggrepel)
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
tb_log_init("35_costdist_conefor")

FIG_SUBDIR     <- "35_costdist"
DIST_KM        <- c(50, 100, 150, 200, 300, 400)
FOCAL_D        <- 100
PROB_AT_THRESH <- 0.5

scenarios <- c("present",
               sprintf("%s_%s", rep(TB_PERIODS, each = length(TB_SSPS)),
                                rep(TB_SSPS,    times = length(TB_PERIODS))))

.bin_path <- function(s) {
  if (s == "present") file.path(TB_OUT_HS_BINARY, "present_wmean.tif")
  else                file.path(TB_OUT_HS_BINARY, sprintf("future_%s.tif", s))
}
.res_path <- function(s) file.path(TB_OUT_RESISTANCE, sprintf("%s.tif", s))

## ---- landscape area --------------------------------------------------------
tr_mask_sf <- sf::st_read(file.path(TB_DATA_ROOT, "TR_mask", "TR_mask.shp"),
                          quiet = TRUE) |> sf::st_transform(TB_CRS_PROJ)
AL_km2 <- as.numeric(sum(sf::st_area(tr_mask_sf))) / 1e6
tb_log(sprintf("A_L = %.0f km2", AL_km2))

## ----------------------------------------------------------------------------
## Patch delineation (returns patches raster restricted to kept ids + table)
## ----------------------------------------------------------------------------
delineate <- function(s, min_km2 = TB_PATCH_MIN_KM2) {
  pp <- .bin_path(s); if (!file.exists(pp)) return(NULL)
  r <- terra::rast(pp); names(r) <- "suit"
  cell_km2 <- prod(terra::res(r)) / 1e6
  pat <- terra::patches(r, directions = 8, zeroAsNA = TRUE)
  fr  <- terra::freq(pat) |> as.data.frame() |>
    dplyr::rename(patch_id = value, n_cells = count) |>
    dplyr::mutate(area_km2 = n_cells * cell_km2) |>
    dplyr::filter(!is.na(patch_id), area_km2 >= min_km2)
  if (!nrow(fr)) return(NULL)
  keep <- terra::ifel(pat %in% fr$patch_id, pat, NA); names(keep) <- "patch"
  ## centroids
  cen <- lapply(fr$patch_id, function(i) {
    crd <- terra::crds(terra::as.points(terra::ifel(keep == i, 1, NA)))
    data.frame(patch_id = i, x = mean(crd[, 1]), y = mean(crd[, 2]))
  }) |> do.call(what = rbind)
  tab <- dplyr::left_join(fr[, c("patch_id", "area_km2")], cen, by = "patch_id")
  list(patches = keep, tab = tab, cell_km2 = cell_km2)
}

## ----------------------------------------------------------------------------
## Least-cost patch-to-patch distance matrix via terra::costDist + zonal(min)
##   D[i, j] = min accumulated cost (cost-m) from any cell of patch i to patch j
## ----------------------------------------------------------------------------
## CENTROID-to-centroid least-cost distance (matches UNICOR's source-point
## convention and the original script-17 centroid graph; isolates the
## resistance-weighting). Origin = the single centroid cell set to cost 0;
## destination cost read at every centroid with terra::extract.
cost_matrix <- function(res_r, cen_xy) {
  n <- nrow(cen_xy)
  cells <- terra::cellFromXY(res_r, cen_xy)
  D <- matrix(NA_real_, n, n)
  for (k in seq_len(n)) {
    orig <- res_r
    orig[cells[k]] <- 0                                  # origin centroid = cost 0
    cd   <- terra::costDist(orig, target = 0, scale = 1, maxiter = 3000)
    D[k, ] <- terra::extract(cd, cen_xy)[, 1]
    if (k %% 20 == 0) tb_log(sprintf("   costDist %d/%d", k, n))
  }
  D <- (D + t(D)) / 2            # enforce symmetry (numerical)
  ## unreachable pairs (separated by NA barriers, e.g. sea) -> large finite cost
  fin <- D[is.finite(D)]
  big <- if (length(fin)) max(fin) * 10 else 1e9
  D[!is.finite(D)] <- big
  diag(D) <- 0
  D / 1000                       # cost-km
}

## ---- PC / IIC on a precomputed cost-distance matrix ------------------------
compute_indices <- function(Dcost_km, area, d_km, R_eff, AL_km2) {
  n <- length(area)
  d_eff <- d_km * R_eff
  alpha <- -log(PROB_AT_THRESH) / d_eff
  ## PC : maximum-product path probability via Dijkstra on -log p
  W <- -log(pmax(exp(-alpha * Dcost_km), .Machine$double.eps)); diag(W) <- 0
  gP <- igraph::graph_from_adjacency_matrix(W, mode = "undirected", weighted = TRUE)
  Pstar <- exp(-igraph::distances(gP, weights = igraph::E(gP)$weight))
  PC <- sum(outer(area, area) * Pstar) / AL_km2^2
  ## IIC : binary graph, edge if Dcost <= d_eff
  Adj <- (Dcost_km <= d_eff) & (Dcost_km > 0)
  gI <- igraph::graph_from_adjacency_matrix(Adj, mode = "undirected")
  nl <- igraph::distances(gI)
  IIC <- sum(outer(area, area) / (1 + nl), na.rm = TRUE) / AL_km2^2
  ## node-removal importance
  dPC <- dIIC <- numeric(n)
  for (k in seq_len(n)) {
    gPk <- igraph::delete_vertices(gP, k)
    PCk <- sum(outer(area[-k], area[-k]) * exp(-igraph::distances(gPk,
              weights = igraph::E(gPk)$weight))) / AL_km2^2
    dPC[k] <- if (PC > 0) 100 * (PC - PCk) / PC else 0
    gIk <- igraph::delete_vertices(gI, k)
    IICk <- sum(outer(area[-k], area[-k]) / (1 + igraph::distances(gIk)),
                na.rm = TRUE) / AL_km2^2
    dIIC[k] <- if (IIC > 0) 100 * (IIC - IICk) / IIC else 0
  }
  list(PC = PC, IIC = IIC, dPC = dPC, dIIC = dIIC)
}

## ----------------------------------------------------------------------------
## 1) PRESENT: build cost matrix, validate vs UNICOR, calibrate R_eff
## ----------------------------------------------------------------------------
tb_log_section("Present cost matrix + calibration")
pres <- delineate("present")
res_pres <- terra::rast(.res_path("present"))
ids_pres <- pres$tab$patch_id
cen_pres <- as.matrix(pres$tab[, c("x", "y")])
Dcost_pres <- cost_matrix(res_pres, cen_pres)
## centroid Euclidean (km) — same geometry as the cost distances, used to
## calibrate the dispersal thresholds and for the Euclidean-vs-cost comparison
Deuc_pres <- as.matrix(stats::dist(cen_pres)) / 1000

## validation against UNICOR cdmatrix (rows follow 13_source_points order)
sp_order <- read.csv(file.path(TB_OUT_TABLES, "13_source_points.csv"))
cdm_file <- file.path(TB_OUT_UNICOR_DIR, "present", "results",
                      "present_sources.cdmatrix.csv")
val_r <- NA_real_
if (file.exists(cdm_file)) {
  cdm <- as.matrix(read.csv(cdm_file, header = FALSE))   # cost-m, in sp_order
  ## bring OUR matrix into the UNICOR (sp_order) ordering, then correlate
  ord <- match(sp_order$patch_id, ids_pres)              # sp row -> Dcost index
  if (nrow(cdm) == length(ord) && !any(is.na(ord))) {
    Dc <- Dcost_pres[ord, ord]
    ut <- upper.tri(Dc)
    a <- Dc[ut]; b <- cdm[ut]
    ok <- is.finite(a) & is.finite(b)
    val_r <- if (sum(ok) > 3) suppressWarnings(cor(a[ok], b[ok])) else NA_real_
    tb_log(sprintf("validation cor(terra costDist, UNICOR cdmatrix) = %s (%d pairs)",
                   ifelse(is.na(val_r), "NA", sprintf("%.4f", val_r)), sum(ok)))
  } else tb_log("UNICOR cdmatrix ordering mismatch; skipping validation", "WARN")
} else tb_log("UNICOR present cdmatrix not found; skipping validation", "WARN")

ut <- upper.tri(Dcost_pres)
ratio <- Dcost_pres[ut] / Deuc_pres[ut]
ratio <- ratio[is.finite(ratio) & Deuc_pres[ut] > 0]
R_eff <- stats::median(ratio, na.rm = TRUE)
tb_log(sprintf("R_eff (median cost/Euclidean, centroid) = %.3f", R_eff))

tb_save_table(data.frame(R_eff = R_eff, validation_cor = val_r,
                         n_pairs = sum(ut), n_patches = length(ids_pres)),
              "35_costdist_calibration")

## ----------------------------------------------------------------------------
## 2) ALL scenarios: cost matrices + indices
## ----------------------------------------------------------------------------
tb_log_section("All scenarios")
idx_rows <- list(); imp_rows <- list()
cost_cache <- list(present = list(D = Dcost_pres, tab = pres$tab))

for (s in scenarios) {
  if (s == "present") { dl <- pres; D <- Dcost_pres }
  else {
    dl <- delineate(s); if (is.null(dl)) { tb_log(sprintf("[%s] no patches", s), "WARN"); next }
    res_s <- terra::rast(.res_path(s))
    tb_log(sprintf("[%s] cost matrix (%d patches)", s, nrow(dl$tab)))
    D <- cost_matrix(res_s, as.matrix(dl$tab[, c("x", "y")]))
    cost_cache[[s]] <- list(D = D, tab = dl$tab)
  }
  area <- dl$tab$area_km2
  for (d in DIST_KM) {
    r <- compute_indices(D, area, d, R_eff, AL_km2)
    idx_rows[[paste(s, d)]] <- data.frame(scenario = s, d_km = d,
      n_patches = length(area), habitat_km2 = sum(area),
      PC = r$PC, IIC = r$IIC)
    imp_rows[[paste(s, d)]] <- data.frame(scenario = s, d_km = d,
      patch_id = dl$tab$patch_id, x = dl$tab$x, y = dl$tab$y,
      area_km2 = area, dPC = r$dPC, dIIC = r$dIIC)
  }
}
idx_df <- do.call(rbind, idx_rows); imp_df <- do.call(rbind, imp_rows)
tb_save_table(idx_df, "35_costdist_indices")
tb_save_table(imp_df, "35_costdist_dpc_dii")
saveRDS(cost_cache, file.path(TB_OUT_RDS, "35_cost_matrices.rds"))

## ----------------------------------------------------------------------------
## 3) FIGURES
## ----------------------------------------------------------------------------
tb_log_section("Figures")
scen_label <- c("present" = "Present",
  "2041_2070_ssp126" = "2070s SSP126", "2041_2070_ssp370" = "2070s SSP370",
  "2041_2070_ssp585" = "2070s SSP585", "2071_2100_ssp126" = "2100s SSP126",
  "2071_2100_ssp370" = "2100s SSP370", "2071_2100_ssp585" = "2100s SSP585")
PAL <- c("Present"="#000000","2070s SSP126"="#56B4E9","2070s SSP370"="#E69F00",
  "2070s SSP585"="#D55E00","2100s SSP126"="#0072B2","2100s SSP370"="#CC79A7",
  "2100s SSP585"="#9E2A2B")
idx_df$lbl <- factor(scen_label[idx_df$scenario], levels = scen_label[scenarios])

mk <- function(yv, ylab, ttl) ggplot(idx_df, aes(d_km, .data[[yv]], color = lbl)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.3) +
  scale_color_manual(values = PAL, name = "Scenario") +
  scale_x_continuous(breaks = DIST_KM) +
  labs(title = ttl, x = "Dispersal threshold (effective km, p = 0.5 at threshold)",
       y = ylab) + theme_trbear_bar(base_size = 12)
tb_save_fig(mk("PC", "PC (cost-distance)",
    "Probability of Connectivity — least-cost (effective) distance"),
    "fig35a_pc_vs_distance", w = 12, h = 7, subdir = FIG_SUBDIR)
tb_save_fig(mk("IIC", "IIC (cost-distance)",
    "Integral Index of Connectivity — least-cost (effective) distance"),
    "fig35b_iic_vs_distance", w = 12, h = 7, subdir = FIG_SUBDIR)

## present dPC: Euclidean (script 17) vs cost-distance
euc <- tryCatch(read.csv(file.path(TB_OUT_TABLES, "17_conefor_dpc_dii.csv")) |>
  dplyr::filter(scenario == "present", d_km == FOCAL_D) |>
  dplyr::select(patch_id, dPC_euc = dPC, dIIC_euc = dIIC), error = function(e) NULL)
cw <- tryCatch(read.csv(file.path(TB_OUT_TABLES, "34_core_crosswalk.csv"))[,
  c("patch_id","core_id")], error = function(e) NULL)
if (!is.null(euc)) {
  cmp <- imp_df |> dplyr::filter(scenario == "present", d_km == FOCAL_D) |>
    dplyr::select(patch_id, dPC_cost = dPC, dIIC_cost = dIIC) |>
    dplyr::left_join(euc, by = "patch_id")
  if (!is.null(cw)) cmp <- dplyr::left_join(cmp, cw, by = "patch_id")
  lim <- range(c(cmp$dPC_euc, cmp$dPC_cost), na.rm = TRUE)
  p35c <- ggplot(cmp, aes(dPC_euc, dPC_cost)) +
    geom_abline(slope = 1, linetype = 2, color = "gray55") +
    geom_point(color = "#0072B2", size = 2.6, alpha = 0.8) +
    { if (!is.null(cw)) ggrepel::geom_text_repel(
        data = dplyr::filter(cmp, dPC_euc > 2 | dPC_cost > 2),
        aes(label = core_id), size = 3, max.overlaps = 20) } +
    coord_equal(xlim = lim, ylim = lim) +
    labs(title = "Per-core dPC: Euclidean vs least-cost distance (present, d = 100 km)",
         x = "dPC (Euclidean, original)", y = "dPC (least-cost, revised)") +
    theme_trbear_bar(base_size = 12)
  tb_save_fig(p35c, "fig35c_euclid_vs_cost", w = 9, h = 8, subdir = FIG_SUBDIR)
}

tb_log_session(); tb_log("35_costdist_conefor DONE")
