## ============================================================================
## 37_patch_threshold.R
## Project: TR_Bear_Connectivity
## Purpose: Sensitivity of the present cost-distance connectivity network to the
##   SOURCE-PATCH minimum-size threshold. The manuscript uses 83 km2 (upper
##   bound of NE-Anatolian bear home ranges). To gauge robustness to this
##   choice, we repeat the delineation at 50, 83 and
##   120 km2 and recompute PC, IIC, the two-core dominance, and n50 over the
##   least-cost (effective) distance graph (present resistance, c = 4).
##
## Outputs:
##   tables/37_patch_threshold.csv   min_km2 x d_km -> n_patches, hab_km2,
##                                                     PC, IIC, n50_PC,
##                                                     dPC_top1, dPC_top2
##   figures/37_patch_threshold/fig37_patch_threshold.png
## ============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("igraph", quietly = TRUE))
    install.packages("igraph", repos = "https://cloud.r-project.org")
  library(terra); library(sf); library(ggplot2); library(dplyr); library(tidyr); library(igraph)
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
tb_log_init("37_patch_threshold")

FIG_SUBDIR     <- "37_patch_threshold"
MIN_KM2        <- c(50, 83, 120)
DIST_KM        <- c(50, 100, 200)
FOCAL_D        <- 100
PROB_AT_THRESH <- 0.5

tr_mask_sf <- sf::st_read(file.path(TB_DATA_ROOT, "TR_mask", "TR_mask.shp"),
                          quiet = TRUE) |> sf::st_transform(TB_CRS_PROJ)
AL_km2 <- as.numeric(sum(sf::st_area(tr_mask_sf))) / 1e6

r <- terra::rast(file.path(TB_OUT_HS_BINARY, "present_wmean.tif")); names(r) <- "suit"
cell_km2 <- prod(terra::res(r)) / 1e6
pat <- terra::patches(r, directions = 8, zeroAsNA = TRUE)
freq_all <- terra::freq(pat) |> as.data.frame() |>
  dplyr::rename(patch_id = value, n_cells = count) |>
  dplyr::mutate(area_km2 = n_cells * cell_km2) |>
  dplyr::filter(!is.na(patch_id))

res_pres <- terra::rast(file.path(TB_OUT_RESISTANCE, "present.tif"))
cal <- tryCatch(read.csv(file.path(TB_OUT_TABLES, "35_costdist_calibration.csv")),
                error = function(e) NULL)

cost_matrix <- function(cen_xy) {                 # centroid-to-centroid (matches script 35)
  n <- nrow(cen_xy); cells <- terra::cellFromXY(res_pres, cen_xy)
  D <- matrix(NA_real_, n, n)
  for (k in seq_len(n)) {
    orig <- res_pres; orig[cells[k]] <- 0
    cd <- terra::costDist(orig, target = 0, scale = 1, maxiter = 3000)
    D[k, ] <- terra::extract(cd, cen_xy)[, 1]
  }
  D <- (D + t(D)) / 2
  fin <- D[is.finite(D)]; big <- if (length(fin)) max(fin) * 10 else 1e9
  D[!is.finite(D)] <- big; diag(D) <- 0; D / 1000
}
patch_centroids <- function(keep, ids) {
  do.call(rbind, lapply(ids, function(i) {
    crd <- terra::crds(terra::as.points(terra::ifel(keep == i, 1, NA)))
    c(mean(crd[, 1]), mean(crd[, 2]))
  }))
}
PC_of <- function(D, area, keepv, d_eff) {
  ar <- area[keepv]; if (length(ar) < 2) return(if (length(ar)) ar^2/AL_km2^2 else 0)
  alpha <- -log(PROB_AT_THRESH) / d_eff
  W <- -log(pmax(exp(-alpha * D[keepv, keepv]), .Machine$double.eps)); diag(W) <- 0
  g <- igraph::graph_from_adjacency_matrix(W, mode = "undirected", weighted = TRUE)
  sum(outer(ar, ar) * exp(-igraph::distances(g, weights = igraph::E(g)$weight))) / AL_km2^2
}
IIC_of <- function(D, area, d_eff) {
  Adj <- (D <= d_eff) & (D > 0)
  g <- igraph::graph_from_adjacency_matrix(Adj, mode = "undirected")
  sum(outer(area, area) / (1 + igraph::distances(g)), na.rm = TRUE) / AL_km2^2
}
dPC_vec <- function(D, area, d_eff) {
  n <- length(area); pc0 <- PC_of(D, area, seq_len(n), d_eff)
  vapply(seq_len(n), function(k)
    if (pc0 > 0) 100*(pc0 - PC_of(D, area, setdiff(seq_len(n), k), d_eff))/pc0 else 0, numeric(1))
}
n50_targeted <- function(D, area, d_eff) {
  n <- length(area); pc0 <- PC_of(D, area, seq_len(n), d_eff)
  ord <- order(dPC_vec(D, area, d_eff), decreasing = TRUE); remaining <- seq_len(n)
  for (i in seq_len(n - 1L)) {
    remaining <- setdiff(remaining, ord[i])
    if (PC_of(D, area, remaining, d_eff) < 0.5 * pc0) return(i)
  }
  NA_integer_
}

rows <- list()
for (mk in MIN_KM2) {
  tb_log_section(sprintf("min = %g km2", mk))
  fr <- freq_all |> dplyr::filter(area_km2 >= mk)
  keep <- terra::ifel(pat %in% fr$patch_id, pat, NA); names(keep) <- "patch"
  ids <- fr$patch_id; area <- fr$area_km2
  tb_log(sprintf("  n_patches = %d (total %.0f km2)", length(ids), sum(area)))
  cen_xy <- patch_centroids(keep, ids)
  D <- cost_matrix(cen_xy)
  R_eff <- if (!is.null(cal)) cal$R_eff[1] else 1
  for (d in DIST_KM) {
    d_eff <- d * R_eff
    dpc <- dPC_vec(D, area, d_eff); top2 <- sort(dpc, decreasing = TRUE)[1:2]
    rows[[paste(mk, d)]] <- data.frame(min_km2 = mk, d_km = d,
      n_patches = length(ids), habitat_km2 = sum(area),
      PC = PC_of(D, area, seq_len(length(ids)), d_eff),
      IIC = IIC_of(D, area, d_eff),
      n50_PC = n50_targeted(D, area, d_eff),
      dPC_top1 = top2[1], dPC_top2 = top2[2])
    tb_log(sprintf("  d=%d done", d))
  }
}
out <- do.call(rbind, rows)
tb_save_table(out, "37_patch_threshold")

## ---- figure ----------------------------------------------------------------
tb_log_section("Figure")
long <- out |> tidyr::pivot_longer(c(PC, IIC, n50_PC), names_to = "metric", values_to = "val")
long$metric <- factor(long$metric, levels = c("PC","IIC","n50_PC"),
                      labels = c("PC","IIC","n50 (cores to halve PC)"))
p <- ggplot(long, aes(factor(min_km2), val, fill = factor(d_km))) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  facet_wrap(~ metric, scales = "free_y") +
  scale_fill_manual(values = c("50"="#56B4E9","100"="#000000","200"="#D55E00"),
                    name = "d (km)") +
  labs(title = "Sensitivity to the source-patch size threshold (present, least-cost graph)",
       x = expression("Minimum source-patch size ("*km^2*")"), y = NULL) +
  theme_trbear_bar(base_size = 12)
tb_save_fig(p, "fig37_patch_threshold", w = 13, h = 6, subdir = FIG_SUBDIR)

tb_log_session(); tb_log("37_patch_threshold DONE")
