## ============================================================================
## 36_c_sensitivity.R
## Project: TR_Bear_Connectivity
## Purpose: Sensitivity of the cost-distance connectivity results to the
##   resistance shape constant c transfer function
##       R(h) = 100 - 99 * (1 - exp(-c*h)) / (1 - exp(-c))
##   The manuscript uses c = 4; we test how PC, IIC and the network
##   robustness (n50) respond to c. Here we re-derive the PRESENT resistance for
##   c in {2, 4, 6, 8}, rebuild the least-cost distance matrix (terra::costDist)
##   and recompute PC, IIC and n50 (cores whose targeted removal halves PC).
##
##   NB: source patches come from the binary suitability map (MAX_TSS) and are
##   therefore IDENTICAL across c; only the matrix resistance (and hence the
##   effective distances among patches) changes with c.
##
## Outputs:
##   tables/36_c_sensitivity.csv     c x d_km -> PC, IIC, n50_PC, R_eff_c
##   figures/36_c_sensitivity/fig36_c_sensitivity.png
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
tb_log_init("36_c_sensitivity")

FIG_SUBDIR     <- "36_c_sensitivity"
C_VALUES       <- c(2, 4, 6, 8)
DIST_KM        <- c(50, 100, 200)
FOCAL_D        <- 100
PROB_AT_THRESH <- 0.5

.hs_to_R <- function(h, c, rmin = TB_RESIST_MIN, rmax = TB_RESIST_MAX)
  rmax - (rmax - rmin) * (1 - exp(-c * h)) / (1 - exp(-c))

tr_mask_sf <- sf::st_read(file.path(TB_DATA_ROOT, "TR_mask", "TR_mask.shp"),
                          quiet = TRUE) |> sf::st_transform(TB_CRS_PROJ)
AL_km2 <- as.numeric(sum(sf::st_area(tr_mask_sf))) / 1e6

## present patches (binary, independent of c)
r <- terra::rast(file.path(TB_OUT_HS_BINARY, "present_wmean.tif")); names(r) <- "suit"
cell_km2 <- prod(terra::res(r)) / 1e6
pat <- terra::patches(r, directions = 8, zeroAsNA = TRUE)
fr  <- terra::freq(pat) |> as.data.frame() |>
  dplyr::rename(patch_id = value, n_cells = count) |>
  dplyr::mutate(area_km2 = n_cells * cell_km2) |>
  dplyr::filter(!is.na(patch_id), area_km2 >= TB_PATCH_MIN_KM2)
keep <- terra::ifel(pat %in% fr$patch_id, pat, NA); names(keep) <- "patch"
ids  <- fr$patch_id; area <- fr$area_km2
cen <- lapply(ids, function(i) {
  crd <- terra::crds(terra::as.points(terra::ifel(keep == i, 1, NA)))
  data.frame(patch_id = i, x = mean(crd[, 1]), y = mean(crd[, 2]))
}) |> do.call(what = rbind)
Deuc <- as.matrix(stats::dist(as.matrix(cen[, c("x", "y")]))) / 1000
tb_log(sprintf("present patches: %d", length(ids)))

## present habitat suitability (continuous) for resistance regeneration
hs <- terra::rast(file.path(TB_OUT_HS_PRESENT, "wmean.tif"))
if (terra::nlyr(hs) > 1) hs <- hs[[1]]

cen_xy <- as.matrix(cen[, c("x", "y")])     # centroid-to-centroid (matches script 35)
cost_matrix <- function(res_r) {
  n <- nrow(cen_xy); cells <- terra::cellFromXY(res_r, cen_xy)
  D <- matrix(NA_real_, n, n)
  for (k in seq_len(n)) {
    orig <- res_r; orig[cells[k]] <- 0
    cd <- terra::costDist(orig, target = 0, scale = 1, maxiter = 3000)
    D[k, ] <- terra::extract(cd, cen_xy)[, 1]
  }
  D <- (D + t(D)) / 2
  fin <- D[is.finite(D)]; big <- if (length(fin)) max(fin) * 10 else 1e9
  D[!is.finite(D)] <- big; diag(D) <- 0; D / 1000
}
PC_of <- function(D, keepv, d_eff) {
  ar <- area[keepv]; if (length(ar) < 2) return(if (length(ar)) ar^2/AL_km2^2 else 0)
  alpha <- -log(PROB_AT_THRESH) / d_eff
  W <- -log(pmax(exp(-alpha * D[keepv, keepv]), .Machine$double.eps)); diag(W) <- 0
  g <- igraph::graph_from_adjacency_matrix(W, mode = "undirected", weighted = TRUE)
  sum(outer(ar, ar) * exp(-igraph::distances(g, weights = igraph::E(g)$weight))) / AL_km2^2
}
IIC_of <- function(D, d_eff) {
  Adj <- (D <= d_eff) & (D > 0)
  g <- igraph::graph_from_adjacency_matrix(Adj, mode = "undirected")
  sum(outer(area, area) / (1 + igraph::distances(g)), na.rm = TRUE) / AL_km2^2
}
dPC_vec <- function(D, d_eff) {
  n <- length(ids); pc0 <- PC_of(D, seq_len(n), d_eff)
  vapply(seq_len(n), function(k)
    if (pc0 > 0) 100*(pc0 - PC_of(D, setdiff(seq_len(n), k), d_eff))/pc0 else 0, numeric(1))
}
n50_targeted <- function(D, d_eff) {
  n <- length(ids); pc0 <- PC_of(D, seq_len(n), d_eff)
  ord <- order(dPC_vec(D, d_eff), decreasing = TRUE)
  remaining <- seq_len(n)
  for (i in seq_len(n - 1L)) {
    remaining <- setdiff(remaining, ord[i])
    if (PC_of(D, remaining, d_eff) < 0.5 * pc0) return(i)
  }
  NA_integer_
}

## R_eff from script 35 (present, c=4); fallback recompute below per c
cal <- tryCatch(read.csv(file.path(TB_OUT_TABLES, "35_costdist_calibration.csv")),
                error = function(e) NULL)

rows <- list()
for (cc in C_VALUES) {
  tb_log_section(sprintf("c = %g", cc))
  res_c <- terra::app(hs, function(h) .hs_to_R(h, cc)); names(res_c) <- "R"
  D <- cost_matrix(res_c)
  ut <- upper.tri(D); ratio <- D[ut] / Deuc[ut]
  R_eff_c <- stats::median(ratio[is.finite(ratio) & Deuc[ut] > 0], na.rm = TRUE)
  ## use the manuscript calibration (present c=4 R_eff) so thresholds are fixed;
  ## fall back to per-c R_eff if script 35 output is unavailable
  R_eff <- if (!is.null(cal)) cal$R_eff[1] else R_eff_c
  for (d in DIST_KM) {
    d_eff <- d * R_eff
    rows[[paste(cc, d)]] <- data.frame(c = cc, d_km = d,
      PC = PC_of(D, seq_len(length(ids)), d_eff),
      IIC = IIC_of(D, d_eff),
      n50_PC = n50_targeted(D, d_eff),
      R_eff_c = R_eff_c)
    tb_log(sprintf("  d=%d PC done", d))
  }
}
out <- do.call(rbind, rows)
tb_save_table(out, "36_c_sensitivity")

## ---- figure ----------------------------------------------------------------
tb_log_section("Figure")
foc <- out |> dplyr::filter(d_km == FOCAL_D)
long <- out |> tidyr::pivot_longer(c(PC, IIC), names_to = "index", values_to = "val")
p1 <- ggplot(long, aes(c, val, color = factor(d_km))) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.3) +
  facet_wrap(~ index, scales = "free_y") +
  scale_color_manual(values = c("50"="#56B4E9","100"="#000000","200"="#D55E00"),
                     name = "d (km)") +
  scale_x_continuous(breaks = C_VALUES) +
  labs(title = "Connectivity sensitivity to resistance shape constant c",
       x = "Resistance shape constant c", y = "Index value") +
  theme_trbear_bar(base_size = 12)
p2 <- ggplot(out, aes(c, n50_PC, color = factor(d_km), group = factor(d_km))) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.6) +
  scale_color_manual(values = c("50"="#56B4E9","100"="#000000","200"="#D55E00"),
                     name = "d (km)") +
  scale_x_continuous(breaks = C_VALUES) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  labs(title = "Network robustness (n50) vs c", x = "Resistance shape constant c",
       y = expression(n[50]~"(cores to halve PC)")) +
  theme_trbear_bar(base_size = 12)
p <- patchwork::wrap_plots(p1, p2, ncol = 1, heights = c(1, 0.8)) +
  patchwork::plot_annotation(tag_levels = "a")
tb_save_fig(p, "fig36_c_sensitivity", w = 11, h = 10, subdir = FIG_SUBDIR)

tb_log_session(); tb_log("36_c_sensitivity DONE")
