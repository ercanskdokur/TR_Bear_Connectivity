## ============================================================================
## 29_network_robustness.R
## Project: TR_Bear_Connectivity
## Purpose: Quantify the STRUCTURAL FRAGILITY of the present connectivity
##   network via sequential node (patch) removal — a percolation / network
##   robustness experiment.
##
##   The present national network concentrates connectivity in a few cores
##   (Pareto). We test how fast global connectivity (PC, IIC) collapses when
##   patches are removed in three orders:
##     (1) TARGETED  — most-important-first (descending present dPC)
##     (2) RANDOM    — uniform random order (B permutations; envelope + median)
##     (3) BY AREA   — largest-first (to show importance != area)
##   A network where targeted removal crashes PC far faster than random removal
##   is fragile and attack-vulnerable.
##
##   Robustness metric R = area under the targeted PC-retention curve
##   (1 = perfectly robust, ~0 = collapses on first removal). We also report
##   n50 = number of cores whose removal drops PC below 50% of its initial value.
##
## Outputs (tables/):
##   29_robustness_curves.csv    strategy x d_km x n_removed -> frac_PC, frac_IIC
##   29_robustness_summary.csv   strategy x d_km -> R_PC, R_IIC, n50_PC
## Figures (figures/29_robustness/):
##   fig29a_robustness_curve.png  PC & IIC retention vs # cores removed (d=100)
##   fig29b_robustness_bydist.png PC retention faceted across dispersal distances
## ============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("igraph", quietly = TRUE))
    install.packages("igraph", repos = "https://cloud.r-project.org")
  library(igraph); library(sf); library(dplyr); library(tidyr)
  library(ggplot2); library(patchwork)
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
tb_log_init("29_network_robustness")

FIG_SUBDIR <- "29_robustness"
DIST_KM    <- c(50, 100, 200)
FOCAL_D    <- 100
B_RANDOM   <- 200L
PROB_AT_THRESH <- 0.5
set.seed(42)

## ----------------------------------------------------------------------------
## Present patches + landscape area
## ----------------------------------------------------------------------------
tb_log_section("Inputs")
dpc <- read.csv(file.path(TB_OUT_TABLES, "35_costdist_dpc_dii.csv"))  # cost-distance
ps  <- dpc |> dplyr::filter(scenario == "present", d_km == FOCAL_D) |>
  dplyr::select(patch_id, x, y, area_km2, dPC)
coords <- as.matrix(ps[, c("x", "y")])
area   <- ps$area_km2
n      <- nrow(ps)
## cost-distance matrix + calibration from script 35
.cm   <- readRDS(file.path(TB_OUT_RDS, "35_cost_matrices.rds"))
R_EFF <- read.csv(file.path(TB_OUT_TABLES, "35_costdist_calibration.csv"))$R_eff[1]

tr_shp <- file.path(TB_DATA_ROOT, "TR_mask", "TR_mask.shp")
AL_km2 <- as.numeric(sum(sf::st_area(
  sf::st_transform(sf::st_read(tr_shp, quiet = TRUE), TB_CRS_PROJ)))) / 1e6
tb_log(sprintf("n_patches = %d | A_L = %.0f km2", n, AL_km2))

.mo       <- match(ps$patch_id, .cm[["present"]]$tab$patch_id)
D_km_full <- .cm[["present"]]$D[.mo, .mo]   # cost-km (resistance-weighted)

## ----------------------------------------------------------------------------
## Index of a subset of patches
## ----------------------------------------------------------------------------
PC_subset <- function(keep, d_km) {
  m <- length(keep)
  if (m == 0) return(0)
  ar <- area[keep]
  if (m == 1) return(ar^2 / AL_km2^2)
  Dk <- D_km_full[keep, keep, drop = FALSE]
  alpha <- -log(PROB_AT_THRESH) / (d_km * R_EFF)
  W <- -log(pmax(exp(-alpha * Dk), .Machine$double.eps)); diag(W) <- 0
  g <- igraph::graph_from_adjacency_matrix(W, mode = "undirected", weighted = TRUE)
  Pstar <- exp(-igraph::distances(g, weights = igraph::E(g)$weight))
  sum(outer(ar, ar) * Pstar) / AL_km2^2
}
IIC_subset <- function(keep, d_km) {
  m <- length(keep)
  if (m == 0) return(0)
  ar <- area[keep]
  if (m == 1) return(ar^2 / AL_km2^2)
  Dk <- D_km_full[keep, keep, drop = FALSE]
  Adj <- (Dk <= d_km * R_EFF) & (Dk > 0)
  g <- igraph::graph_from_adjacency_matrix(Adj, mode = "undirected")
  nl <- igraph::distances(g)
  sum(outer(ar, ar) / (1 + nl), na.rm = TRUE) / AL_km2^2
}

## retention curve for a given removal order (vector of patch indices removed
## in sequence); returns data.frame n_removed (0..n-1) -> frac_PC, frac_IIC
retention_curve <- function(order_idx, d_km, pc0, iic0) {
  remaining <- seq_len(n)
  out <- data.frame(n_removed = 0L, frac_PC = 1, frac_IIC = 1)
  for (i in seq_len(n - 1L)) {
    remaining <- setdiff(remaining, order_idx[i])
    out <- rbind(out, data.frame(
      n_removed = i,
      frac_PC   = PC_subset(remaining, d_km) / pc0,
      frac_IIC  = IIC_subset(remaining, d_km) / iic0))
  }
  out
}

## ----------------------------------------------------------------------------
## Run all strategies x distances
## ----------------------------------------------------------------------------
tb_log_section("Removal experiments")

curve_rows <- list(); summ_rows <- list()
for (d in DIST_KM) {
  pc0  <- PC_subset(seq_len(n), d)
  iic0 <- IIC_subset(seq_len(n), d)
  tb_log(sprintf("[d=%d] PC0=%.4g IIC0=%.4g", d, pc0, iic0))

  ## dPC order is distance-specific
  dpc_d <- dpc |> dplyr::filter(scenario == "present", d_km == d)
  dpc_d <- dpc_d[match(ps$patch_id, dpc_d$patch_id), ]
  ord_targeted <- order(dpc_d$dPC, decreasing = TRUE)
  ord_area     <- order(area, decreasing = TRUE)

  c_t <- retention_curve(ord_targeted, d, pc0, iic0); c_t$strategy <- "Targeted (dPC)"
  c_a <- retention_curve(ord_area,     d, pc0, iic0); c_a$strategy <- "By area"

  ## random envelope
  rand_pc <- matrix(NA_real_, nrow = n, ncol = B_RANDOM)
  rand_iic <- matrix(NA_real_, nrow = n, ncol = B_RANDOM)
  for (b in seq_len(B_RANDOM)) {
    cc <- retention_curve(sample(seq_len(n)), d, pc0, iic0)
    rand_pc[, b]  <- cc$frac_PC
    rand_iic[, b] <- cc$frac_IIC
  }
  c_r <- data.frame(
    n_removed = 0:(n - 1L),
    frac_PC   = apply(rand_pc,  1, median),
    frac_IIC  = apply(rand_iic, 1, median),
    pc_lo     = apply(rand_pc,  1, quantile, 0.025),
    pc_hi     = apply(rand_pc,  1, quantile, 0.975),
    strategy  = "Random (median, 95% env)")
  c_t$pc_lo <- c_t$pc_hi <- NA; c_a$pc_lo <- c_a$pc_hi <- NA

  cd <- dplyr::bind_rows(c_t, c_a, c_r); cd$d_km <- d
  curve_rows[[as.character(d)]] <- cd

  ## summary metrics
  R_of <- function(v) mean(v)                       # area under retention curve
  n50  <- function(v) { w <- which(v < 0.5); if (length(w)) min(w) - 1L else NA_integer_ }
  for (st in unique(cd$strategy)) {
    sub <- cd[cd$strategy == st, ]
    summ_rows[[paste(d, st)]] <- data.frame(
      d_km = d, strategy = st,
      R_PC = R_of(sub$frac_PC), R_IIC = R_of(sub$frac_IIC),
      n50_PC = n50(sub$frac_PC))
  }
}
curves <- dplyr::bind_rows(curve_rows)
summ   <- dplyr::bind_rows(summ_rows)
tb_save_table(curves, "29_robustness_curves")
tb_save_table(summ,   "29_robustness_summary")

## ----------------------------------------------------------------------------
## FIGURES
## ----------------------------------------------------------------------------
tb_log_section("Figures")
PAL_ST <- c("Targeted (dPC)" = "#9E2A2B",
            "By area"        = "#E69F00",
            "Random (median, 95% env)" = "#0072B2")

cd <- curves |> dplyr::filter(d_km == FOCAL_D)
env <- cd |> dplyr::filter(strategy == "Random (median, 95% env)")

mk_panel <- function(yvar, ylab, ttl) {
  ggplot(cd, aes(n_removed, .data[[yvar]], color = strategy)) +
    { if (yvar == "frac_PC")
        geom_ribbon(data = env, aes(ymin = pc_lo, ymax = pc_hi),
                    fill = "#0072B2", alpha = 0.15, color = NA) } +
    geom_line(linewidth = 1.1) +
    geom_hline(yintercept = 0.5, linetype = 3, color = "gray40") +
    scale_color_manual(values = PAL_ST, name = "Removal order") +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(title = ttl, x = "Number of core patches removed", y = ylab) +
    theme_trbear_bar(base_size = 12)
}
p_pc  <- mk_panel("frac_PC",  "PC retained (% of initial)",  "A. Probability of Connectivity")
p_iic <- mk_panel("frac_IIC", "IIC retained (% of initial)", "B. Integral Index of Connectivity")
p29a <- (p_pc + p_iic) +
  patchwork::plot_annotation(
    title = sprintf("Network robustness to sequential core loss (present, d = %d km)", FOCAL_D),
    subtitle = "Targeted removal collapsing far faster than random = a fragile, attack-vulnerable network.",
    theme = theme(plot.title = element_text(face = "bold", size = 15, color = TB_COLOR_FRAME)))
tb_save_fig(p29a, "fig29a_robustness_curve", w = 16, h = 7, subdir = FIG_SUBDIR)

## ---- fig29b: PC retention across dispersal distances -----------------------
cd2 <- curves
cd2$d_lbl <- factor(sprintf("d = %d km", cd2$d_km),
                    levels = sprintf("d = %d km", DIST_KM))
env2 <- cd2 |> dplyr::filter(strategy == "Random (median, 95% env)")
p29b <- ggplot(cd2, aes(n_removed, frac_PC, color = strategy)) +
  geom_ribbon(data = env2, aes(ymin = pc_lo, ymax = pc_hi),
              fill = "#0072B2", alpha = 0.15, color = NA) +
  geom_line(linewidth = 1.0) +
  geom_hline(yintercept = 0.5, linetype = 3, color = "gray40") +
  scale_color_manual(values = PAL_ST, name = "Removal order") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  facet_wrap(~ d_lbl) +
  labs(title = "PC retention under core loss across dispersal distances",
       x = "Number of core patches removed", y = "PC retained (% of initial)") +
  theme_trbear_bar(base_size = 12)
tb_save_fig(p29b, "fig29b_robustness_bydist", w = 16, h = 6, subdir = FIG_SUBDIR)

tb_log_session()
tb_log("29_network_robustness DONE")
