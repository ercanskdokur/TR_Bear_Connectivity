## ============================================================================
## 17_conefor.R
## Project: TR_Bear_Connectivity
## Purpose: Compute Conefor-style connectivity indices (PC, IIC) and per-patch
##   importance (dPC, dIIC) for the present scenario and each of the 6 future
##   scenarios, at 6 dispersal threshold distances:
##       50, 100, 150, 200, 300, 400 km   (covering female → long-range male)
##
## Approach (pure R, no external Conefor binary required):
##   - For each scenario, read binary HS raster (suitable / unsuitable)
##   - Identify connected suitable patches with terra::patches; keep all
##     patches above TB_PATCH_MIN_KM2 (consistent with UNICOR source selection)
##   - Get patch centroid (x, y) and area a_i
##   - For each distance threshold d:
##       p_ij = exp(-alpha · dist_ij), with alpha = -ln(0.5) / d
##         (so p = 0.5 at d, p ≈ 1 for short distances, p → 0 for long ones)
##       PC  = Σ_i Σ_j a_i a_j P*_ij  / A_L²
##         where P*_ij is the maximum-product probability along any path
##         (computed via Dijkstra on −log p_ij weights)
##       IIC = Σ_i Σ_j a_i a_j / (1 + nl_ij) / A_L²
##         where nl_ij is the number of links along the shortest path in the
##         unweighted graph of all pairs with dist_ij ≤ d
##       A_L = total landscape area (TR landmass, km²)
##       dPC_k = 100 · (PC - PC_minus_k) / PC   (patch importance)
##       dIIC_k = 100 · (IIC - IIC_minus_k) / IIC
##
##   Reference: Saura & Pascual-Hortal (2007) Landsc Urban Plan 83:91-103;
##              Saura & Torné (2009) Environ Model Softw 24:135-139.
##
## Outputs (tables/):
##   17_conefor_indices.csv      scenario × distance → PC, IIC, n_patches,
##                                                     total habitat km²
##   17_conefor_dpc_dii.csv      scenario × distance × patch_id → a_km2,
##                                                                 dPC, dIIC
##
## Figures (figures/17_conefor/):
##   fig17a_pc_vs_distance.png   PC vs distance, lines per scenario
##   fig17b_iic_vs_distance.png  IIC vs distance, lines per scenario
##   fig17c_top10_dpc.png        Top-10 dPC patches map (present, d = 100 km)
## ============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("igraph", quietly = TRUE)) {
    install.packages("igraph", repos = "https://cloud.r-project.org")
  }
  library(terra); library(sf); library(ggplot2); library(dplyr); library(tidyr)
  library(igraph); library(patchwork); library(rnaturalearth); library(tidyterra)
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
tb_log_init("17_conefor")

FIG_SUBDIR     <- "17_conefor"
DIST_KM        <- c(50, 100, 150, 200, 300, 400)
PROB_AT_THRESH <- 0.5

scenarios <- c("present",
               sprintf("%s_%s", rep(TB_PERIODS, each = length(TB_SSPS)),
                                rep(TB_SSPS,    times = length(TB_PERIODS))))

.bin_path <- function(s) {
  if (s == "present") file.path(TB_OUT_HS_BINARY, "present_wmean.tif")
  else                file.path(TB_OUT_HS_BINARY, sprintf("future_%s.tif", s))
}

## ----------------------------------------------------------------------------
## 1) Total landscape area A_L (km²) — TR landmass
## ----------------------------------------------------------------------------
tb_log_section("Landscape area")

tr_mask_shp <- file.path(TB_DATA_ROOT, "TR_mask", "TR_mask.shp")
tr_mask_sf  <- if (file.exists(tr_mask_shp))
  sf::st_read(tr_mask_shp, quiet = TRUE) |> sf::st_transform(TB_CRS_PROJ) else NULL

if (!is.null(tr_mask_sf)) {
  AL_km2 <- as.numeric(sum(sf::st_area(tr_mask_sf))) / 1e6
  tb_log(sprintf("A_L (TR landmass) = %.0f km²", AL_km2))
} else {
  ## Fallback: use raster extent
  r0 <- terra::rast(.bin_path("present"))
  AL_km2 <- prod(terra::res(r0)) * terra::ncell(r0) / 1e6
  tb_log(sprintf("A_L (raster extent fallback) = %.0f km²", AL_km2), "WARN")
}

## ----------------------------------------------------------------------------
## 2) Per-scenario patch extraction
## ----------------------------------------------------------------------------
tb_log_section("Patch extraction")

extract_patches <- function(s) {
  pp <- .bin_path(s)
  if (!file.exists(pp)) {
    tb_log(sprintf("[%s] binary HS missing", s), "WARN"); return(NULL)
  }
  r <- terra::rast(pp); names(r) <- "suit"
  cell_km2 <- prod(terra::res(r)) / 1e6
  pat <- terra::patches(r, directions = 8, zeroAsNA = TRUE)
  fr  <- terra::freq(pat) |> as.data.frame() |>
    dplyr::rename(patch_id = value, n_cells = count) |>
    dplyr::mutate(area_km2 = n_cells * cell_km2) |>
    dplyr::filter(!is.na(patch_id), area_km2 >= TB_PATCH_MIN_KM2)
  if (!nrow(fr)) return(NULL)

  ## Centroid of each kept patch
  out <- lapply(fr$patch_id, function(i) {
    pts <- terra::as.points(terra::ifel(pat == i, 1, NA))
    if (nrow(pts) == 0) return(NULL)
    crd <- terra::crds(pts)
    data.frame(patch_id = i,
               x = mean(crd[, 1]), y = mean(crd[, 2]),
               n_cells = nrow(crd))
  })
  out <- do.call(rbind, out)
  out$area_km2 <- out$n_cells * cell_km2
  out$scenario <- s
  out
}

patches_by_scen <- lapply(scenarios, extract_patches)
names(patches_by_scen) <- scenarios

for (s in scenarios) {
  ps <- patches_by_scen[[s]]
  if (is.null(ps)) next
  tb_log(sprintf("[%s] n_patches=%d | total_km2=%.0f | A_L_used=%.0f",
                 s, nrow(ps), sum(ps$area_km2), AL_km2))
}

## ----------------------------------------------------------------------------
## 3) Index computation — pure R via igraph
## ----------------------------------------------------------------------------
tb_log_section("PC / IIC / dPC / dIIC")

## Maximum-product shortest paths (PC): use -log(p) as edge weight
compute_pc <- function(coords, area, d_km, AL_km2) {
  n <- nrow(coords)
  if (n < 2) return(list(PC = sum(area^2) / AL_km2^2,
                         dPC = rep(0, n)))
  D <- as.matrix(stats::dist(coords))            # metres
  D_km <- D / 1000
  alpha <- -log(PROB_AT_THRESH) / d_km
  W <- -log(pmax(exp(-alpha * D_km), .Machine$double.eps))  # -log(p_ij)
  diag(W) <- 0
  g <- igraph::graph_from_adjacency_matrix(W, mode = "undirected", weighted = TRUE)
  ## shortest paths give min(-log p) = -log(P*)
  L <- igraph::distances(g, weights = igraph::E(g)$weight)
  Pstar <- exp(-L)
  PC_num <- sum(outer(area, area) * Pstar)
  PC <- PC_num / AL_km2^2

  dPC <- numeric(n)
  for (k in seq_len(n)) {
    sub <- (seq_len(n))[-k]
    g_k <- igraph::delete_vertices(g, k)
    L_k <- igraph::distances(g_k)
    Pstar_k <- exp(-L_k)
    num_k <- sum(outer(area[sub], area[sub]) * Pstar_k)
    PC_k <- num_k / AL_km2^2
    dPC[k] <- if (PC > 0) 100 * (PC - PC_k) / PC else 0
  }
  list(PC = PC, dPC = dPC)
}

## Integral index of connectivity (IIC): binary graph at d threshold,
## edge if dist_ij ≤ d. nl_ij = unweighted shortest path length.
compute_iic <- function(coords, area, d_km, AL_km2) {
  n <- nrow(coords)
  if (n < 2) return(list(IIC = sum(area^2) / AL_km2^2,
                         dIIC = rep(0, n)))
  D_km <- as.matrix(stats::dist(coords)) / 1000
  Adj <- (D_km <= d_km) & (D_km > 0)
  g <- igraph::graph_from_adjacency_matrix(Adj, mode = "undirected")
  nl <- igraph::distances(g)
  ## i==j: nl=0, so a_i^2 / 1 contributes
  IIC_num <- sum(outer(area, area) / (1 + nl), na.rm = TRUE)
  IIC <- IIC_num / AL_km2^2

  dIIC <- numeric(n)
  for (k in seq_len(n)) {
    sub <- (seq_len(n))[-k]
    g_k <- igraph::delete_vertices(g, k)
    nl_k <- igraph::distances(g_k)
    num_k <- sum(outer(area[sub], area[sub]) / (1 + nl_k), na.rm = TRUE)
    IIC_k <- num_k / AL_km2^2
    dIIC[k] <- if (IIC > 0) 100 * (IIC - IIC_k) / IIC else 0
  }
  list(IIC = IIC, dIIC = dIIC)
}

idx_rows <- list()
imp_rows <- list()
for (s in scenarios) {
  ps <- patches_by_scen[[s]]
  if (is.null(ps) || nrow(ps) < 1) next
  coords <- as.matrix(ps[, c("x", "y")])
  area   <- ps$area_km2

  for (d in DIST_KM) {
    tb_log(sprintf("[%s | d=%d km] n=%d patches", s, d, nrow(ps)))
    pc_res  <- compute_pc(coords, area, d, AL_km2)
    iic_res <- compute_iic(coords, area, d, AL_km2)

    idx_rows[[paste(s, d, sep = "_")]] <- data.frame(
      scenario     = s,
      d_km         = d,
      n_patches    = nrow(ps),
      habitat_km2  = sum(area),
      PC           = pc_res$PC,
      IIC          = iic_res$IIC)

    imp_rows[[paste(s, d, sep = "_")]] <- data.frame(
      scenario  = s,
      d_km      = d,
      patch_id  = ps$patch_id,
      x         = ps$x, y = ps$y,
      area_km2  = ps$area_km2,
      dPC       = pc_res$dPC,
      dIIC      = iic_res$dIIC)
  }
}
idx_df <- do.call(rbind, idx_rows)
imp_df <- do.call(rbind, imp_rows)
tb_save_table(idx_df, "17_conefor_indices")
tb_save_table(imp_df, "17_conefor_dpc_dii")

## ----------------------------------------------------------------------------
## 4) FIGURES
## ----------------------------------------------------------------------------
tb_log_section("Figures")

scen_label <- c(
  "present"          = "Present",
  "2041_2070_ssp126" = "2070s SSP126",
  "2041_2070_ssp370" = "2070s SSP370",
  "2041_2070_ssp585" = "2070s SSP585",
  "2071_2100_ssp126" = "2100s SSP126",
  "2071_2100_ssp370" = "2100s SSP370",
  "2071_2100_ssp585" = "2100s SSP585")

PAL_SCEN <- c(
  "Present"      = "#000000",
  "2070s SSP126" = "#56B4E9",
  "2070s SSP370" = "#E69F00",
  "2070s SSP585" = "#D55E00",
  "2100s SSP126" = "#0072B2",
  "2100s SSP370" = "#CC79A7",
  "2100s SSP585" = "#9E2A2B")

idx_df$scen_lbl <- factor(scen_label[idx_df$scenario],
                          levels = scen_label[scenarios])

p17a <- ggplot(idx_df, aes(d_km, PC, color = scen_lbl, group = scen_lbl)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.4) +
  scale_color_manual(values = PAL_SCEN, name = "Scenario") +
  scale_x_continuous(breaks = DIST_KM) +
  labs(title    = "Probability of Connectivity (PC) across dispersal distances",
       x = "Dispersal threshold distance (km, p = 0.5 at threshold)",
       y = "PC (dimensionless)") +
  theme_trbear_bar(base_size = 12)
tb_save_fig(p17a, "fig17a_pc_vs_distance", w = 12, h = 7, subdir = FIG_SUBDIR)

p17b <- ggplot(idx_df, aes(d_km, IIC, color = scen_lbl, group = scen_lbl)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.4) +
  scale_color_manual(values = PAL_SCEN, name = "Scenario") +
  scale_x_continuous(breaks = DIST_KM) +
  labs(title    = "Integral Index of Connectivity (IIC) across dispersal distances",
       x = "Dispersal threshold distance (km, link cutoff)",
       y = "IIC (dimensionless)") +
  theme_trbear_bar(base_size = 12)
tb_save_fig(p17b, "fig17b_iic_vs_distance", w = 12, h = 7, subdir = FIG_SUBDIR)

## ---- fig17c: Top-10 dPC patches map at d = 100 km, present ------------------
present_imp <- imp_df |>
  filter(scenario == "present", d_km == 100) |>
  arrange(desc(dPC))

if (nrow(present_imp) >= 1) {
  pres_bin <- terra::rast(.bin_path("present"))
  if (!is.null(tr_mask_sf))
    pres_bin <- terra::mask(pres_bin, terra::vect(tr_mask_sf))
  pres_bin_fac <- terra::as.factor(terra::as.int(pres_bin))
  levels(pres_bin_fac) <- data.frame(id = c(0, 1),
                                      class = c("Unsuitable", "Suitable"))
  names(pres_bin_fac) <- "class"

  top10 <- present_imp[1:min(10, nrow(present_imp)), ]
  top10$rank <- seq_len(nrow(top10))
  top10_sf <- sf::st_as_sf(top10, coords = c("x", "y"), crs = TB_CRS_PROJ)

  world_sf <- tryCatch(
    rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") |>
      sf::st_transform(TB_CRS_PROJ),
    error = function(e) NULL)

  e <- terra::ext(pres_bin); pad <- 30000
  xl <- c(e$xmin - pad, e$xmax + pad)
  yl <- c(e$ymin - pad, e$ymax + pad)

  p17c <- ggplot()
  if (!is.null(world_sf)) p17c <- p17c +
    geom_sf(data = world_sf, fill = "#E8E8E8",
            color = "#7C8A93", linewidth = 0.4)
  p17c <- p17c +
    tidyterra::geom_spatraster(data = pres_bin_fac, na.rm = TRUE) +
    scale_fill_manual(values = TB_PAL_BINARY, na.translate = FALSE,
                      name = "Habitat")
  if (!is.null(tr_mask_sf)) p17c <- p17c +
    geom_sf(data = tr_mask_sf, fill = NA,
            color = TB_COLOR_FRAME, linewidth = 0.5)
  p17c <- p17c +
    geom_sf(data = top10_sf, aes(size = dPC),
            fill = "#9E2A2B", color = "black",
            shape = 21, stroke = 0.6, alpha = 0.9) +
    scale_size_continuous(name = "dPC (%)",
                          range = c(3, 9),
                          breaks = pretty(top10$dPC, 5)) +
    coord_sf(xlim = xl, ylim = yl, datum = sf::st_crs(4326), expand = FALSE) +
    tb_map_decorations() +
    labs(title = "Top-10 patches by dPC importance — present, d = 100 km") +
    theme_trbear()
  tb_save_fig(p17c, "fig17c_top10_dpc", w = 14, h = 9, subdir = FIG_SUBDIR)
}

tb_log_session()
tb_log("17_conefor DONE")
