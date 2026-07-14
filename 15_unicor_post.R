## ============================================================================
## 15_unicor_post.R
## Project: TR_Bear_Connectivity
## Purpose: Post-process UNICOR outputs (7 scenarios) into publication figures.
##
## Improvements over earlier versions:
##   1. SHARED normalization — all 7 KDE rasters divided by the GLOBAL maximum
##      (max across all scenarios), so panel colours are directly comparable.
##   2. Sqrt-stretched + percentile-binned visual scale — corridors stand out
##      instead of being washed out by the long zero-dominated tail.
##   3. Grid alignment — futures are resampled onto present grid before
##      computing delta (UNICOR pads extent slightly differently per run).
##   4. Top-X% corridor highlight figures (binary, 7-panel + overlay).
##
## UNICOR output extensions (Arc/Info ASCII grid, custom suffix):
##   .kdepaths       — kernel-density of cumulative LCPs (the corridor surface)
##   .addedpaths.txt — cumulative LCP count per cell
##   .levels         — categorical (quantile-binned) KDE
##   .cdmatrix.csv   — cost-distance matrix between source points
## ============================================================================

suppressPackageStartupMessages({
  library(terra); library(sf); library(ggplot2); library(dplyr); library(tidyr)
  library(tidyterra); library(patchwork); library(rnaturalearth); library(scales)
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
tb_log_init("15_unicor_post")

FIG_SUBDIR <- "15_unicor"
TOP_PCT    <- 0.05   ## top-X% corridor threshold (fraction of cells)

## ----------------------------------------------------------------------------
## 1) Locate + load UNICOR outputs
## ----------------------------------------------------------------------------
tb_log_section("Locate UNICOR outputs")

scenarios <- c("present",
               sprintf("%s_%s", rep(TB_PERIODS, each = length(TB_SSPS)),
                                rep(TB_SSPS,    times = length(TB_PERIODS))))

.find_one <- function(dir, pattern) {
  hits <- list.files(dir, pattern = pattern, full.names = TRUE,
                     ignore.case = TRUE)
  if (length(hits)) hits[1] else NA_character_
}

scen_paths <- lapply(scenarios, function(s) {
  res_dir <- file.path(TB_OUT_UNICOR_DIR, s, "results")
  if (!dir.exists(res_dir)) {
    tb_log(sprintf("MISSING results dir: %s", res_dir), "WARN")
    return(list(scen = s, kde = NA, path = NA, lev = NA, cdm = NA))
  }
  list(
    scen = s,
    kde  = .find_one(res_dir, "\\.kdepaths$"),
    path = .find_one(res_dir, "\\.addedpaths\\.txt$"),
    lev  = .find_one(res_dir, "\\.levels$"),
    cdm  = .find_one(res_dir, "cdmatrix.*\\.csv$")
  )
})
names(scen_paths) <- scenarios

for (s in scenarios) {
  sp <- scen_paths[[s]]
  tb_log(sprintf("[%s] kde=%s | path=%s", s,
                 basename(sp$kde  %||% "NA"),
                 basename(sp$path %||% "NA")))
}

.read_aaigrid <- function(src) {
  if (is.na(src) || !file.exists(src)) return(NULL)
  proxy <- tempfile(fileext = ".asc")
  file.copy(src, proxy, overwrite = TRUE)
  r <- terra::rast(proxy)
  if (is.na(terra::crs(r)) || terra::crs(r) == "") terra::crs(r) <- TB_CRS_PROJ
  ## UNICOR sometimes leaves tiny-negative numerical artefacts — clip to 0.
  r <- terra::ifel(r < 0, 0, r)
  names(r) <- "value"
  r
}

tb_log_section("Load rasters")
kde_raw  <- lapply(scen_paths, function(sp) .read_aaigrid(sp$kde))
path_raw <- lapply(scen_paths, function(sp) .read_aaigrid(sp$path))
names(kde_raw)  <- scenarios
names(path_raw) <- scenarios

## ----------------------------------------------------------------------------
## 2) Grid alignment — resample everything onto present's grid
## ----------------------------------------------------------------------------
tb_log_section("Grid alignment")

ref_kde <- kde_raw[["present"]]
if (is.null(ref_kde)) {
  tb_log("present KDE missing — cannot align grids", "ERROR")
  tb_log_session(); quit(status = 1)
}

for (s in scenarios) {
  if (s == "present") next
  if (is.null(kde_raw[[s]])) next
  if (!terra::compareGeom(ref_kde, kde_raw[[s]], stopOnError = FALSE)) {
    tb_log(sprintf("resample %s KDE → present grid", s))
    kde_raw[[s]] <- terra::resample(kde_raw[[s]], ref_kde, method = "bilinear")
    if (!is.null(path_raw[[s]])) {
      path_raw[[s]] <- terra::resample(path_raw[[s]], ref_kde, method = "bilinear")
    }
  }
}

## ----------------------------------------------------------------------------
## 3) Global max → SHARED normalization
## ----------------------------------------------------------------------------
tb_log_section("Shared normalization")

global_max <- max(sapply(kde_raw, function(r) {
  if (is.null(r)) NA_real_ else max(terra::values(r, mat = FALSE), na.rm = TRUE)
}), na.rm = TRUE)
tb_log(sprintf("global KDE max across all 7 scenarios: %.4g", global_max))

kde_norm <- lapply(kde_raw, function(r) {
  if (is.null(r)) return(NULL)
  rn <- r / global_max
  names(rn) <- "value"
  rn
})

## ----------------------------------------------------------------------------
## 4) Summary table — raw + normalized + corridor top-5%
## ----------------------------------------------------------------------------
tb_log_section("Summary")

cell_area_km2 <- prod(terra::res(ref_kde)) / 1e6

global_top_thr <- {
  v_all <- unlist(lapply(kde_raw, function(r) {
    if (is.null(r)) NULL else terra::values(r, mat = FALSE)
  }))
  v_all <- v_all[!is.na(v_all) & v_all > 0]
  quantile(v_all, 1 - TOP_PCT, names = FALSE)
}
tb_log(sprintf("top-%.0f%% threshold (global): %.4g", 100 * TOP_PCT, global_top_thr))

summary_rows <- lapply(scenarios, function(s) {
  rr <- kde_raw[[s]]
  if (is.null(rr)) return(data.frame(scenario = s))
  v  <- terra::values(rr, mat = FALSE); v <- v[!is.na(v)]
  vn <- v / global_max
  n_tot <- length(v)
  n_corr <- sum(v >= global_top_thr)
  data.frame(
    scenario          = s,
    n_cells           = n_tot,
    raw_max           = max(v),
    raw_mean          = mean(v),
    raw_median        = median(v),
    norm_mean         = mean(vn),
    norm_median       = median(vn),
    pct_corridor_top5 = 100 * n_corr / n_tot,
    area_corridor_km2 = n_corr * cell_area_km2
  )
})
summary_df <- do.call(rbind, summary_rows)
tb_save_table(summary_df, "15_connectivity_summary")

## ----------------------------------------------------------------------------
## 5) Basemap, mask, source points for plotting
## ----------------------------------------------------------------------------
sources_gpkg <- file.path(TB_OUT_VECTORS, "sources.gpkg")
sources_sf <- if (file.exists(sources_gpkg))
  sf::st_read(sources_gpkg, quiet = TRUE) else NULL

world_sf <- tryCatch(
  rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") |>
    sf::st_transform(TB_CRS_PROJ),
  error = function(e) NULL)
tr_mask_shp <- file.path(TB_DATA_ROOT, "TR_mask", "TR_mask.shp")
tr_mask_sf <- if (file.exists(tr_mask_shp))
  sf::st_read(tr_mask_shp, quiet = TRUE) |> sf::st_transform(TB_CRS_PROJ) else NULL

e <- terra::ext(ref_kde); ext_pad <- 30000
xlim_p <- c(e$xmin - ext_pad, e$xmax + ext_pad)
ylim_p <- c(e$ymin - ext_pad, e$ymax + ext_pad)

.clip_to_tr <- function(r) {
  if (is.null(tr_mask_sf) || is.null(r)) return(r)
  terra::mask(r, terra::vect(tr_mask_sf))
}

.title_for <- function(k) {
  if (k == "present") return("Present")
  parts <- strsplit(k, "_")[[1]]
  per <- paste(parts[1], parts[2], sep = "_")
  ssp <- parts[3]
  sprintf("%s — %s", TB_PERIOD_LABELS[per], TB_SSP_LABELS[ssp])
}

## Sqrt transform so low-but-positive corridors are visible without log-NaN.
## Limits run from a small floor up to 1 (shared normalized scale).
.kde_scale <- function(name = "Connectivity\n(KDE, shared-max norm)") {
  scale_fill_viridis_c(
    option = "magma", direction = -1,
    trans  = scales::pseudo_log_trans(sigma = 0.01),
    na.value = "transparent",
    limits = c(0, 1),
    breaks = c(0, 0.01, 0.05, 0.1, 0.25, 0.5, 1),
    labels = c("0", "0.01", "0.05", "0.10", "0.25", "0.50", "1.00"),
    name   = name)
}

## ----------------------------------------------------------------------------
## 6) FIGURES
## ----------------------------------------------------------------------------

## ---- fig15a: present connectivity main map (continuous, normalized) ---------
tb_log_section("fig15a present")

if (!is.null(kde_norm[["present"]])) {
  r <- .clip_to_tr(kde_norm[["present"]])
  p15a <- ggplot()
  if (!is.null(world_sf)) p15a <- p15a +
    geom_sf(data = world_sf, fill = "#E8E8E8",
            color = "#7C8A93", linewidth = 0.4)
  p15a <- p15a +
    tidyterra::geom_spatraster(data = r, na.rm = TRUE) +
    .kde_scale()
  if (!is.null(tr_mask_sf)) p15a <- p15a +
    geom_sf(data = tr_mask_sf, fill = NA,
            color = TB_COLOR_FRAME, linewidth = 0.5)
  if (!is.null(sources_sf)) p15a <- p15a +
    geom_sf(data = sources_sf, fill = "white", color = "#000000",
            shape = 21, size = 1.4, stroke = 0.4)
  p15a <- p15a +
    coord_sf(xlim = xlim_p, ylim = ylim_p, datum = sf::st_crs(4326),
             expand = FALSE) +
    tb_map_decorations() +
    labs(title    = "Brown bear connectivity across Türkiye for present-day",
         subtitle = sprintf(
           "UNICOR KDE of cumulative least-cost paths between %d source patch centroids (≥ %d km²).",
           nrow(sources_sf %||% data.frame()), TB_PATCH_MIN_KM2)) +
    theme_trbear()
  tb_save_fig(p15a, "fig15a_present_connectivity", w = 14, h = 9, subdir = FIG_SUBDIR)
}

## ---- fig15b: 7-panel SHARED normalization -----------------------------------
tb_log_section("fig15b 7-panel KDE")

.panel_kde <- function(r, title, show_legend = FALSE) {
  rc <- .clip_to_tr(r)
  p <- ggplot()
  if (!is.null(world_sf)) p <- p + geom_sf(data = world_sf,
                                            fill = "#E8E8E8",
                                            color = "#7C8A93",
                                            linewidth = 0.3)
  p <- p +
    tidyterra::geom_spatraster(data = rc, na.rm = TRUE) +
    .kde_scale(name = "Connectivity")
  if (!is.null(tr_mask_sf)) p <- p +
    geom_sf(data = tr_mask_sf, fill = NA, color = TB_COLOR_FRAME,
            linewidth = 0.35)
  p <- p +
    coord_sf(xlim = xlim_p, ylim = ylim_p, datum = sf::st_crs(4326),
             expand = FALSE) +
    labs(title = title) +
    theme_trbear(base_size = 11) +
    theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold",
                                    margin = ggplot2::margin(b = 3)),
          plot.margin = ggplot2::margin(3, 3, 3, 3))
  if (!show_legend) p <- p + theme(legend.position = "none")
  p
}

plots_b <- lapply(seq_along(scenarios), function(i) {
  r <- kde_norm[[scenarios[i]]]
  if (is.null(r)) return(patchwork::plot_spacer())
  .panel_kde(r, .title_for(scenarios[i]),
             show_legend = (i == length(scenarios)))
})

## Layout: 2 rows × 4 cols with present in top-left, futures filling the rest
plots_b8 <- c(plots_b[1], plots_b[2:7], list(patchwork::plot_spacer()))
p15b <- patchwork::wrap_plots(plots_b8, ncol = 4, nrow = 2, guides = "collect") +
  patchwork::plot_annotation(
    title    = "Brown bear connectivity – present-day and future scenarios",
    theme = theme(plot.title    = element_text(face = "bold", size = 18,
                                               color = TB_COLOR_FRAME)))
tb_save_fig(p15b, "fig15b_kde_7panel", w = 24, h = 12, subdir = FIG_SUBDIR)

## ---- fig15c: 6-panel delta KDE (future − present), ABSOLUTE values ----------
tb_log_section("fig15c delta")

if (!is.null(kde_raw[["present"]])) {
  fut_scens <- scenarios[scenarios != "present"]

  ## Symmetric limit from the largest absolute delta value seen
  delta_max <- max(sapply(fut_scens, function(s) {
    if (is.null(kde_raw[[s]])) return(0)
    d <- kde_raw[[s]] - kde_raw[["present"]]
    max(abs(terra::values(d, mat = FALSE)), na.rm = TRUE)
  }), na.rm = TRUE)
  tb_log(sprintf("delta abs-max across 6 futures: %.4g", delta_max))

  plots_c <- lapply(seq_along(fut_scens), function(i) {
    s <- fut_scens[i]
    if (is.null(kde_raw[[s]])) return(patchwork::plot_spacer())
    d <- kde_raw[[s]] - kde_raw[["present"]]
    names(d) <- "delta"
    dc <- .clip_to_tr(d)
    p <- ggplot()
    if (!is.null(world_sf)) p <- p + geom_sf(data = world_sf,
                                              fill = "#E8E8E8",
                                              color = "#7C8A93",
                                              linewidth = 0.3)
    p <- p +
      tidyterra::geom_spatraster(data = dc, na.rm = TRUE) +
      scale_fill_gradient2(low = "#D55E00", mid = "#F7F7F7", high = "#009E73",
                           midpoint = 0,
                           limits = c(-delta_max, delta_max),
                           oob = scales::squish, na.value = "transparent",
                           name = "Δ KDE\n(absolute,\nfuture − present)")
    if (!is.null(tr_mask_sf)) p <- p +
      geom_sf(data = tr_mask_sf, fill = NA, color = TB_COLOR_FRAME,
              linewidth = 0.35)
    p <- p +
      coord_sf(xlim = xlim_p, ylim = ylim_p, datum = sf::st_crs(4326),
               expand = FALSE) +
      labs(title = .title_for(s)) +
      theme_trbear(base_size = 11) +
      theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
            plot.margin = ggplot2::margin(3, 3, 3, 3))
    if (i != length(fut_scens)) p <- p + theme(legend.position = "none")
    p
  })

  p15c <- patchwork::wrap_plots(plots_c, ncol = 3, nrow = 2, guides = "collect") +
    patchwork::plot_annotation(
      title    = "Connectivity change (present vs future scenarios)",
      subtitle = "Δ KDE in absolute units (future – present, after grid alignment)",
      theme = theme(plot.title    = element_text(face = "bold", size = 18,
                                                 color = TB_COLOR_FRAME),
                    plot.subtitle = element_text(size = 12, color = TB_COLOR_AXIS)))
  tb_save_fig(p15c, "fig15c_connectivity_delta_6panel", w = 22, h = 13, subdir = FIG_SUBDIR)
}

## ---- fig15d: present cumulative-paths raster (corridor density) -------------
tb_log_section("fig15d present paths")

if (!is.null(path_raw[["present"]])) {
  rp <- path_raw[["present"]]
  v  <- terra::values(rp, mat = FALSE); v <- v[!is.na(v) & v > 0]
  rp_n <- if (length(v) > 0) rp / max(v) else rp * 0
  names(rp_n) <- "value"
  rpc <- .clip_to_tr(rp_n)

  p15d <- ggplot()
  if (!is.null(world_sf)) p15d <- p15d +
    geom_sf(data = world_sf, fill = "#E8E8E8",
            color = "#7C8A93", linewidth = 0.4)
  p15d <- p15d +
    tidyterra::geom_spatraster(data = rpc, na.rm = TRUE) +
    scale_fill_viridis_c(option = "viridis", direction = 1,
                         trans = scales::pseudo_log_trans(sigma = 0.005),
                         limits = c(0, 1),
                         breaks = c(0, 0.01, 0.05, 0.25, 1),
                         labels = c("0","0.01","0.05","0.25","1"),
                         na.value = "transparent",
                         name = "Path density\n(0–1)")
  if (!is.null(tr_mask_sf)) p15d <- p15d +
    geom_sf(data = tr_mask_sf, fill = NA,
            color = TB_COLOR_FRAME, linewidth = 0.5)
  if (!is.null(sources_sf)) p15d <- p15d +
    geom_sf(data = sources_sf, fill = "white", color = "#000000",
            shape = 21, size = 1.4, stroke = 0.4)
  p15d <- p15d +
    coord_sf(xlim = xlim_p, ylim = ylim_p, datum = sf::st_crs(4326),
             expand = FALSE) +
    tb_map_decorations() +
    labs(title    = "Cumulative least-cost paths across Türkiye for present day",
         subtitle = NULL) +
    theme_trbear()
  tb_save_fig(p15d, "fig15d_present_path", w = 14, h = 9, subdir = FIG_SUBDIR)
}

## ---- fig15e: top-5% corridor binary, present --------------------------------
tb_log_section("fig15e present top-5% corridor")

if (!is.null(kde_raw[["present"]])) {
  pres_top <- kde_raw[["present"]] >= global_top_thr
  pres_top <- terra::as.int(pres_top)
  pres_top_fac <- terra::as.factor(pres_top)
  levels(pres_top_fac) <- data.frame(
    id    = c(0, 1),
    class = c("Matrix", sprintf("Top-%.0f%% corridor", 100 * TOP_PCT)))
  names(pres_top_fac) <- "class"
  pres_top_c <- .clip_to_tr(pres_top_fac)

  p15e <- ggplot()
  if (!is.null(world_sf)) p15e <- p15e +
    geom_sf(data = world_sf, fill = "#E8E8E8",
            color = "#7C8A93", linewidth = 0.4)
  p15e <- p15e +
    tidyterra::geom_spatraster(data = pres_top_c, na.rm = TRUE) +
    scale_fill_manual(values = setNames(c("#F2F2F2", "#9E2A2B"),
                                         c("Matrix", sprintf("Top-%.0f%% corridor", 100 * TOP_PCT))),
                      na.translate = FALSE, name = NULL)
  if (!is.null(tr_mask_sf)) p15e <- p15e +
    geom_sf(data = tr_mask_sf, fill = NA,
            color = TB_COLOR_FRAME, linewidth = 0.5)
  if (!is.null(sources_sf)) p15e <- p15e +
    geom_sf(data = sources_sf, fill = "white", color = "#000000",
            shape = 21, size = 1.4, stroke = 0.4)
  p15e <- p15e +
    coord_sf(xlim = xlim_p, ylim = ylim_p, datum = sf::st_crs(4326),
             expand = FALSE) +
    tb_map_decorations() +
    labs(title    = sprintf("Top-%.0f%% connectivity corridors across Türkiye for present-day", 100 * TOP_PCT),
         subtitle = sprintf("Cells with KDE ≥ global %.0fth percentile (threshold = %.3g).",
                            100 * (1 - TOP_PCT), global_top_thr)) +
    theme_trbear()
  tb_save_fig(p15e, "fig15e_top5_corridor_present", w = 14, h = 9, subdir = FIG_SUBDIR)
}

## ---- fig15f: top-5% corridor 7-panel (corridor persistence) -----------------
tb_log_section("fig15f top-5% corridor 7-panel")

.panel_top <- function(r, title) {
  bin <- terra::as.int(r >= global_top_thr)
  bin <- terra::as.factor(bin)
  levels(bin) <- data.frame(id = c(0, 1),
                            class = c("Matrix", "Top-5% corridor"))
  names(bin) <- "class"
  bc <- .clip_to_tr(bin)
  p <- ggplot()
  if (!is.null(world_sf)) p <- p + geom_sf(data = world_sf,
                                            fill = "#E8E8E8",
                                            color = "#7C8A93",
                                            linewidth = 0.3)
  p <- p +
    tidyterra::geom_spatraster(data = bc, na.rm = TRUE) +
    scale_fill_manual(values = setNames(c("#F2F2F2", "#9E2A2B"),
                                         c("Matrix", "Top-5% corridor")),
                      na.translate = FALSE, name = NULL)
  if (!is.null(tr_mask_sf)) p <- p +
    geom_sf(data = tr_mask_sf, fill = NA,
            color = TB_COLOR_FRAME, linewidth = 0.3)
  p +
    coord_sf(xlim = xlim_p, ylim = ylim_p, datum = sf::st_crs(4326),
             expand = FALSE) +
    labs(title = title) +
    theme_trbear(base_size = 11) +
    theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold",
                                    margin = ggplot2::margin(b = 3)),
          plot.margin = ggplot2::margin(3, 3, 3, 3),
          legend.position = "none")
}

plots_f <- lapply(scenarios, function(s) {
  if (is.null(kde_raw[[s]])) return(patchwork::plot_spacer())
  .panel_top(kde_raw[[s]], .title_for(s))
})

plots_f8 <- c(plots_f[1], plots_f[2:7], list(patchwork::plot_spacer()))
p15f <- patchwork::wrap_plots(plots_f8, ncol = 4, nrow = 2) +
  patchwork::plot_annotation(
    title    = sprintf("Top-%.0f%% corridors – present vs future scenarios",
                       100 * TOP_PCT),
    theme = theme(plot.title    = element_text(face = "bold", size = 18,
                                               color = TB_COLOR_FRAME)))
tb_save_fig(p15f, "fig15f_top5_corridor_7panel", w = 24, h = 12, subdir = FIG_SUBDIR)

## ----------------------------------------------------------------------------
## 7) Persist the shared global threshold + max for downstream scripts
## ----------------------------------------------------------------------------
tb_log_section("Persist metadata")

tb_save_table(
  data.frame(metric = c("global_kde_max", "top_pct_threshold", "top_pct_fraction"),
             value  = c(global_max,        global_top_thr,    TOP_PCT)),
  "15_global_kde_meta")

tb_log_session()
tb_log("15_unicor_post DONE")
