## ============================================================================
## 12_resistance.R
## Project: TR_Bear_Connectivity
## Purpose: Convert habitat suitability rasters into landscape resistance for
##   UNICOR connectivity modelling. Formula follows Trainor et al. (2013),
##   parameterised per Shokri et al. (2021):
##
##       R(h) = 100 − 99 · ((1 − exp(−c·h)) / (1 − exp(−c))),  c = 4
##
##   h ∈ [0, 1] is habitat suitability; c controls curve shape (smaller c =
##   more linear; larger c = sharper resistance drop as h increases). c = 4
##   reflects bears' ability to use moderately suitable habitat during
##   dispersal. h = 1 → R = 1; h = 0 → R = 100.
##
##   Inputs:
##     - derived/hs_present/wmean.tif                       (present W_MEAN HS)
##     - derived/hs_gcm_avg/<period>_<ssp>.tif              (6 GCM-averaged HS)
##   Outputs:
##     - derived/resistance/present.tif
##     - derived/resistance/<period>_<ssp>.tif              (× 6)
##     - tables/12_resistance_summary.csv
##     - figures/12_resistance/
##         fig12a_resistance_present.png
##         fig12b_resistance_7panel.png   (present + 6 futures, log-scale)
##         fig12c_resistance_curve.png    (HS → R transfer function)
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
tb_log_init("12_resistance")

FIG_SUBDIR <- "12_resistance"

## ----------------------------------------------------------------------------
## 1) Transfer function
## ----------------------------------------------------------------------------
tb_log_section("Resistance transfer function")

.hs_to_R <- function(h,
                     c    = TB_RESIST_C,
                     rmin = TB_RESIST_MIN,
                     rmax = TB_RESIST_MAX) {
  rmax - (rmax - rmin) * (1 - exp(-c * h)) / (1 - exp(-c))
}

tb_log(sprintf("c = %.1f, range [%.0f, %.0f]", TB_RESIST_C, TB_RESIST_MIN, TB_RESIST_MAX))
tb_log(sprintf("R(h=1.0) = %.4f", .hs_to_R(1.0)))
tb_log(sprintf("R(h=0.7) = %.4f", .hs_to_R(0.7)))
tb_log(sprintf("R(h=0.5) = %.4f", .hs_to_R(0.5)))
tb_log(sprintf("R(h=0.3) = %.4f", .hs_to_R(0.3)))
tb_log(sprintf("R(h=0.1) = %.4f", .hs_to_R(0.1)))
tb_log(sprintf("R(h=0.0) = %.4f", .hs_to_R(0.0)))

## ----------------------------------------------------------------------------
## 2) Build resistance rasters: present + 6 GCM-averaged futures
## ----------------------------------------------------------------------------
tb_log_section("Build resistance rasters")

scen_list <- c(
  present = file.path(TB_OUT_HS_PRESENT, "wmean.tif"),
  setNames(
    file.path(TB_OUT_HS_AVG, sprintf("%s_%s.tif", rep(TB_PERIODS, each = length(TB_SSPS)),
                                                  rep(TB_SSPS,    times = length(TB_PERIODS)))),
    sprintf("%s_%s", rep(TB_PERIODS, each = length(TB_SSPS)),
                     rep(TB_SSPS,    times = length(TB_PERIODS)))
  )
)
tb_log(sprintf("scenarios queued: %d", length(scen_list)))

cell_area_km2 <- NULL
res_rasters   <- vector("list", length(scen_list))
summary_rows  <- vector("list", length(scen_list))

tb_tic()
for (i in seq_along(scen_list)) {
  k   <- names(scen_list)[i]
  src <- scen_list[[i]]
  if (!file.exists(src)) {
    tb_log(sprintf("SKIP %s — missing input: %s", k, src), "WARN")
    next
  }
  hs <- terra::rast(src)
  if (terra::nlyr(hs) > 1) hs <- hs[[1]]
  if (is.null(cell_area_km2)) cell_area_km2 <- prod(terra::res(hs)) / 1e6

  R <- terra::app(hs, .hs_to_R)
  names(R) <- "R"

  fn <- file.path(TB_OUT_RESISTANCE, sprintf("%s.tif", k))
  terra::writeRaster(R, fn, overwrite = TRUE, datatype = "FLT4S",
                     gdal = c("COMPRESS=DEFLATE","PREDICTOR=2","TILED=YES"))
  tb_log(sprintf("wrote %s", fn))

  res_rasters[[i]] <- R

  v <- terra::values(R, mat = FALSE); v <- v[!is.na(v)]
  summary_rows[[i]] <- data.frame(
    scenario  = k,
    n_cells   = length(v),
    R_mean    = mean(v),
    R_median  = median(v),
    R_q25     = quantile(v, 0.25, names = FALSE),
    R_q75     = quantile(v, 0.75, names = FALSE),
    R_max     = max(v),
    pct_easy_R_lt_10     = 100 * mean(v < 10),
    pct_moderate_10_50   = 100 * mean(v >= 10 & v < 50),
    pct_high_50_80       = 100 * mean(v >= 50 & v < 80),
    pct_barrier_ge_80    = 100 * mean(v >= 80),
    area_easy_km2        = sum(v < 10) * cell_area_km2,
    area_barrier_km2     = sum(v >= 80) * cell_area_km2
  )
}
tb_toc("resistance rasters built")

summary_df <- do.call(rbind, summary_rows)
names(res_rasters) <- names(scen_list)
tb_save_table(summary_df, "12_resistance_summary")

## ----------------------------------------------------------------------------
## FIGURES
## ----------------------------------------------------------------------------
tb_log_section("Figures")

world_sf <- tryCatch(
  rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") |>
    sf::st_transform(TB_CRS_PROJ),
  error = function(e) NULL)
tr_mask_shp <- file.path(TB_DATA_ROOT, "TR_mask", "TR_mask.shp")
tr_mask_sf <- if (file.exists(tr_mask_shp))
  sf::st_read(tr_mask_shp, quiet = TRUE) |> sf::st_transform(TB_CRS_PROJ) else NULL

template <- res_rasters[["present"]]
e <- terra::ext(template); ext_pad <- 30000
xlim_p <- c(e$xmin - ext_pad, e$xmax + ext_pad)
ylim_p <- c(e$ymin - ext_pad, e$ymax + ext_pad)

.tb_clip_to_tr <- function(r) {
  if (is.null(tr_mask_sf)) return(r)
  terra::mask(r, terra::vect(tr_mask_sf))
}

## ---- fig12a: present resistance (large, log-scale) --------------------------
tb_log_section("fig12a present")

R_pres <- .tb_clip_to_tr(res_rasters[["present"]])

p12a <- ggplot()
if (!is.null(world_sf)) p12a <- p12a +
  geom_sf(data = world_sf, fill = "#E8E8E8",
          color = "#7C8A93", linewidth = 0.4)
p12a <- p12a +
  tidyterra::geom_spatraster(data = R_pres, na.rm = TRUE) +
  scale_fill_viridis_c(option = "rocket", direction = -1, trans = "log10",
                       breaks = c(1, 5, 10, 25, 50, 80, 100),
                       labels = c("1", "5", "10", "25", "50", "80", "100"),
                       name = "Resistance\n(log scale)",
                       limits = c(1, 100),
                       na.value = "transparent",
                       guide = guide_colorbar(barwidth = 0.8,
                                              barheight = 12,
                                              ticks.colour = "black"))
if (!is.null(tr_mask_sf)) p12a <- p12a +
  geom_sf(data = tr_mask_sf, fill = NA,
          color = TB_COLOR_FRAME, linewidth = 0.5)
p12a <- p12a +
  coord_sf(xlim = xlim_p, ylim = ylim_p, datum = sf::st_crs(4326),
           expand = FALSE) +
  tb_map_decorations() +
  labs(title    = "Landscape resistance for brown bear movement across Türkiye for present-day",
       subtitle = NULL) +
  theme_trbear()
tb_save_fig(p12a, "fig12a_resistance_present", w = 13, h = 8.5, subdir = FIG_SUBDIR)

## ---- fig12b: 6-panel resistance (futures only, shared log scale) ------------
tb_log_section("fig12b 6-panel (futures)")

panel_order <- sprintf("%s_%s", rep(TB_PERIODS, each = length(TB_SSPS)),
                                  rep(TB_SSPS,    times = length(TB_PERIODS)))

.tb_panel_R <- function(r, title, show_legend = FALSE) {
  r_c <- .tb_clip_to_tr(r)
  p <- ggplot()
  if (!is.null(world_sf)) p <- p + geom_sf(data = world_sf,
                                            fill = "#E8E8E8",
                                            color = "#7C8A93",
                                            linewidth = 0.3)
  p <- p +
    tidyterra::geom_spatraster(data = r_c, na.rm = TRUE) +
    scale_fill_viridis_c(option = "rocket", direction = -1, trans = "log10",
                         limits = c(1, 100),
                         breaks = c(1, 5, 10, 25, 50, 80, 100),
                         labels = c("1", "5", "10", "25", "50", "80", "100"),
                         name = "Resistance\n(log scale)",
                         na.value = "transparent",
                         guide = guide_colorbar(barwidth = 1.0,
                                                barheight = 14,
                                                ticks.colour = "black"))
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

.title_for <- function(k) {
  if (k == "present") return("Present")
  parts <- strsplit(k, "_")[[1]]
  per <- paste(parts[1], parts[2], sep = "_")
  ssp <- parts[3]
  sprintf("%s — %s",
          TB_PERIOD_LABELS[per],
          TB_SSP_LABELS[ssp])
}

plots_b <- lapply(seq_along(panel_order), function(i) {
  k <- panel_order[i]
  r <- res_rasters[[k]]
  if (is.null(r)) return(patchwork::plot_spacer())
  .tb_panel_R(r, .title_for(k), show_legend = (i == length(panel_order)))
})

p12b <- patchwork::wrap_plots(plots_b, ncol = 3, nrow = 2, guides = "collect") +
  patchwork::plot_annotation(
    title    = "Landscape resistance throughout future scenarios",
    theme = theme(plot.title    = element_text(face = "bold", size = 18,
                                               color = TB_COLOR_FRAME)))
tb_save_fig(p12b, "fig12b_resistance_7panel", w = 22, h = 13, subdir = FIG_SUBDIR)

## ---- fig12c: transfer function curve ----------------------------------------
tb_log_section("fig12c transfer curve")

curve_df <- data.frame(HS = seq(0, 1, length.out = 501))
curve_df$R <- .hs_to_R(curve_df$HS)

p12c <- ggplot(curve_df, aes(x = HS, y = R)) +
  geom_area(fill = "#F08080", alpha = 0.35) +
  geom_line(color = "#9E2A2B", linewidth = 1.2) +
  geom_hline(yintercept = c(10, 50, 80),
             linetype = "dashed", color = "gray50", linewidth = 0.4) +
  annotate("text", x = 0.02, y = 12, label = "R = 10  (~easy)",   hjust = 0, size = 3.2, color = "gray35") +
  annotate("text", x = 0.02, y = 55, label = "R = 50  (high)",    hjust = 0, size = 3.2, color = "gray35") +
  annotate("text", x = 0.02, y = 88, label = "R = 80  (barrier)", hjust = 0, size = 3.2, color = "gray35") +
  scale_y_log10(breaks = c(1, 5, 10, 25, 50, 80, 100),
                limits = c(0.9, 120),
                expand = expansion(mult = c(0, 0.02))) +
  scale_x_continuous(breaks = seq(0, 1, 0.1), expand = c(0, 0)) +
  labs(title    = sprintf("Habitat suitability → resistance transfer function (c = %d)", TB_RESIST_C),
       subtitle = NULL,
       x = "Habitat suitability (W_MEAN)", y = "Landscape resistance (R)") +
  theme_trbear_bar()

tb_save_fig(p12c, "fig12c_resistance_curve", w = 11, h = 7, subdir = FIG_SUBDIR)

tb_log_session()
tb_log("12_resistance DONE")
