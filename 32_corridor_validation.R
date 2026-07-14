## ============================================================================
## 32_corridor_validation.R
## Project: TR_Bear_Connectivity
## Purpose: INDEPENDENT validation of the predicted corridor network. The
##   corridors were derived purely from habitat-suitability -> resistance ->
##   least-cost paths; they never "saw" the occurrence or conflict data. If they
##   capture real bear movement, then (a) independent presence records and
##   (b) movement-related conflicts (esp. road accidents) should fall inside the
##   top-5% corridor far more than expected by chance.
##
## Method:
##   - Present top-5% corridor mask from the UNICOR KDE surface.
##   - Observed: fraction of conflict points (21_conflict_clean.csv) and
##     presence points (presence_clean.gpkg) inside the corridor.
##   - Null model (B = 999): resample the same number of points among
##       (A) all TR land cells, and
##       (B) suitable-habitat cells only (stricter — controls for the fact that
##           both bears and corridors concentrate in good habitat).
##     -> p-value and enrichment ratio (observed / null-mean).
##   - By conflict activity type: binomial enrichment vs the suitable-habitat
##     null expectation; road accidents are predicted to be the most
##     corridor-associated (movement mortality, not settlement foraging).
##
## Outputs (tables/):
##   32_corridor_validation.csv        group -> n, obs_frac, null mean/CI, enrich, p
##   32_validation_by_activity.csv     activity -> n, obs_frac, enrich, p
## Figures (figures/32_validate/):
##   fig32a_validation_null.png        null histograms + observed lines
##   fig32b_validation_by_activity.png enrichment by activity type
## ============================================================================

suppressPackageStartupMessages({
  library(terra); library(sf); library(dplyr); library(tidyr); library(ggplot2)
  library(patchwork)
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
tb_log_init("32_corridor_validation")

FIG_SUBDIR <- "32_validate"
TOP_PCT    <- 0.05
B_NULL     <- 999L
set.seed(42)

## ----------------------------------------------------------------------------
## 1) Present corridor mask
## ----------------------------------------------------------------------------
tb_log_section("Corridor mask")
.read_aaigrid <- function(src) {
  proxy <- tempfile(fileext = ".asc"); file.copy(src, proxy, overwrite = TRUE)
  r <- terra::rast(proxy)
  if (is.na(terra::crs(r)) || terra::crs(r) == "") terra::crs(r) <- TB_CRS_PROJ
  r <- terra::ifel(r < 0, 0, r); names(r) <- "value"; r
}
kde_src <- list.files(file.path(TB_OUT_UNICOR_DIR, "present", "results"),
                      pattern = "\\.kdepaths$", full.names = TRUE)[1]
kde <- .read_aaigrid(kde_src)
vv  <- terra::values(kde)[, 1]
pos <- which(!is.na(vv) & vv > 0)
thr <- stats::quantile(vv[pos], 1 - TOP_PCT, names = FALSE)
corr <- terra::ifel(kde >= thr, 1, 0); names(corr) <- "corr"
tb_log(sprintf("corridor threshold=%.4g | corridor cells=%d", thr, sum(vv >= thr, na.rm = TRUE)))

## land + suitable masks aligned to corridor grid
land <- terra::rast(file.path(TB_OUT_RASTERS, "tr_landmask.tif"))
if (!terra::compareGeom(corr, land, stopOnError = FALSE))
  land <- terra::resample(land, corr, method = "near")
suit <- terra::rast(file.path(TB_OUT_HS_BINARY, "present_wmean.tif"))
if (!terra::compareGeom(corr, suit, stopOnError = FALSE))
  suit <- terra::resample(suit, corr, method = "near")

corr_v <- terra::values(corr)[, 1]
land_cells <- which(!is.na(terra::values(land)[, 1]) & terra::values(land)[, 1] >= 1 & !is.na(corr_v))
suit_cells <- which(terra::values(suit)[, 1] == 1 & !is.na(corr_v))
exp_p_land <- mean(corr_v[land_cells] == 1)
exp_p_suit <- mean(corr_v[suit_cells] == 1)
tb_log(sprintf("expected corridor fraction: land=%.4f | suitable=%.4f", exp_p_land, exp_p_suit))

## ----------------------------------------------------------------------------
## 2) Observed points
## ----------------------------------------------------------------------------
tb_log_section("Observed points")
conf <- read.csv(file.path(TB_OUT_TABLES, "21_conflict_clean.csv"))
conf_sf <- sf::st_as_sf(conf, coords = c("Long","Lat"), crs = TB_CRS_WGS) |>
  sf::st_transform(TB_CRS_PROJ)
conf_xy <- sf::st_coordinates(conf_sf)

pres_gpkg <- file.path(TB_OUT_VECTORS, "presence_clean.gpkg")
pres_xy <- if (file.exists(pres_gpkg)) {
  p <- sf::st_read(pres_gpkg, quiet = TRUE) |> sf::st_transform(TB_CRS_PROJ)
  sf::st_coordinates(p)
} else NULL

frac_in <- function(xy) {
  e <- terra::extract(corr, xy)[, 1]
  mean(e == 1, na.rm = TRUE)
}

## ----------------------------------------------------------------------------
## 3) Resampling null
## ----------------------------------------------------------------------------
tb_log_section("Resampling null model")
null_frac <- function(n, domain_cells) {
  vapply(seq_len(B_NULL), function(b) {
    cc <- sample(domain_cells, n, replace = FALSE)
    mean(corr_v[cc] == 1)
  }, numeric(1))
}

groups <- list(
  list(name = "Conflict points", xy = conf_xy),
  if (!is.null(pres_xy)) list(name = "Presence points", xy = pres_xy) else NULL)
groups <- groups[!sapply(groups, is.null)]

val_rows <- list(); null_long <- list()
for (gp in groups) {
  n_g <- nrow(gp$xy); obs <- frac_in(gp$xy)
  for (dom in c("land","suitable")) {
    dc <- if (dom == "land") land_cells else suit_cells
    nd <- null_frac(n_g, dc)
    p  <- (1 + sum(nd >= obs)) / (B_NULL + 1)
    val_rows[[paste(gp$name, dom)]] <- data.frame(
      group = gp$name, null_domain = dom, n = n_g,
      obs_frac = obs, null_mean = mean(nd),
      null_lo = quantile(nd, 0.025), null_hi = quantile(nd, 0.975),
      enrichment = obs / mean(nd), p_value = p)
    null_long[[paste(gp$name, dom)]] <- data.frame(
      group = gp$name, null_domain = dom, null_frac = nd, obs_frac = obs)
  }
}
val_df <- dplyr::bind_rows(val_rows)
tb_save_table(val_df, "32_corridor_validation")
null_df <- dplyr::bind_rows(null_long)

## ----------------------------------------------------------------------------
## 4) By activity (binomial vs suitable-habitat expectation)
## ----------------------------------------------------------------------------
tb_log_section("By activity")
conf$in_corr <- terra::extract(corr, conf_xy)[, 1]
act_df <- conf |>
  dplyr::filter(!is.na(in_corr)) |>
  dplyr::group_by(activity) |>
  dplyr::summarise(n = dplyr::n(), k = sum(in_corr == 1), .groups = "drop") |>
  dplyr::mutate(
    obs_frac   = k / n,
    expected   = exp_p_suit,
    enrichment = obs_frac / exp_p_suit,
    p_value    = mapply(function(k, n) stats::binom.test(k, n, exp_p_suit,
                                                         alternative = "greater")$p.value, k, n)) |>
  dplyr::arrange(dplyr::desc(enrichment))
tb_save_table(act_df, "32_validation_by_activity")

## ----------------------------------------------------------------------------
## 5) FIGURES
## ----------------------------------------------------------------------------
tb_log_section("Figures")
lab_df <- val_df |>
  dplyr::mutate(txt = sprintf("obs = %.1f%%\n%.1f× expected\np %s",
                              100 * obs_frac, enrichment,
                              ifelse(p_value < 0.001, "< 0.001", sprintf("= %.3f", p_value))))
p32a <- ggplot(null_df, aes(null_frac)) +
  geom_histogram(bins = 40, fill = "#BBD4E6", color = "white", linewidth = 0.1) +
  geom_vline(data = val_df, aes(xintercept = obs_frac),
             color = "#9E2A2B", linewidth = 1.1) +
  geom_text(data = lab_df, aes(x = obs_frac, y = Inf, label = txt),
            hjust = 1.05, vjust = 1.3, size = 3, color = "#9E2A2B") +
  facet_grid(group ~ null_domain, scales = "free",
             labeller = labeller(null_domain = c(land = "Null: all TR land",
                                                  suitable = "Null: suitable habitat only"))) +
  scale_x_continuous(labels = scales::percent) +
  labs(title = "Independent validation of predicted corridors",
       subtitle = "Histogram = fraction of randomised points falling in the top-5% corridor; red line = observed.",
       x = "Fraction of points inside top-5% corridor", y = "Null replicates") +
  theme_trbear_bar(base_size = 12)
tb_save_fig(p32a, "fig32a_validation_null", w = 14, h = 8, subdir = FIG_SUBDIR)

## ---- fig32b: enrichment by activity ---------------------------------------
act_plot <- act_df |>
  dplyr::mutate(activity = factor(activity, levels = rev(activity)),
                sig = ifelse(p_value < 0.05, "p < 0.05", "n.s."))
p32b <- ggplot(act_plot, aes(enrichment, activity, fill = sig)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 1, linetype = 2, color = "gray40") +
  geom_text(aes(label = sprintf("%.1f× (n=%d)", enrichment, n)),
            hjust = -0.1, size = 3.2, color = TB_COLOR_AXIS) +
  scale_fill_manual(values = c("p < 0.05" = "#9E2A2B", "n.s." = "#BDBDBD"), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = "Corridor association of conflicts by activity type",
       subtitle = "Enrichment = observed / expected (suitable-habitat null). > 1 = over-represented inside corridors.",
       x = "Enrichment (× expected)", y = NULL) +
  theme_trbear_bar(base_size = 12)
tb_save_fig(p32b, "fig32b_validation_by_activity", w = 13, h = 7, subdir = FIG_SUBDIR)

tb_log_session()
tb_log("32_corridor_validation DONE")
