## ============================================================================
## 16_landscapemetrics.R
## Project: TR_Bear_Connectivity
## Purpose: Class-level landscape metrics on UNICOR corridor masks.
##   Two corridor thresholds are evaluated per scenario, with GLOBAL thresholds
##   pooled across all 7 KDE rasters so the masks are directly comparable:
##     - Top-5% corridor backbone   (~13.05 KDE units; matches 15_unicor_post)
##     - Top-1% corridor strict core (re-computed from pooled values)
##   Habitat baseline: present_wmean binary HS (suitable / unsuitable).
##
## Metrics (per scenario × threshold):
##   CA          Class Area (km²) — total corridor extent
##   NP          Number of Patches
##   AREA_MN     Mean patch area (ha)
##   AREA_CV     Patch-size coefficient of variation (heterogeneity)
##   LPI         Largest patch index (% of landscape)
##   ENN_MN      Mean Euclidean nearest-neighbour distance (m)
##   COHESION    Patch cohesion index — physical connectedness
##   AI          Aggregation index — like-adjacency
##   CLUMPY      Clumpiness deviation from random
##
## Outputs (tables/):
##   16_lsm_corridor_top5.csv     class metrics on 5% mask, per scenario
##   16_lsm_corridor_top1.csv     class metrics on 1% mask, per scenario
##   16_lsm_habitat_present.csv   class metrics on present binary HS (baseline)
##   16_lsm_summary_long.csv      long-format combined table for plotting
##
## Figures (figures/16_lsm/):
##   fig16a_metric_trends.png     6-panel: CA, NP, AREA_MN, ENN_MN, COHESION, AI
##                                vs scenario (5% vs 1% lines)
##   fig16b_size_distribution.png patch-area histogram per scenario (present-5%)
## ============================================================================

suppressPackageStartupMessages({
  for (.p in c("landscapemetrics", "patchwork")) {
    if (!requireNamespace(.p, quietly = TRUE)) {
      install.packages(.p, repos = "https://cloud.r-project.org")
    }
  }
  library(terra); library(sf); library(ggplot2); library(dplyr); library(tidyr)
  library(landscapemetrics); library(patchwork)
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
tb_log_init("16_landscapemetrics")

FIG_SUBDIR <- "16_lsm"
TOP_PCTS   <- c(top5 = 0.05, top1 = 0.01)

scenarios <- c("present",
               sprintf("%s_%s", rep(TB_PERIODS, each = length(TB_SSPS)),
                                rep(TB_SSPS,    times = length(TB_PERIODS))))

## ----------------------------------------------------------------------------
## 1) Locate + load aligned KDE rasters (re-uses 15_unicor_post loader logic)
## ----------------------------------------------------------------------------
tb_log_section("Load UNICOR KDE rasters")

.find_one <- function(dir, pattern) {
  hits <- list.files(dir, pattern = pattern, full.names = TRUE,
                     ignore.case = TRUE)
  if (length(hits)) hits[1] else NA_character_
}

.read_aaigrid <- function(src) {
  if (is.na(src) || !file.exists(src)) return(NULL)
  proxy <- tempfile(fileext = ".asc")
  file.copy(src, proxy, overwrite = TRUE)
  r <- terra::rast(proxy)
  if (is.na(terra::crs(r)) || terra::crs(r) == "") terra::crs(r) <- TB_CRS_PROJ
  r <- terra::ifel(r < 0, 0, r)
  names(r) <- "value"
  r
}

kde_raw <- lapply(scenarios, function(s) {
  res_dir <- file.path(TB_OUT_UNICOR_DIR, s, "results")
  .read_aaigrid(.find_one(res_dir, "\\.kdepaths$"))
})
names(kde_raw) <- scenarios

if (is.null(kde_raw[["present"]])) {
  tb_log("present KDE missing — abort", "ERROR"); tb_log_session(); quit(status = 1)
}

## Align to present
ref_kde <- kde_raw[["present"]]
for (s in scenarios) {
  if (s == "present") next
  if (is.null(kde_raw[[s]])) next
  if (!terra::compareGeom(ref_kde, kde_raw[[s]], stopOnError = FALSE)) {
    tb_log(sprintf("resample %s → present grid", s))
    kde_raw[[s]] <- terra::resample(kde_raw[[s]], ref_kde, method = "bilinear")
  }
}

## ---- Clip to Türkiye mask --------------------------------------------------
tr_mask_shp <- file.path(TB_DATA_ROOT, "TR_mask", "TR_mask.shp")
tr_mask_sf  <- if (file.exists(tr_mask_shp))
  sf::st_read(tr_mask_shp, quiet = TRUE) |> sf::st_transform(TB_CRS_PROJ) else NULL
.clip <- function(r) if (is.null(tr_mask_sf)) r else terra::mask(r, terra::vect(tr_mask_sf))

kde_raw <- lapply(kde_raw, .clip)

## ----------------------------------------------------------------------------
## 2) Global thresholds (pooled across all 7 scenarios)
## ----------------------------------------------------------------------------
tb_log_section("Global thresholds")

v_all <- unlist(lapply(kde_raw, function(r) {
  if (is.null(r)) NULL else terra::values(r, mat = FALSE)
}))
v_all <- v_all[!is.na(v_all) & v_all > 0]

thresholds <- sapply(TOP_PCTS, function(p) {
  quantile(v_all, 1 - p, names = FALSE)
})
tb_log(sprintf("top-5%% threshold = %.4g | top-1%% threshold = %.4g",
               thresholds["top5"], thresholds["top1"]))
write.csv(data.frame(label   = names(thresholds),
                     fraction = TOP_PCTS,
                     threshold = thresholds),
          file.path(TB_OUT_TABLES, "16_corridor_thresholds.csv"),
          row.names = FALSE)

## ----------------------------------------------------------------------------
## 3) Compute class metrics for each scenario × threshold
## ----------------------------------------------------------------------------
tb_log_section("Class metrics — corridor masks")

LSM_CLASS <- c("lsm_c_ca", "lsm_c_np", "lsm_c_pland", "lsm_c_area_mn",
               "lsm_c_area_cv", "lsm_c_lpi", "lsm_c_enn_mn",
               "lsm_c_cohesion", "lsm_c_ai", "lsm_c_clumpy")

.metrics_for_mask <- function(r_bin, label_class = "corridor") {
  ## r_bin: numeric SpatRaster with values 0/1/NA, where 1 = focal class
  res <- landscapemetrics::calculate_lsm(
    r_bin, what = LSM_CLASS, directions = 8, verbose = FALSE)
  ## Keep only focal class (value == 1)
  res <- dplyr::filter(res, class == 1)
  res$class_label <- label_class
  res
}

corridor_rows <- list()
for (s in scenarios) {
  r <- kde_raw[[s]]
  if (is.null(r)) next
  for (lab in names(TOP_PCTS)) {
    thr <- thresholds[[lab]]
    bin <- terra::ifel(r >= thr, 1, 0)
    bin <- terra::mask(bin, r)
    m <- .metrics_for_mask(bin, label_class = lab)
    m$scenario <- s
    m$threshold <- lab
    corridor_rows[[paste(s, lab, sep = "_")]] <- m
    tb_log(sprintf("[%s | %s] CA=%.0f NP=%d ENN_MN=%.0f COHESION=%.2f",
                   s, lab,
                   m$value[m$metric == "ca"],
                   m$value[m$metric == "np"],
                   m$value[m$metric == "enn_mn"],
                   m$value[m$metric == "cohesion"]))
  }
}
corridor_df <- do.call(rbind, corridor_rows)

corridor_wide <- corridor_df |>
  select(scenario, threshold, metric, value) |>
  tidyr::pivot_wider(names_from = metric, values_from = value)

tb_save_table(filter(corridor_wide, threshold == "top5"), "16_lsm_corridor_top5")
tb_save_table(filter(corridor_wide, threshold == "top1"), "16_lsm_corridor_top1")
tb_save_table(corridor_df, "16_lsm_summary_long")

## ----------------------------------------------------------------------------
## 4) Habitat-patch baseline (present binary HS)
## ----------------------------------------------------------------------------
tb_log_section("Class metrics — present binary HS")

pres_bin_path <- file.path(TB_OUT_HS_BINARY, "present_wmean.tif")
if (file.exists(pres_bin_path)) {
  pres_bin <- terra::rast(pres_bin_path) |> .clip()
  hab_m <- .metrics_for_mask(pres_bin, label_class = "habitat")
  hab_m$scenario <- "present"
  hab_m$threshold <- "habitat"
  tb_save_table(hab_m, "16_lsm_habitat_present")
  tb_log(sprintf("[habitat baseline] CA=%.0f NP=%d AREA_MN=%.1f COHESION=%.2f",
                 hab_m$value[hab_m$metric == "ca"],
                 hab_m$value[hab_m$metric == "np"],
                 hab_m$value[hab_m$metric == "area_mn"],
                 hab_m$value[hab_m$metric == "cohesion"]))
} else {
  tb_log("present_wmean.tif missing — habitat baseline skipped", "WARN")
}

## ----------------------------------------------------------------------------
## 5) FIGURES
## ----------------------------------------------------------------------------
tb_log_section("Figures")

## Order scenarios for plotting
scen_order <- scenarios
scen_label <- c(
  "present"          = "Present",
  "2041_2070_ssp126" = "2070s\nSSP126",
  "2041_2070_ssp370" = "2070s\nSSP370",
  "2041_2070_ssp585" = "2070s\nSSP585",
  "2071_2100_ssp126" = "2100s\nSSP126",
  "2071_2100_ssp370" = "2100s\nSSP370",
  "2071_2100_ssp585" = "2100s\nSSP585"
)

corridor_df$scenario <- factor(corridor_df$scenario, levels = scen_order)
corridor_df$threshold <- factor(corridor_df$threshold,
                                levels = c("top5", "top1"),
                                labels = c("Top-5% (backbone)",
                                           "Top-1% (strict core)"))

PAL_THR <- c("Top-5% (backbone)"    = "#9E2A2B",
             "Top-1% (strict core)" = "#1F3A93")

.metric_panel <- function(metric_name, ytitle, log_y = FALSE) {
  d <- filter(corridor_df, metric == metric_name)
  p <- ggplot(d, aes(scenario, value, color = threshold,
                     group = threshold)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.4) +
    scale_color_manual(values = PAL_THR, name = "Corridor mask") +
    scale_x_discrete(labels = scen_label) +
    labs(title = toupper(metric_name), y = ytitle, x = NULL) +
    theme_trbear_bar(base_size = 11) +
    theme(axis.text.x = element_text(size = 9),
          plot.title = element_text(size = 12, face = "bold"))
  if (log_y) p <- p + scale_y_log10(labels = scales::label_comma())
  p
}

p_ca       <- .metric_panel("ca",       "Class area (ha)", log_y = TRUE)
p_np       <- .metric_panel("np",       "Number of patches")
p_pland    <- .metric_panel("pland",    "PLAND — % of landscape")
p_lpi      <- .metric_panel("lpi",      "LPI — largest patch index (%)")
p_area_mn  <- .metric_panel("area_mn",  "Mean patch area (ha)", log_y = TRUE)
p_enn_mn   <- .metric_panel("enn_mn",   "Mean NN distance (m)")
p_cohesion <- .metric_panel("cohesion", "Cohesion")
p_ai       <- .metric_panel("ai",       "Aggregation index")
p_clumpy   <- .metric_panel("clumpy",   "CLUMPY — clumpiness")

## 3 × 3 grid: structure | extent | configuration
fig16a <- (p_ca   + p_pland   + p_lpi)    /
          (p_np   + p_area_mn + p_enn_mn) /
          (p_cohesion + p_ai  + p_clumpy) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(
    title    = "Corridor landscape metrics – present vs future scenarios",
    subtitle = "Class-level metrics on UNICOR top-5% (backbone) and top-1% (strict core) KDE masks.\nRow 1: corridor extent (CA, PLAND, LPI).  Row 2: patch count, mean size, isolation (NP, AREA_MN, ENN_MN).  Row 3: spatial configuration (COHESION, AI, CLUMPY).",
    theme = theme(plot.title    = element_text(face = "bold", size = 16,
                                                color = TB_COLOR_FRAME),
                  plot.subtitle = element_text(size = 10.5, color = TB_COLOR_AXIS)))
tb_save_fig(fig16a, "fig16a_metric_trends", w = 18, h = 14, subdir = FIG_SUBDIR)

## ---- fig16b: patch-area distribution — present, both thresholds -------------
.patch_areas <- function(r_bin) {
  r_bin <- terra::ifel(r_bin == 1, 1, NA)
  p <- terra::patches(r_bin, directions = 8, zeroAsNA = TRUE)
  fr <- terra::freq(p) |> as.data.frame()
  if (!nrow(fr)) return(numeric())
  cell_km2 <- prod(terra::res(r_bin)) / 1e6
  fr$count * cell_km2
}

r_pres <- kde_raw[["present"]]
sizes_top5 <- .patch_areas(terra::ifel(r_pres >= thresholds[["top5"]], 1, 0))
sizes_top1 <- .patch_areas(terra::ifel(r_pres >= thresholds[["top1"]], 1, 0))

sizes_df <- rbind(
  data.frame(area_km2 = sizes_top5,
             threshold = "Top-5% (backbone)"),
  data.frame(area_km2 = sizes_top1,
             threshold = "Top-1% (strict core)"))
sizes_df$threshold <- factor(sizes_df$threshold,
                             levels = c("Top-5% (backbone)",
                                        "Top-1% (strict core)"))

fig16b <- ggplot(sizes_df, aes(area_km2, fill = threshold)) +
  geom_histogram(bins = 40, color = "white", alpha = 0.85) +
  scale_fill_manual(values = PAL_THR, name = "Corridor mask") +
  scale_x_log10(labels = scales::label_comma()) +
  facet_wrap(~ threshold, scales = "free_y") +
  labs(title    = "Corridor patch-size distribution — present scenario",
       subtitle = sprintf("Top-5%% n = %d patches | Top-1%% n = %d patches",
                          length(sizes_top5), length(sizes_top1)),
       x = "Patch area (km², log scale)", y = "Number of patches") +
  theme_trbear_bar(base_size = 12) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"))
tb_save_fig(fig16b, "fig16b_size_distribution", w = 14, h = 7, subdir = FIG_SUBDIR)

tb_log_session()
tb_log("16_landscapemetrics DONE")
