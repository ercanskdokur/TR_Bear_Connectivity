## ============================================================================
## 28_dpc_rankshift.R
## Project: TR_Bear_Connectivity
## Purpose: Does the IDENTITY and RANKING of the most important core patches
##   change under climate change?  ("Importance reshuffling")
##
##   The present-day network is dominated by a handful of cores (Pareto). As
##   climate shrinks/shifts suitable habitat northward, the relative dPC
##   importance of individual cores may re-rank: some southern cores collapse,
##   some northern cores rise. We track each PRESENT core forward into each of
##   the 6 future scenarios by spatial overlap and follow its dPC and rank.
##
## Method:
##   - 17_conefor.R already computed per-patch dPC for present + 6 futures at
##     6 dispersal distances (17_conefor_dpc_dii.csv). BUT patch_ids are
##     scenario-specific (terra::patches is run independently per scenario).
##   - We re-run terra::patches on the SAME binary rasters used by 17 (so ids
##     match the CSV exactly) and match every future patch to the present patch
##     it overlaps most (modal overlap). This yields a stable present-anchored
##     identity for each core across scenarios.
##   - A present core with no overlapping future patch (>= TB_PATCH_MIN_KM2) is
##     "Lost" in that scenario.
##
## Focal dispersal distance: 100 km (female-to-subadult dispersal; the distance
##   used throughout the connectivity figures). Table reports all distances.
##
## Outputs (tables/):
##   28_dpc_rankshift_long.csv   present_patch x scenario x d_km -> matched dPC,
##                                rank, status
##   28_dpc_rankshift_d100.csv   wide, d=100, top present cores
## Figures (figures/28_rankshift/):
##   fig28a_dpc_rank_bump.png    bump chart of rank across scenarios (top cores)
##   fig28b_dpc_value_traj.png   dPC value trajectory (log y) for top cores
## ============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("ggrepel", quietly = TRUE))
    install.packages("ggrepel", repos = "https://cloud.r-project.org")
  library(terra); library(sf); library(dplyr); library(tidyr)
  library(ggplot2); library(ggrepel)
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
tb_log_init("28_dpc_rankshift")

FIG_SUBDIR <- "28_rankshift"
FOCAL_D    <- 100
TOP_N      <- 10

scenarios <- c("present",
               sprintf("%s_%s", rep(TB_PERIODS, each = length(TB_SSPS)),
                                rep(TB_SSPS,    times = length(TB_PERIODS))))
scen_label <- c(
  "present"          = "Present",
  "2041_2070_ssp126" = "2070s SSP126",
  "2041_2070_ssp370" = "2070s SSP370",
  "2041_2070_ssp585" = "2070s SSP585",
  "2071_2100_ssp126" = "2100s SSP126",
  "2071_2100_ssp370" = "2100s SSP370",
  "2071_2100_ssp585" = "2100s SSP585")

.bin_path <- function(s) {
  if (s == "present") file.path(TB_OUT_HS_BINARY, "present_wmean.tif")
  else                file.path(TB_OUT_HS_BINARY, sprintf("future_%s.tif", s))
}

## ----------------------------------------------------------------------------
## 1) dPC table from script 17
## ----------------------------------------------------------------------------
tb_log_section("Read 35_costdist_dpc_dii.csv")
dpc <- read.csv(file.path(TB_OUT_TABLES, "35_costdist_dpc_dii.csv"))  # cost-distance
stopifnot(all(c("scenario","d_km","patch_id","dPC","dIIC","area_km2") %in% names(dpc)))

## Canonical core labels (C01..C93) from 34_core_crosswalk.R; fall back to raw id
.cw_file <- file.path(TB_OUT_TABLES, "34_core_crosswalk.csv")
core_cw <- if (file.exists(.cw_file)) read.csv(.cw_file)[, c("patch_id","core_id")] else
  data.frame(patch_id = unique(dpc$patch_id), core_id = paste0("id", unique(dpc$patch_id)))
core_of <- function(pid) {
  out <- core_cw$core_id[match(pid, core_cw$patch_id)]
  ifelse(is.na(out), paste0("id", pid), out)
}

## ----------------------------------------------------------------------------
## 2) Patch rasters (ids match the CSV because same terra::patches on same tif)
## ----------------------------------------------------------------------------
tb_log_section("Patch rasters + cross-scenario matching")

make_patch_raster <- function(s) {
  pp <- .bin_path(s)
  if (!file.exists(pp)) { tb_log(sprintf("[%s] binary missing", s), "WARN"); return(NULL) }
  r <- terra::rast(pp)
  terra::patches(r, directions = 8, zeroAsNA = TRUE)
}

pat_pres <- make_patch_raster("present")
names(pat_pres) <- "pres_id"

## present anchor table (id, area, centroid) at FOCAL_D
anchor <- dpc |>
  dplyr::filter(scenario == "present", d_km == FOCAL_D) |>
  dplyr::select(pres_id = patch_id, x, y, area_km2)
tb_log(sprintf("present anchors (cores >= %d km2): %d", TB_PATCH_MIN_KM2, nrow(anchor)))

## For a future scenario: modal overlap of each present patch onto future patches
match_to_present <- function(s) {
  if (s == "present") {
    return(data.frame(scenario = s, pres_id = anchor$pres_id,
                      fut_id = anchor$pres_id))
  }
  patf <- make_patch_raster(s)
  if (is.null(patf)) return(NULL)
  names(patf) <- "fut_id"
  if (!terra::compareGeom(pat_pres, patf, stopOnError = FALSE)) {
    patf <- terra::resample(patf, pat_pres, method = "near")
    names(patf) <- "fut_id"
  }
  ct <- terra::crosstab(c(pat_pres, patf), long = TRUE)
  names(ct) <- c("pres_id", "fut_id", "n")
  ct <- ct[!is.na(ct$pres_id) & !is.na(ct$fut_id) & ct$n > 0, ]
  if (!nrow(ct)) return(NULL)
  best <- ct |>
    dplyr::group_by(pres_id) |>
    dplyr::slice_max(n, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::transmute(scenario = s,
                     pres_id = as.integer(as.character(pres_id)),
                     fut_id  = as.integer(as.character(fut_id)))
  best
}

match_all <- do.call(rbind, lapply(scenarios, match_to_present))

## ----------------------------------------------------------------------------
## 3) Attach dPC/dIIC of the matched future patch, per scenario x distance
## ----------------------------------------------------------------------------
tb_log_section("Join matched dPC and rank")

DIST_ALL <- sort(unique(dpc$d_km))

long <- lapply(DIST_ALL, function(d) {
  base <- expand.grid(pres_id = anchor$pres_id, scenario = scenarios,
                      stringsAsFactors = FALSE)
  base <- dplyr::left_join(base, match_all, by = c("pres_id", "scenario"))
  ## join dPC of matched future patch (scenario-specific patch_id == fut_id)
  fut_vals <- dpc |>
    dplyr::filter(d_km == d) |>
    dplyr::select(scenario, fut_id = patch_id, dPC, dIIC, area_km2_fut = area_km2)
  out <- dplyr::left_join(base, fut_vals, by = c("scenario", "fut_id"))
  out$d_km <- d
  out$dPC  <- ifelse(is.na(out$dPC), 0, out$dPC)      # unmatched/lost -> 0
  out$dIIC <- ifelse(is.na(out$dIIC), 0, out$dIIC)
  out$status <- ifelse(is.na(out$fut_id), "Lost", "Present")
  out
}) |> dplyr::bind_rows()

## rank within scenario x distance (1 = most important); Lost cores get worst rank
long <- long |>
  dplyr::group_by(scenario, d_km) |>
  dplyr::mutate(rank = dplyr::min_rank(dplyr::desc(dPC))) |>
  dplyr::ungroup()

## present-rank for selecting & colouring top cores
pres_rank <- long |>
  dplyr::filter(scenario == "present", d_km == FOCAL_D) |>
  dplyr::transmute(pres_id, present_dPC = dPC, present_rank = rank,
                   core_id = core_of(pres_id))
long <- dplyr::left_join(long, pres_rank, by = "pres_id")

tb_save_table(long, "28_dpc_rankshift_long")

wide_d100 <- long |>
  dplyr::filter(d_km == FOCAL_D) |>
  dplyr::select(core_id, pres_id, present_rank, present_dPC, scenario, dPC, rank, status) |>
  tidyr::pivot_wider(names_from = scenario,
                     values_from = c(dPC, rank, status)) |>
  dplyr::arrange(present_rank)
tb_save_table(wide_d100, "28_dpc_rankshift_d100")

## ----------------------------------------------------------------------------
## 4) FIGURES (focal distance)
## ----------------------------------------------------------------------------
tb_log_section("Figures")

top_ids <- pres_rank |> dplyr::arrange(present_rank) |>
  dplyr::slice_head(n = TOP_N) |> dplyr::pull(pres_id)

## ordered core labels (C01..) by present rank for legend ordering
core_levels <- pres_rank |>
  dplyr::filter(pres_id %in% top_ids) |>
  dplyr::arrange(present_rank) |>
  dplyr::pull(core_id)
plot_df <- long |>
  dplyr::filter(d_km == FOCAL_D, pres_id %in% top_ids) |>
  dplyr::mutate(scen_lbl = factor(scen_label[scenario],
                                  levels = scen_label[scenarios]),
                core_lbl = factor(core_id, levels = core_levels))

pal_core <- setNames(
  colorRampPalette(c("#9E2A2B","#D55E00","#E69F00","#009E73","#0072B2","#3B0F70"))(length(top_ids)),
  core_levels)

## ---- fig28a: bump chart (rank) --------------------------------------------
end_lab <- plot_df |> dplyr::filter(scenario == "present")
p28a <- ggplot(plot_df, aes(scen_lbl, rank, group = core_lbl, color = core_lbl)) +
  geom_line(linewidth = 1.1, alpha = 0.85) +
  geom_point(aes(shape = status), size = 3) +
  scale_shape_manual(values = c("Present" = 16, "Lost" = 4), name = "Core fate") +
  scale_color_manual(values = pal_core, name = "Core (present rank)") +
  scale_y_reverse(breaks = scales::breaks_width(2)) +
  ggrepel::geom_text_repel(data = end_lab,
                           aes(label = core_id),
                           nudge_x = -0.35, direction = "y", size = 3.2,
                           segment.color = NA, show.legend = FALSE) +
  labs(title = sprintf("Re-ranking of core-patch importance under climate change (dPC, d = %d km)", FOCAL_D),
       subtitle = "Each line = one present-day core tracked forward by spatial overlap. ✕ = core lost (no suitable successor patch).",
       x = NULL, y = "Importance rank (1 = most important)") +
  theme_trbear_bar(base_size = 12) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
tb_save_fig(p28a, "fig28a_dpc_rank_bump", w = 13, h = 8, subdir = FIG_SUBDIR)

## ---- fig28b: dPC value trajectory (log) -----------------------------------
plot_df_v <- plot_df |> dplyr::mutate(dPC_plot = ifelse(dPC <= 0, NA, dPC))
p28b <- ggplot(plot_df_v, aes(scen_lbl, dPC_plot, group = core_lbl, color = core_lbl)) +
  geom_line(linewidth = 1.0, alpha = 0.85) +
  geom_point(aes(shape = status), size = 2.6) +
  scale_shape_manual(values = c("Present" = 16, "Lost" = 4), name = "Core fate") +
  scale_color_manual(values = pal_core, name = "Core (present rank)") +
  scale_y_log10() +
  labs(title = sprintf("Absolute dPC importance of top-%d present cores across scenarios", TOP_N),
       subtitle = sprintf("d = %d km. Log scale. Dropping lines = cores whose connectivity contribution collapses.", FOCAL_D),
       x = NULL, y = "dPC (%, log scale)") +
  theme_trbear_bar(base_size = 12) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
tb_save_fig(p28b, "fig28b_dpc_value_traj", w = 13, h = 8, subdir = FIG_SUBDIR)

tb_log_session()
tb_log("28_dpc_rankshift DONE")
