## ============================================================================
## 23_tables_master.R
## Project: TR_Bear_Connectivity
## Purpose: Master summary tables consolidating outputs from
##   scripts 06–21. Produces:
##     (1) ONE wide CSV: scenario × all metrics (24 columns)
##     (2) FOUR thematic CSVs (manuscript Table 1–4 candidates):
##          T1 — Habitat suitability (HS km², gain/loss)
##          T2 — Corridor topology (UNICOR area, landscapemetrics, PC, IIC)
##          T3 — Protected-area coverage (per-layer + summary, in/out km²)
##          T4 — Road conflict (paved km in corridor, pinch %)
##
## Inputs (from tables/):
##   08_present_summary.csv, 09_future_each_summary.csv,
##   10_gcm_avg_summary.csv, 11_gainloss_summary.csv,
##   12_resistance_summary.csv, 15_connectivity_summary.csv,
##   16_lsm_corridor_top5.csv, 16_lsm_corridor_top1.csv,
##   17_conefor_indices.csv, 17_conefor_dpc_dii.csv,
##   18_pa_overlay_summary.csv, 18_pa_by_layer.csv,
##   19_roads_overlay.csv,
##   21_xtab_bioregion_activity.csv, 21_xtab_pa_activity.csv,
##   21_road_dist_by_activity.csv
##
## Outputs (tables/):
##   23_master_wide.csv
##   23_T1_habitat.csv
##   23_T2_corridor.csv
##   23_T3_pa_coverage.csv
##   23_T4_roads.csv
## ============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr)
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
tb_log_init("23_tables_master")

scenarios <- c("present",
               sprintf("%s_%s", rep(TB_PERIODS, each = length(TB_SSPS)),
                                rep(TB_SSPS,    times = length(TB_PERIODS))))

scen_label <- c(
  "present"          = "Present",
  "2041_2070_ssp126" = "2070s SSP1-2.6",
  "2041_2070_ssp370" = "2070s SSP3-7.0",
  "2041_2070_ssp585" = "2070s SSP5-8.5",
  "2071_2100_ssp126" = "2100s SSP1-2.6",
  "2071_2100_ssp370" = "2100s SSP3-7.0",
  "2071_2100_ssp585" = "2100s SSP5-8.5")

.read_t <- function(name) {
  fn <- file.path(TB_OUT_TABLES, name)
  if (!file.exists(fn)) {
    tb_log(sprintf("MISSING table: %s", fn), "WARN"); return(NULL)
  }
  readr::read_csv(fn, show_col_types = FALSE)
}

## ---- Load all inputs --------------------------------------------------------
tb_log_section("Load inputs")

pa_sum   <- .read_t("18_pa_overlay_summary.csv")
pa_layer <- .read_t("18_pa_by_layer.csv")
roads    <- .read_t("19_roads_overlay.csv")
res_sum  <- .read_t("12_resistance_summary.csv")
con_sum  <- .read_t("15_connectivity_summary.csv")
lsm5     <- .read_t("16_lsm_corridor_top5.csv")
lsm1     <- .read_t("16_lsm_corridor_top1.csv")
con_idx  <- .read_t("17_conefor_indices.csv")
gainloss <- .read_t("11_gainloss_summary.csv")

## Conefor at d = 100 km is the headline (matches manuscript)
con_idx_100 <- if (!is.null(con_idx))
  con_idx |> dplyr::filter(d_km == 100) |>
    dplyr::select(scenario, n_patches, habitat_km2,
                  PC_100km = PC, IIC_100km = IIC) else NULL

## ---- T1: Habitat ------------------------------------------------------------
tb_log_section("T1 habitat")

T1 <- data.frame(scenario = scenarios)
T1$scenario_label <- scen_label[T1$scenario]
if (!is.null(pa_sum))
  T1 <- T1 |> dplyr::left_join(
    pa_sum |> dplyr::select(scenario, hs_total_km2),
    by = "scenario")
if (!is.null(gainloss))
  T1 <- T1 |> dplyr::left_join(
    gainloss |> dplyr::select(any_of(c("scenario","gain_km2","loss_km2",
                                        "stable_km2","gain_pct","loss_pct"))),
    by = "scenario")
T1 <- T1 |>
  dplyr::mutate(
    pct_of_TR = round(100 * hs_total_km2 / sum(hs_total_km2[scenario == "present"]), 1))
tb_save_table(T1, "23_T1_habitat")

## ---- T2: Corridor topology --------------------------------------------------
tb_log_section("T2 corridor")

T2 <- data.frame(scenario = scenarios)
T2$scenario_label <- scen_label[T2$scenario]
if (!is.null(con_sum))
  T2 <- T2 |> dplyr::left_join(
    con_sum |> dplyr::select(scenario,
                              area_corridor_km2,
                              pct_corridor_top5,
                              raw_max_kde = raw_max),
    by = "scenario")
if (!is.null(lsm5))
  T2 <- T2 |> dplyr::left_join(
    lsm5 |> dplyr::transmute(scenario,
                              top5_NP   = np,
                              top5_AREA_MN = area_mn,
                              top5_ENN_MN  = enn_mn,
                              top5_COHESION = cohesion),
    by = "scenario")
if (!is.null(lsm1))
  T2 <- T2 |> dplyr::left_join(
    lsm1 |> dplyr::transmute(scenario,
                              top1_NP   = np,
                              top1_ENN_MN  = enn_mn,
                              top1_COHESION = cohesion),
    by = "scenario")
if (!is.null(con_idx_100))
  T2 <- T2 |> dplyr::left_join(con_idx_100, by = "scenario")
tb_save_table(T2, "23_T2_corridor")

## ---- T3: PA coverage --------------------------------------------------------
tb_log_section("T3 PA coverage")

T3 <- data.frame(scenario = scenarios)
T3$scenario_label <- scen_label[T3$scenario]
if (!is.null(pa_sum))
  T3 <- T3 |> dplyr::left_join(
    pa_sum |> dplyr::select(scenario,
                             hs_total_km2,
                             hs_in_pa_km2,
                             hs_pct_in_pa,
                             hs_gap_km2,
                             corr_pct_in_pa),
    by = "scenario")
tb_save_table(T3, "23_T3_pa_coverage")

## Also save per-layer breakdown summary (one row per layer, present only)
if (!is.null(pa_layer)) {
  layer_pres <- pa_layer |> dplyr::filter(scenario == "present") |>
    dplyr::arrange(dplyr::desc(hs_in_layer_km2))
  tb_save_table(layer_pres, "23_T3b_pa_by_layer_present")
}

## ---- T4: Roads --------------------------------------------------------------
tb_log_section("T4 roads")

T4 <- data.frame(scenario = scenarios)
T4$scenario_label <- scen_label[T4$scenario]
if (!is.null(roads))
  T4 <- T4 |> dplyr::left_join(
    roads |> dplyr::select(scenario,
                            corridor_total_km2,
                            corridor_with_road,
                            pct_corr_with_road,
                            road_km_in_corridor,
                            pct_roads_in_corr),
    by = "scenario")
tb_save_table(T4, "23_T4_roads")

## ---- Wide master ------------------------------------------------------------
tb_log_section("Master wide")

master <- T1 |>
  dplyr::left_join(T2 |> dplyr::select(-scenario_label), by = "scenario") |>
  dplyr::left_join(T3 |> dplyr::select(-scenario_label, -hs_total_km2),
                   by = "scenario") |>
  dplyr::left_join(T4 |> dplyr::select(-scenario_label), by = "scenario") |>
  dplyr::mutate(across(where(is.numeric), ~ round(.x, 4)))
tb_save_table(master, "23_master_wide")

## ---- Side bonuses: conflict overlays (if 21 already ran) --------------------
con_clean <- .read_t("21_conflict_clean.csv")
if (!is.null(con_clean)) {
  br_act <- .read_t("21_xtab_bioregion_activity.csv")
  pa_act <- .read_t("21_xtab_pa_activity.csv")
  rd_act <- .read_t("21_road_dist_by_activity.csv")
  if (!is.null(br_act)) tb_save_table(br_act, "23_T5_bioregion_activity")
  if (!is.null(pa_act)) tb_save_table(pa_act, "23_T5_pa_activity")
  if (!is.null(rd_act)) tb_save_table(rd_act, "23_T5_road_dist_activity")
  tb_log("21 conflict overlays copied as T5 set")
}

## ---- Log summary preview ----------------------------------------------------
tb_log_section("Master preview")
print(as.data.frame(master |> dplyr::select(scenario_label, hs_total_km2,
                                              area_corridor_km2, hs_pct_in_pa,
                                              pct_corr_with_road,
                                              PC_100km)))

tb_log_session()
tb_log("23_tables_master DONE")
