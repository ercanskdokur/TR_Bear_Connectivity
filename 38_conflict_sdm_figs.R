## ============================================================================
## 38_conflict_sdm_figs.R
## Project: TR_Bear_Connectivity
## Purpose: Document the human-bear conflict (HBC) ensemble model with the SAME
##   diagnostics already shown for the bear ENM (manuscript Fig. S3-S6 / Table
##   S2-S3), so the conflict surface that feeds road-crossing prioritisation and
##   the compound-risk layer is fully documented.
##
##   Reads the parallel ENMTML conflict run (enmtml_conflict_result):
##     - Evaluation_Table.txt                  -> per-algorithm AUC/TSS (S3-like)
##     - Algorithm/*/.../VariableImportance.txt -> consensus importance (S4-like)
##     - Thresholds_Ensemble.txt               -> MAX_TSS threshold (S5-like)
##     - Ensemble/W_MEAN/Ursus_arctos_conflict.tif -> suitability map (S6-like)
##
## Outputs:
##   tables/38_conflict_eval.csv          per-algorithm cross-validated metrics
##   tables/38_conflict_importance.csv    consensus predictor importance
##   figures/38_conflict_sdm/
##     fig38a_conflict_performance.png     per-algorithm TSS & AUC
##     fig38b_conflict_importance.png      consensus variable importance
##     fig38c_conflict_suitability.png     continuous HBC risk surface (W_MEAN)
##     fig38d_conflict_binary.png          binary HBC risk (MAX_TSS)
## ============================================================================

suppressPackageStartupMessages({
  library(terra); library(sf); library(ggplot2); library(dplyr); library(tidyr)
  library(tidyterra); library(rnaturalearth)
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
tb_log_init("38_conflict_sdm_figs")

FIG_SUBDIR <- "38_conflict_sdm"
CR <- TB_OUT_ENMTML_CONFLICT
ALGS <- c("BIO","GLM","GAM","SVM","RDF","BRT","MXD","MAH")

## ----------------------------------------------------------------------------
## 1) Per-algorithm performance (S3-like)
## ----------------------------------------------------------------------------
tb_log_section("Performance table")
ev <- read.delim(file.path(CR, "Evaluation_Table.txt"), check.names = FALSE)
ev_alg <- ev |> dplyr::filter(Algorithm %in% ALGS) |>
  dplyr::transmute(Algorithm, AUC, TSS, Boyce, Jaccard, Sorensen, OR,
                   AUC_SD, TSS_SD) |>
  dplyr::arrange(dplyr::desc(TSS))
tb_save_table(ev_alg, "38_conflict_eval")

pa_df <- ev_alg |> tidyr::pivot_longer(c(AUC, TSS), names_to = "metric", values_to = "val")
pa_df$Algorithm <- factor(pa_df$Algorithm, levels = ev_alg$Algorithm)
p38a <- ggplot(pa_df, aes(Algorithm, val, fill = metric)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_hline(yintercept = c(0.5, 0.7), linetype = 3, color = "gray55") +
  scale_fill_manual(values = c("AUC" = "#0072B2", "TSS" = "#D55E00"), name = NULL) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.03))) +
  labs(title = "Human-bear conflict ensemble: per-algorithm discrimination",
       subtitle = "Bootstrap cross-validation (10 replicates, 70/30 split)",
       x = NULL, y = "Score") +
  theme_trbear_bar(base_size = 12)
tb_save_fig(p38a, "fig38a_conflict_performance", w = 11, h = 6, subdir = FIG_SUBDIR)

## ----------------------------------------------------------------------------
## 2) Consensus variable importance (S4-like)
## ----------------------------------------------------------------------------
tb_log_section("Variable importance")
imp_list <- lapply(setdiff(ALGS, "MAH"), function(al) {
  f <- file.path(CR, "Algorithm", al, "Response Curves & Variable Importance",
                 "VariableImportance.txt")
  if (!file.exists(f)) return(NULL)
  d <- read.delim(f, check.names = FALSE)
  nc <- ncol(d)   # last two cols are always Variables, Overall (robust to row-name shift)
  data.frame(Algorithm = al,
             Variables = as.character(d[[nc - 1]]),
             imp = suppressWarnings(as.numeric(d[[nc]])))
})
imp <- do.call(rbind, imp_list)
imp_sum <- imp |> dplyr::group_by(Variables) |>
  dplyr::summarise(mean_imp = mean(imp, na.rm = TRUE),
                   sd_imp = sd(imp, na.rm = TRUE),
                   n_algos = dplyr::n(), .groups = "drop") |>
  dplyr::arrange(dplyr::desc(mean_imp))
tb_save_table(imp_sum, "38_conflict_importance")

imp_sum$Variables <- factor(imp_sum$Variables, levels = rev(imp_sum$Variables))
p38b <- ggplot(imp_sum, aes(mean_imp, Variables)) +
  geom_col(fill = "#9E2A2B", alpha = 0.85) +
  geom_errorbarh(aes(xmin = pmax(0, mean_imp - sd_imp), xmax = mean_imp + sd_imp),
                 height = 0.3, color = "gray35") +
  labs(title = "Human-bear conflict model: consensus predictor importance",
       subtitle = "Mean permutation importance across algorithms (± SD)",
       x = "Normalised importance", y = NULL) +
  theme_trbear_bar(base_size = 12)
tb_save_fig(p38b, "fig38b_conflict_importance", w = 10, h = 8, subdir = FIG_SUBDIR)

## ----------------------------------------------------------------------------
## 3) Suitability maps (S6-like + binary)
## ----------------------------------------------------------------------------
tb_log_section("Suitability maps")
wm <- terra::rast(file.path(CR, "Ensemble", "W_MEAN", "Ursus_arctos_conflict.tif"))
if (terra::nlyr(wm) > 1) wm <- wm[[1]]
if (max(terra::values(wm, mat = FALSE), na.rm = TRUE) > 1.5) wm <- wm / 1000  # ENMTML 0-1000 scale

world_sf <- tryCatch(rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") |>
  sf::st_transform(TB_CRS_PROJ), error = function(e) NULL)
tr_sf <- sf::st_read(file.path(TB_DATA_ROOT, "TR_mask", "TR_mask.shp"), quiet = TRUE) |>
  sf::st_transform(TB_CRS_PROJ)
wm <- terra::mask(wm, terra::vect(tr_sf))
e <- terra::ext(wm); pad <- 30000
xl <- c(e$xmin - pad, e$xmax + pad); yl <- c(e$ymin - pad, e$ymax + pad)

base_map <- function() {
  p <- ggplot()
  if (!is.null(world_sf)) p <- p + geom_sf(data = world_sf, fill = "#E8E8E8",
    color = "#7C8A93", linewidth = 0.4)
  p
}
p38c <- base_map() +
  tidyterra::geom_spatraster(data = wm, na.rm = TRUE) +
  scale_fill_viridis_c(option = "inferno", direction = -1, na.value = "transparent",
                       limits = c(0, 1), name = "HBC risk\n(W_MEAN)") +
  geom_sf(data = tr_sf, fill = NA, color = TB_COLOR_FRAME, linewidth = 0.5) +
  coord_sf(xlim = xl, ylim = yl, datum = sf::st_crs(4326), expand = FALSE) +
  tb_map_decorations() +
  labs(title = "Human-bear conflict risk across Türkiye (present, W_MEAN ensemble)") +
  theme_trbear()
tb_save_fig(p38c, "fig38c_conflict_suitability", w = 14, h = 9, subdir = FIG_SUBDIR)

## binary at MAX_TSS (ensemble)
thr <- tryCatch({
  td <- read.delim(file.path(CR, "Thresholds_Ensemble.txt"), check.names = FALSE)
  wrow <- td[grepl("W_MEAN|WMEA", td[[2]], ignore.case = TRUE), ]
  as.numeric(wrow$THR[1] %||% wrow[["MAX_TSS"]][1] %||% NA)
}, error = function(e) NA)
if (is.na(thr)) thr <- as.numeric(stats::quantile(terra::values(wm, mat = FALSE),
                                                  0.7, na.rm = TRUE))
tb_log(sprintf("conflict MAX_TSS threshold = %.3f", thr))
binr <- terra::as.factor(terra::as.int(wm >= thr))
levels(binr) <- data.frame(id = c(0, 1), class = c("Low risk", "High risk"))
names(binr) <- "class"
p38d <- base_map() +
  tidyterra::geom_spatraster(data = binr, na.rm = TRUE) +
  scale_fill_manual(values = c("Low risk" = "#F2F2F2", "High risk" = "#9E2A2B"),
                    na.translate = FALSE, name = NULL) +
  geom_sf(data = tr_sf, fill = NA, color = TB_COLOR_FRAME, linewidth = 0.5) +
  coord_sf(xlim = xl, ylim = yl, datum = sf::st_crs(4326), expand = FALSE) +
  tb_map_decorations() +
  labs(title = "Binary human-bear conflict risk (MAX_TSS threshold)") +
  theme_trbear()
tb_save_fig(p38d, "fig38d_conflict_binary", w = 14, h = 9, subdir = FIG_SUBDIR)

tb_log_session(); tb_log("38_conflict_sdm_figs DONE")
