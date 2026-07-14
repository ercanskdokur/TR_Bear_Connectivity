## ============================================================================
## 07_postprocess_eval.R
## Project: TR_Bear_Connectivity
## Purpose: Consolidate ENMTML algorithm evaluation + variable importance into
##   publication-ready tables and figures.
## Inputs (TB_OUT_ENMTML):
##   - Evaluation_Table.txt
##   - Thresholds_Algorithms.txt, Thresholds_Ensemble.txt
##   - Algorithm/<algo>/Response Curves & Variable Importance/VariableImportance.txt
## Outputs:
##   tables/  07_evaluation_summary.csv
##            07_variable_importance_long.csv
##            07_variable_importance_consensus.csv
##            07_thresholds.csv
##   figures/07_eval/  fig07a_metric_dotplot.png
##                     fig07b_varimp_heatmap.png
##                     fig07c_top_predictors_ranked.png
##                     fig07d_ensemble_thresholds.png
## ============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(forcats)
  library(patchwork); library(stringr); library(scales)
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
tb_log_init("07_postprocess_eval")

SP    <- "Ursus_arctos"
ALGOS <- TB_ENM_ALGORITHMS
FIG_SUBDIR <- "07_eval"
TAB_SUBDIR <- NULL  ## tables/ root

## ----------------------------------------------------------------------------
## 1) Evaluation table — consolidate
## ----------------------------------------------------------------------------
tb_log_section("Evaluation table")

eval_raw <- read.table(file.path(TB_OUT_ENMTML, "Evaluation_Table.txt"),
                       header = TRUE, sep = "\t", stringsAsFactors = FALSE)
tb_log(sprintf("Evaluation_Table rows: %d", nrow(eval_raw)))

metric_cols    <- c("AUC","TSS","Kappa","Jaccard","Sorensen","Boyce","Fpb","OR")
metric_sd_cols <- paste0(metric_cols, "_SD")

eval_long <- eval_raw |>
  select(Algorithm, all_of(c(metric_cols, metric_sd_cols))) |>
  pivot_longer(cols = -Algorithm,
               names_to = "stat", values_to = "value") |>
  mutate(
    metric = stringr::str_remove(stat, "_SD$"),
    kind   = ifelse(grepl("_SD$", stat), "sd", "mean")
  ) |>
  select(-stat) |>
  pivot_wider(names_from = kind, values_from = value)

## NOTE: ENMTML reports Kappa == TSS for every algorithm in this run (BOOT/MAX_TSS).
## Drop Kappa from the summary to avoid duplicating a column.
eval_summary <- eval_long |>
  filter(metric != "Kappa") |>
  pivot_wider(names_from = metric, values_from = c(mean, sd), names_glue = "{metric}_{.value}") |>
  arrange(desc(TSS_mean))

tb_save_table(eval_summary, "07_evaluation_summary", subdir = TAB_SUBDIR)

## ----------------------------------------------------------------------------
## 2) Variable importance — load 8 files, consolidate
## ----------------------------------------------------------------------------
tb_log_section("Variable importance")

## NOTE: ENMTML writes the importance column with algorithm-specific name
##  (BRT: "Importance", RDF: "IncNodePurity", others: "Overall"). Pull the 4th
##  column regardless of header to keep things consistent.
vi_list <- lapply(ALGOS, function(a) {
  fn <- file.path(TB_OUT_ENMTML, "Algorithm", a,
                  "Response Curves & Variable Importance",
                  "VariableImportance.txt")
  if (!file.exists(fn)) { tb_log(sprintf("MISSING: %s", fn), "WARN"); return(NULL) }
  d <- read.table(fn, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  data.frame(Algorithm  = a,
             Variables  = as.character(d$Variables),
             Importance = as.numeric(d[[4]]),
             stringsAsFactors = FALSE)
})
vi_long <- do.call(rbind, vi_list)

## GLM uses polynomial expansion (I(X^2)); other algorithms don't. To keep the
## per-algorithm normalization fair and consensus interpretable, fold the
## squared terms back into their base variable (sum I(X^2) + X within GLM)
## before any further analysis.
vi_long <- vi_long |>
  mutate(BaseVar = sub("^I\\((.*)\\^2\\)$", "\\1", Variables)) |>
  group_by(Algorithm, BaseVar) |>
  summarise(Importance = sum(Importance, na.rm = TRUE), .groups = "drop") |>
  rename(Variables = BaseVar)

## Normalize per-algorithm so each algo sums to 1 (different algos report
## importance in different intrinsic units — variable rank is comparable,
## absolute magnitude is not).
vi_long <- vi_long |>
  group_by(Algorithm) |>
  mutate(Importance_norm = Importance / sum(Importance, na.rm = TRUE)) |>
  ungroup()

tb_save_table(vi_long, "07_variable_importance_long", subdir = TAB_SUBDIR)

## Treat "variable not used by this algorithm" (e.g., GAM variable selection)
## as importance = 0 across all known predictors, not NA — otherwise the
## consensus mean is biased upward by ignoring zero-importance algorithms.
all_vars_seen <- sort(unique(vi_long$Variables))
vi_full <- tidyr::expand_grid(Algorithm = ALGOS,
                              Variables = all_vars_seen) |>
  left_join(vi_long, by = c("Algorithm", "Variables")) |>
  mutate(
    Importance      = ifelse(is.na(Importance),      0, Importance),
    Importance_norm = ifelse(is.na(Importance_norm), 0, Importance_norm)
  )

vi_consensus <- vi_full |>
  group_by(Variables) |>
  summarise(
    mean_imp_norm = mean(Importance_norm, na.rm = TRUE),
    sd_imp_norm   = sd(Importance_norm,   na.rm = TRUE),
    n_algos_used  = sum(Importance_norm > 0),
    .groups = "drop"
  ) |>
  arrange(desc(mean_imp_norm))

tb_save_table(vi_consensus, "07_variable_importance_consensus", subdir = TAB_SUBDIR)

## ----------------------------------------------------------------------------
## 3) Thresholds table
## ----------------------------------------------------------------------------
tb_log_section("Thresholds")

## Thresholds_Algorithms.txt holds 6 rows per algorithm (MAX_TSS, LPT,
## MAX_KAPPA, SENSITIVITY, JACCARD, SORENSEN). Only MAX_TSS is used in this
## pipeline — keep that row and drop the rest.
thr_algo <- read.table(file.path(TB_OUT_ENMTML, "Thresholds_Algorithms.txt"),
                       header = TRUE, sep = "\t", stringsAsFactors = FALSE)
thr_algo <- thr_algo[thr_algo$THR == "MAX_TSS", ]
thr_ens  <- read.table(file.path(TB_OUT_ENMTML, "Thresholds_Ensemble.txt"),
                       header = TRUE, sep = "\t", stringsAsFactors = FALSE)
thr_ens$Algorithm <- paste0("Ensemble_", thr_ens$Ensemble)
thr_ens$Ensemble  <- NULL

common <- intersect(names(thr_algo), names(thr_ens))
thr_all <- rbind(thr_algo[, common], thr_ens[, common])
tb_save_table(thr_all, "07_thresholds", subdir = TAB_SUBDIR)

## ============================================================================
## FIGURES
## ============================================================================

## ---- fig07a: metric dotplot with error bars (4 panel: AUC/TSS/Boyce/Kappa) ---
tb_log_section("fig07a metric dotplot")

## Replace Kappa (== TSS in this run) with Sorensen — independent information.
m_keep <- c("TSS","AUC","Boyce","Sorensen")
dot_df <- eval_long |>
  filter(metric %in% m_keep) |>
  mutate(
    metric    = factor(metric, levels = m_keep),
    Algorithm = factor(Algorithm, levels = eval_summary$Algorithm)
  )

p07a <- ggplot(dot_df, aes(y = Algorithm, x = mean, color = Algorithm)) +
  geom_vline(data = data.frame(metric = m_keep,
                               ref = c(0.7, 0.8, 0.5, 0.7)),
             aes(xintercept = ref), linetype = "dashed",
             color = "gray60", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = mean - sd, xmax = mean + sd),
                 height = 0.25, linewidth = 0.6) +
  geom_point(size = 3.5) +
  facet_wrap(~ metric, scales = "free_x", nrow = 1) +
  scale_color_manual(values = setNames(
    TB_PAL_OKABE_ITO[seq_along(eval_summary$Algorithm)],
    eval_summary$Algorithm)) +
  guides(color = "none") +
  labs(title    = "Algorithm performance — ENMTML BOOT (10 reps)",
       subtitle = "Mean ± SD across bootstrap replicates. Dashed: conventional acceptability thresholds.",
       x = NULL, y = NULL) +
  theme_trbear_bar()

tb_save_fig(p07a, "fig07a_metric_dotplot", w = 13, h = 6, subdir = FIG_SUBDIR)

## ---- fig07b: variable importance heatmap + consensus barplot ----------------
## ENMTML 1.0.0's permutation-based variable importance returns all-zero values
## for GAM (spline-based formulation is incompatible with the default routine);
## we drop the GAM column from the heatmap and from the consensus to avoid
## misrepresenting a missing measurement as zero importance.  This matches the
## consensus already saved in 07_variable_importance_consensus.csv (n_algos = 7).
tb_log_section("fig07b varimp heatmap (GAM excluded)")

ALGOS_KEPT <- setdiff(ALGOS, "GAM")
var_order  <- vi_consensus$Variables  ## descending mean importance
vi_full_p  <- vi_full |>
  filter(Algorithm %in% ALGOS_KEPT) |>
  mutate(
    Variables = factor(Variables, levels = rev(var_order)),
    Algorithm = factor(Algorithm, levels = ALGOS_KEPT)
  )

p07b_heat <- ggplot(vi_full_p,
                    aes(x = Algorithm, y = Variables, fill = Importance_norm)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_viridis_c(option = "mako", direction = -1,
                       name = "Normalized\nimportance",
                       labels = scales::percent_format(accuracy = 1)) +
  labs(title = "Per algorithm", x = NULL, y = NULL) +
  theme_trbear_bar() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        legend.position = "bottom",
        legend.key.width = grid::unit(1.2, "cm"))

p07b_bar <- ggplot(vi_consensus |>
                     mutate(Variables = factor(Variables, levels = rev(var_order))),
                   aes(x = mean_imp_norm, y = Variables)) +
  geom_col(fill = "#0072B2", alpha = 0.85) +
  geom_errorbarh(aes(xmin = pmax(mean_imp_norm - sd_imp_norm, 0),
                     xmax = mean_imp_norm + sd_imp_norm),
                 height = 0.25, color = "#263238") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(title = "Consensus (mean ± SD)",
       x = "Importance", y = NULL) +
  theme_trbear_bar() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

p07b <- (p07b_heat + p07b_bar) +
  patchwork::plot_layout(widths = c(2.6, 1)) +
  patchwork::plot_annotation(
    title    = "Variable importance across algorithms (GAM excluded)",
    subtitle = "Per-algorithm normalized (column sums to 100%); polynomial terms folded into base variable.",
    caption  = "GAM omitted because ENMTML's permutation-based importance routine returns zeros for spline-based models.",
    theme    = theme(plot.title    = element_text(face = "bold", size = 14,
                                                    color = TB_COLOR_FRAME),
                     plot.subtitle = element_text(size = 11,
                                                    color = TB_COLOR_FRAME),
                     plot.caption  = element_text(size = 9,
                                                    color = "grey40",
                                                    hjust = 0))
  )
tb_save_fig(p07b, "fig07b_varimp_heatmap", w = 14, h = 11, subdir = FIG_SUBDIR)

## ---- fig07c: top predictors consensus ranked --------------------------------
tb_log_section("fig07c top predictors ranked")

vi_top <- vi_consensus |> slice_head(n = 15)

p07c <- ggplot(vi_top,
               aes(x = mean_imp_norm,
                   y = fct_reorder(Variables, mean_imp_norm))) +
  geom_col(aes(fill = mean_imp_norm), show.legend = FALSE) +
  geom_errorbarh(aes(xmin = pmax(mean_imp_norm - sd_imp_norm, 0),
                     xmax = mean_imp_norm + sd_imp_norm),
                 height = 0.25, color = "#263238", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f%%", mean_imp_norm * 100)),
            hjust = -0.2, size = 3.6, color = "#263238") +
  scale_fill_viridis_c(option = "mako", direction = -1) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.18))) +
  labs(title = "Top 15 predictors (consensus across 8 algorithms)",
       subtitle = "Bars: mean normalized importance; error bars: ±SD across algorithms",
       x = "Mean importance", y = NULL) +
  theme_trbear_bar()

tb_save_fig(p07c, "fig07c_top_predictors_ranked", w = 11, h = 9, subdir = FIG_SUBDIR)

## ---- fig07d: ensemble vs per-algo thresholds --------------------------------
tb_log_section("fig07d thresholds")

thr_p <- thr_all |>
  mutate(Type = ifelse(grepl("^Ensemble_", Algorithm), "Ensemble", "Per-algorithm"),
         Algorithm = factor(Algorithm, levels = c(ALGOS, paste0("Ensemble_", c("MEAN","W_MEAN")))))

p07d <- ggplot(thr_p, aes(x = THR_VALUE,
                          y = fct_reorder(Algorithm, THR_VALUE),
                          fill = Type)) +
  geom_col(width = 0.65, alpha = 0.9) +
  geom_text(aes(label = sprintf("%.3f", THR_VALUE)),
            hjust = -0.15, size = 3.6, color = "#263238") +
  scale_fill_manual(values = c("Per-algorithm" = "#56B4E9",
                               "Ensemble"      = "#D55E00")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "MAX_TSS thresholds",
       subtitle = "Per-algorithm (ENMTML internal) vs. Ensemble (computed via bg-sampling)",
       x = "Threshold (HS)", y = NULL, fill = NULL) +
  theme_trbear_bar() +
  theme(legend.position = "top")

tb_save_fig(p07d, "fig07d_ensemble_thresholds", w = 11, h = 7, subdir = FIG_SUBDIR)

tb_log_session()
tb_log("07_postprocess_eval DONE")
