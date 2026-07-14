## ============================================================================
## fig6_costdist.R  — regenerate main Fig. 6 (climate-change response 4-panel)
## with panel C (PC ratio) computed on COST-DISTANCE PC (35_costdist_indices.csv)
## instead of Euclidean (17_conefor_indices.csv). CSV-only; no rasters needed.
## Output: figures/fig6_costdist/fig6_scenario_change_costdist.png
## ============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(scales) })

.find <- function() {
  env_dir <- Sys.getenv("TB_SCRIPTS", unset = "")
  if (nzchar(env_dir) && file.exists(file.path(env_dir, "00_paths.R"))) return(env_dir)
  if (file.exists("00_paths.R")) return(getwd())
  stop("00_paths.R not found")
}
setwd(.find()); source("00_paths.R"); source("00_helpers.R")

FIG_SUBDIR <- "fig6_costdist"
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
SCEN_PAL <- setNames(c("#2C2C2C", "#4E79A7", "#F28E2B", "#E15759",
                       "#76B7B2", "#B07AA1", "#9C2A2B"),
                     scen_label[scenarios])

## ---- data (CSV only) -------------------------------------------------------
pa_summary <- read.csv(file.path(TB_OUT_TABLES, "18_pa_overlay_summary.csv"))
con_sum    <- read.csv(file.path(TB_OUT_TABLES, "15_connectivity_summary.csv"))
pa_summary <- merge(pa_summary, con_sum[, c("scenario", "area_corridor_km2")],
                    by = "scenario", all.x = TRUE)
pa_summary$corr_total_km2 <- pa_summary$area_corridor_km2

## *** COST-DISTANCE PC ***
con_idx <- read.csv(file.path(TB_OUT_TABLES, "35_costdist_indices.csv"))

pa_summary$scen_lbl <- factor(scen_label[pa_summary$scenario],
                              levels = scen_label[scenarios])
con_100 <- con_idx[con_idx$d_km == 100, ]
con_100$scen_lbl <- factor(scen_label[con_100$scenario],
                           levels = scen_label[scenarios])

## ---- panels ----------------------------------------------------------------
p2a <- ggplot(pa_summary, aes(scen_lbl, hs_total_km2, fill = scen_lbl)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = SCEN_PAL, guide = "none") +
  labs(title = "A. Suitable habitat (km²)", x = NULL, y = "km²") +
  scale_y_continuous(labels = scales::label_comma()) +
  theme_trbear_bar(base_size = 11) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

p2b <- ggplot(pa_summary, aes(scen_lbl, corr_total_km2, fill = scen_lbl)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = SCEN_PAL, guide = "none") +
  labs(title = "B. Top-5% corridor area (km²)", x = NULL, y = "km²") +
  scale_y_continuous(labels = scales::label_comma()) +
  theme_trbear_bar(base_size = 11) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

PC_pres <- con_100$PC[con_100$scenario == "present"]
con_100$PC_rel <- con_100$PC / PC_pres
p2c <- ggplot(con_100, aes(scen_lbl, PC_rel, group = 1)) +
  geom_hline(yintercept = 1, linetype = 2, color = "#7C8A93", linewidth = 0.4) +
  geom_line(linewidth = 0.9, color = "#1F3A93") +
  geom_point(size = 3, color = "#1F3A93") +
  geom_text(aes(label = ifelse(100 * PC_rel < 1,
                               sprintf("%.1f%%", 100 * PC_rel),
                               sprintf("%.0f%%", 100 * PC_rel))),
            vjust = -0.9, size = 3.2, color = "#1F3A93") +
  scale_y_continuous(labels = scales::label_percent(accuracy = 1),
                     limits = c(0, 1.12), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "C. Conefor PC at d = 100 km (relative to present)",
       x = NULL, y = "PC ratio (present = 100%)") +
  theme_trbear_bar(base_size = 11) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

p2d <- ggplot(pa_summary, aes(scen_lbl, hs_in_pa_km2, fill = scen_lbl)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.2) +
  scale_fill_manual(values = SCEN_PAL, guide = "none") +
  geom_text(aes(label = sprintf("%.1f%%", hs_pct_in_pa)),
            vjust = -0.6, size = 3.0, color = TB_COLOR_AXIS) +
  scale_y_continuous(labels = scales::label_comma(),
                     expand = expansion(mult = c(0, 0.12))) +
  labs(title = "D. Suitable habitat inside PA (km²)", x = NULL, y = "km²") +
  theme_trbear_bar(base_size = 11) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

fig2 <- (p2a + p2b) / (p2c + p2d) +
  patchwork::plot_annotation(
    title = "Climate-change response: habitat, corridor, connectivity and protection",
    theme = theme(plot.title = element_text(face = "bold", size = 15,
                                             color = TB_COLOR_FRAME)))
tb_save_fig(fig2, "fig6_scenario_change_costdist", w = 16, h = 11, subdir = FIG_SUBDIR)

cat("\n[PC ratios, cost-distance, d=100]\n")
print(data.frame(scenario = con_100$scenario,
                 PC = signif(con_100$PC, 4),
                 PC_rel_pct = round(100 * con_100$PC_rel, 1)))
cat("[OK] fig6_scenario_change_costdist.png\n")
