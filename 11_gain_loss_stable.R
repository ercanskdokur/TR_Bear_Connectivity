## ============================================================================
## 11_gain_loss_stable.R
## Project: TR_Bear_Connectivity
## Purpose: Categorical range change vs present, for each of the 6 GCM-averaged
##   future scenarios. 4 classes:
##     0 = Stable unsuitable (both unsuit)
##     1 = Loss              (present suit, future unsuit)
##     2 = Gain              (present unsuit, future suit)
##     3 = Stable suitable   (both suit)
##   Outputs: per-scenario categorical raster + summary table + stacked bar.

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
tb_log_init("11_gain_loss_stable")

SP <- "Ursus_arctos"
FIG_SUBDIR <- "11_gainloss"

CLASS_LEVELS <- c("Unsuitable","Loss","Gain","Suitable")
CLASS_VALUES <- c(0, 1, 2, 3)
## More vivid (colour-blind-safe) palette:
##   Unsuitable = warm grey, Loss = vivid orange-red,
##   Gain = teal-green, Suitable = saturated blue.
CLASS_COLOR  <- c("Unsuitable" = "#9E9E9E",
                  "Loss"       = "#E15759",
                  "Gain"       = "#59A14F",
                  "Suitable"   = "#4E79A7")

## ----------------------------------------------------------------------------
## 1) Inputs
## ----------------------------------------------------------------------------
tb_log_section("Load present + 6 future binaries")

present_bin <- terra::rast(file.path(TB_OUT_HS_BINARY, "present_wmean.tif"))
names(present_bin) <- "present"

scen_avg <- expand.grid(period = TB_PERIODS, ssp = TB_SSPS,
                        stringsAsFactors = FALSE)
scen_avg <- scen_avg[order(scen_avg$period, scen_avg$ssp), ]
rownames(scen_avg) <- NULL
scen_avg$scen_key <- sprintf("%s_%s", scen_avg$period, scen_avg$ssp)

## ----------------------------------------------------------------------------
## 2) Build gain/loss raster per scenario + centroid stats
## ----------------------------------------------------------------------------
tb_log_section("Build gain/loss rasters")

cell_area_km2 <- prod(terra::res(present_bin)) / 1e6

gl_rasters <- vector("list", nrow(scen_avg))
gl_summary <- vector("list", nrow(scen_avg))

for (i in seq_len(nrow(scen_avg))) {
  k <- scen_avg$scen_key[i]
  fn <- file.path(TB_OUT_HS_BINARY, sprintf("future_%s.tif", k))
  if (!file.exists(fn)) { tb_log(sprintf("missing %s", fn), "WARN"); next }
  fut <- terra::rast(fn); names(fut) <- "future"

  ## Align in case of trivial differences (should be identical grid)
  if (!terra::compareGeom(present_bin, fut, stopOnError = FALSE)) {
    tb_log("geom mismatch — resampling future to present grid", "WARN")
    fut <- terra::resample(fut, present_bin, method = "near")
  }

  ## Combined class:
  ##   0=SU (both 0), 1=Loss (p=1,f=0), 2=Gain (p=0,f=1), 3=SS (both 1)
  cls <- present_bin * 2 - fut  ## not quite — use logical formula
  cls <- terra::ifel(present_bin == 0 & fut == 0, 0,
            terra::ifel(present_bin == 1 & fut == 0, 1,
              terra::ifel(present_bin == 0 & fut == 1, 2,
                terra::ifel(present_bin == 1 & fut == 1, 3, NA))))
  names(cls) <- "change"

  ## Save categorical raster
  fn_out <- file.path(TB_OUT_GAINLOSS, sprintf("%s.tif", k))
  terra::writeRaster(cls, fn_out, overwrite = TRUE, datatype = "INT1U",
                     gdal = c("COMPRESS=DEFLATE","TILED=YES"))

  gl_rasters[[i]] <- cls

  ## Tabulate
  v <- terra::values(cls, mat = FALSE); v <- v[!is.na(v)]
  tab <- table(factor(v, levels = CLASS_VALUES))
  tot <- sum(tab)
  gl_summary[[i]] <- data.frame(
    scenario = k,
    period   = scen_avg$period[i],
    ssp      = scen_avg$ssp[i],
    area_SU_km2   = as.numeric(tab["0"]) * cell_area_km2,
    area_Loss_km2 = as.numeric(tab["1"]) * cell_area_km2,
    area_Gain_km2 = as.numeric(tab["2"]) * cell_area_km2,
    area_SS_km2   = as.numeric(tab["3"]) * cell_area_km2,
    pct_SU   = 100 * as.numeric(tab["0"]) / tot,
    pct_Loss = 100 * as.numeric(tab["1"]) / tot,
    pct_Gain = 100 * as.numeric(tab["2"]) / tot,
    pct_SS   = 100 * as.numeric(tab["3"]) / tot
  )
}

gl_summary_df <- do.call(rbind, gl_summary)
tb_save_table(gl_summary_df, "11_gainloss_summary")

## ----------------------------------------------------------------------------
## FIGURES
## ----------------------------------------------------------------------------
world_sf <- tryCatch(
  rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") |>
    sf::st_transform(TB_CRS_PROJ),
  error = function(e) NULL)
tr_mask_shp <- file.path(TB_DATA_ROOT, "TR_mask", "TR_mask.shp")
tr_mask_sf <- if (file.exists(tr_mask_shp))
  sf::st_read(tr_mask_shp, quiet = TRUE) |> sf::st_transform(TB_CRS_PROJ) else NULL

template <- gl_rasters[[which(!sapply(gl_rasters, is.null))[1]]]
e <- terra::ext(template); ext_pad <- 30000
xlim_p <- c(e$xmin - ext_pad, e$xmax + ext_pad)
ylim_p <- c(e$ymin - ext_pad, e$ymax + ext_pad)

## ---- fig11a: 6-panel categorical maps ---------------------------------------
tb_log_section("fig11a 6-panel categorical")

.tb_panel_cls <- function(r, title, show_legend = FALSE) {
  rfac <- r
  levels(rfac) <- data.frame(value = CLASS_VALUES, class = CLASS_LEVELS)
  p <- ggplot()
  if (!is.null(world_sf)) p <- p + geom_sf(data = world_sf,
                                            fill = TB_FILL_LAND,
                                            color = TB_COLOR_LAND,
                                            linewidth = 0.25)
  p <- p +
    tidyterra::geom_spatraster(data = rfac, na.rm = TRUE) +
    scale_fill_manual(values = CLASS_COLOR,
                      breaks = CLASS_LEVELS,
                      na.translate = FALSE, drop = FALSE,
                      name = "Class")
  if (!is.null(tr_mask_sf)) p <- p +
    geom_sf(data = tr_mask_sf, fill = NA, color = TB_COLOR_FRAME,
            linewidth = 0.35)
  p <- p +
    coord_sf(xlim = xlim_p, ylim = ylim_p, datum = sf::st_crs(4326),
             expand = FALSE) +
    labs(title = title) +
    theme_trbear(base_size = 11) +
    theme(plot.title = element_text(size = 12, hjust = 0.5,
                                    margin = ggplot2::margin(b = 3)),
          plot.margin = ggplot2::margin(3, 3, 3, 3))
  if (!show_legend) p <- p + theme(legend.position = "none")
  p
}

plots_a <- vector("list", nrow(scen_avg))
for (i in seq_len(nrow(scen_avg))) {
  r <- gl_rasters[[i]]
  if (is.null(r)) { plots_a[[i]] <- patchwork::plot_spacer(); next }
  ttl <- sprintf("%s — %s",
                 TB_PERIOD_LABELS[scen_avg$period[i]],
                 TB_SSP_LABELS[scen_avg$ssp[i]])
  plots_a[[i]] <- .tb_panel_cls(r, ttl, show_legend = (i == nrow(scen_avg)))
}
ord_a <- with(scen_avg, order(period, ssp))
plots_a <- plots_a[ord_a]

p11a <- patchwork::wrap_plots(plots_a, ncol = 3, nrow = 2, guides = "collect") +
  patchwork::plot_annotation(
    title    = "Habitat changes classes – GCMs averaged futures",
    subtitle = "Loss = suit→unsuit, Gain = unsuit→suit; thresholded at present MAX_TSS.",
    theme = theme(plot.title    = element_text(face = "bold", size = 18,
                                               color = TB_COLOR_FRAME),
                  plot.subtitle = element_text(size = 12, color = TB_COLOR_AXIS)))
tb_save_fig(p11a, "fig11a_gainloss_6panel", w = 22, h = 13, subdir = FIG_SUBDIR)

## ---- fig11b: stacked bar (% area per class × scenario) ----------------------
tb_log_section("fig11b stacked bar")

bar_df <- gl_summary_df |>
  select(scenario, period, ssp, pct_SU, pct_Loss, pct_Gain, pct_SS) |>
  pivot_longer(cols = starts_with("pct_"),
               names_to = "class", values_to = "pct") |>
  mutate(
    class = recode(class,
                   pct_SU   = "Unsuitable",
                   pct_Loss = "Loss",
                   pct_Gain = "Gain",
                   pct_SS   = "Suitable"),
    period_lbl = TB_PERIOD_LABELS[period],
    ssp_lbl    = TB_SSP_LABELS[ssp],
    scen_lbl   = sprintf("%s\n%s",
                         gsub(" \\(.*\\)", "", period_lbl),
                         gsub(" \\(.*\\)", "", ssp_lbl))
  )

## Order class levels so that the LARGEST mean proportion is at the BOTTOM of
## the stack (default ggplot stacking is first level on bottom).  This puts
## Unsuitable (~80%) at the bottom and Gain (<0.2%) at the top.
class_order <- bar_df |>
  group_by(class) |>
  summarise(mean_pct = mean(pct, na.rm = TRUE), .groups = "drop") |>
  arrange(desc(mean_pct)) |>
  pull(class)
bar_df$class <- factor(bar_df$class, levels = class_order)

bar_df$scen_lbl <- factor(bar_df$scen_lbl,
                          levels = unique(bar_df$scen_lbl[order(bar_df$period, bar_df$ssp)]))

## Sort by class factor order, then cumsum gives bottom-up y positions that
## match position_stack(reverse = TRUE) stacking (Unsuitable bottom, Gain top).
## We compute y_pos manually because position_stack(vjust, reverse=TRUE) on
## geom_text put labels at the wrong segment midpoints in this layout.
bar_df <- bar_df |>
  dplyr::arrange(scen_lbl, class) |>
  dplyr::group_by(scen_lbl) |>
  dplyr::mutate(y_pos = cumsum(pct) - pct/2) |>
  dplyr::ungroup()

p11b <- ggplot(bar_df, aes(x = scen_lbl, y = pct, fill = class)) +
  geom_col(width = 0.65, color = "white", linewidth = 0.3,
           position = position_stack(reverse = TRUE)) +
  geom_text(data = subset(bar_df, class %in% c("Loss","Gain","Suitable") & pct >= 0.5),
            aes(y = y_pos, label = sprintf("%.1f%%", pct)),
            color = "white", size = 3.3, fontface = "bold") +
  scale_fill_manual(values = CLASS_COLOR, name = "Class",
                    breaks = class_order) +
  scale_y_continuous(labels = scales::label_percent(scale = 1),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Range change composition",
       x = NULL, y = "Area share") +
  theme_trbear_bar() +
  theme(legend.position = "top",
        axis.text.x = element_text(size = 10))

tb_save_fig(p11b, "fig11b_gainloss_stacked_bar", w = 13, h = 7, subdir = FIG_SUBDIR)

tb_log_session()
tb_log("11_gain_loss_stable DONE")
