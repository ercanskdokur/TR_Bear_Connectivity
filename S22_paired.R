## ============================================================================
## S22_paired.R  — Fig. S22 (two-panel, paired dPC vs dIIC)
## Project: TR_Bear_Connectivity (connectivity manuscript)
## Purpose: Show area-weighted (dPC) versus topological (dIIC) patch importance
##   side by side, making the PC-IIC dichotomy (H5) visible at a glance. Cores
##   that rank highly under dIIC but not dPC (C10, C40, C25) are highlighted as
##   stepping stones; cores dominant under both (C01, C02) are labelled in both
##   panels. Present scenario, d = 100 km.
## Output: figures/S22_paired/figS22_dpc_diic_paired.png
## ============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("ggrepel", quietly = TRUE))
    install.packages("ggrepel", repos = "https://cloud.r-project.org")
  library(terra); library(sf); library(ggplot2); library(dplyr)
  library(tidyterra); library(patchwork); library(rnaturalearth); library(ggrepel)
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
tb_log_init("S22_paired")

FIG_SUBDIR <- "S22_paired"
FOCAL_D    <- 100
STEP_IDS   <- c("C40", "C39", "C37")     # stepping stones to highlight (cost-distance)
DOM_IDS    <- c("C01", "C02")            # dominant cores

## ---- data ------------------------------------------------------------------
dpc <- read.csv(file.path(TB_OUT_TABLES, "35_costdist_dpc_dii.csv")) |>   # cost-distance
  dplyr::filter(scenario == "present", d_km == FOCAL_D)
cw  <- read.csv(file.path(TB_OUT_TABLES, "34_core_crosswalk.csv"))[, c("patch_id","core_id")]
dpc <- dplyr::left_join(dpc, cw, by = "patch_id")

bin <- terra::rast(file.path(TB_OUT_HS_BINARY, "present_wmean.tif"))
tr_sf <- sf::st_read(file.path(TB_DATA_ROOT, "TR_mask", "TR_mask.shp"), quiet = TRUE) |>
  sf::st_transform(TB_CRS_PROJ)
bin <- terra::mask(bin, terra::vect(tr_sf))
bf  <- terra::as.factor(terra::as.int(bin))
levels(bf) <- data.frame(id = c(0, 1), class = c("Unsuitable", "Suitable"))
names(bf) <- "class"

world_sf <- tryCatch(rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") |>
                       sf::st_transform(TB_CRS_PROJ), error = function(e) NULL)
e <- terra::ext(bin); pad <- 30000
xl <- c(e$xmin - pad, e$xmax + pad); yl <- c(e$ymin - pad, e$ymax + pad)

## ---- panel builder ---------------------------------------------------------
mk_panel <- function(index, accent, title, label_ids) {
  top <- dpc |> dplyr::arrange(dplyr::desc(.data[[index]])) |> dplyr::slice_head(n = 10)
  top$grp <- ifelse(top$core_id %in% STEP_IDS, "dIIC > dPC", "Top-10 core")
  top_sf  <- sf::st_as_sf(top, coords = c("x","y"), crs = TB_CRS_PROJ)
  lab     <- top[top$core_id %in% label_ids, ]
  lab_sf  <- sf::st_as_sf(lab, coords = c("x","y"), crs = TB_CRS_PROJ)

  p <- ggplot()
  if (!is.null(world_sf)) p <- p +
    geom_sf(data = world_sf, fill = "#E8E8E8", color = "#7C8A93", linewidth = 0.3)
  p <- p +
    tidyterra::geom_spatraster(data = bf, na.rm = TRUE) +
    scale_fill_manual(values = TB_PAL_BINARY, na.translate = FALSE, name = "Habitat") +
    geom_sf(data = tr_sf, fill = NA, color = TB_COLOR_FRAME, linewidth = 0.4) +
    geom_sf(data = top_sf, aes(size = .data[[index]], color = grp),
            shape = 16, alpha = 0.78) +
    geom_sf(data = top_sf, aes(size = .data[[index]]),
            shape = 21, color = "black", fill = NA, stroke = 0.4) +
    scale_color_manual(values = c("Top-10 core" = accent,
                                  "dIIC > dPC" = "#D4A017"),
                       name = NULL) +
    scale_size_continuous(name = sprintf("%s (%%)", index), range = c(3, 15)) +
    ggrepel::geom_text_repel(data = lab_sf, aes(label = core_id, geometry = geometry),
                             stat = "sf_coordinates", size = 3.4, fontface = "bold",
                             color = "black", box.padding = 0.5, min.segment.length = 0,
                             seed = 1) +
    coord_sf(xlim = xl, ylim = yl, datum = sf::st_crs(4326), expand = FALSE) +
    annotation_scale(location = "bl", width_hint = 0.18, style = "ticks",
                     text_col = TB_COLOR_FRAME, line_col = TB_COLOR_FRAME) +
    labs(title = title) +
    theme_trbear(base_size = 11)
  p
}

pa <- mk_panel("dPC",  "#0072B2",
               "Area-weighted importance (top-10 dPC)", DOM_IDS)
pb <- mk_panel("dIIC", "#9E2A2B",
               "Topological importance (top-10 dIIC)", c(DOM_IDS, STEP_IDS))

p <- (pa / pb) +
  patchwork::plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 14))
tb_save_fig(p, "figS22_dpc_diic_paired", w = 10, h = 12.5, subdir = FIG_SUBDIR)
tb_log_session(); tb_log("S22_paired DONE")
