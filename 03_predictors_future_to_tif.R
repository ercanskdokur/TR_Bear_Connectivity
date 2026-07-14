## ============================================================================
## 03_predictors_future_to_tif.R
## Project: TR_Bear_Connectivity — ENMTML pipeline
## Purpose: Build 18 future scenario sub-folders ENMTML will consume.
##   For each (period × SSP × GCM) scenario:
##     (1) Read 19 Bio TIFs from Predictors_TIF/future/<scenario>/ (dynamic).
##     (2) Reproject + mask to target grid.
##     (3) Write to data/predictors_enmtml/future/<scenario>/Bio01..Bio19.tif.
##     (4) COPY 11 static layers (4 topo + 7 hum) from present ENMTML dir
##         into the same sub-folder (ENMTML needs all predictors in each sub-folder).
##   18 scenarios total: 2 periods × 3 SSPs × 3 GCMs.
## Inputs:
##   data/Predictors_TIF/future/<scenario>/Bio*.tif  (source)
##   data/predictors_enmtml/present/*.tif            (built by script 02)
##   outputs/rasters/target_grid.tif
##   outputs/rasters/tr_landmask.tif
## Outputs:
##   data/predictors_enmtml/future/<scenario>/*.tif  (30 layers per scenario)
##   outputs/figures/03_future/<scenario>_overview.png
## ============================================================================

suppressPackageStartupMessages({
  library(terra); library(sf); library(ggplot2); library(tidyterra); library(patchwork)
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
cat(sprintf("[bootstrap] wd = %s\n", getwd()))
source("00_paths.R"); source("00_helpers.R")
tb_log_init("03_predictors_future_to_tif")
tb_pkg_versions(c("terra","sf"))

## ============================================================================
## 1. LOAD TARGET GRID + STATIC PREDICTORS LIST
## ============================================================================
tb_log_section("1. LOAD TARGET + STATIC LIST")
if (!file.exists(file.path(TB_OUT_RASTERS, "target_grid.tif")))
  stop("target_grid.tif missing — run 02_predictors_present_to_tif.R first")

target  <- terra::rast(file.path(TB_OUT_RASTERS, "target_grid.tif"))
tr_mask <- terra::rast(file.path(TB_OUT_RASTERS, "tr_landmask.tif"))
grid_info <- readRDS(file.path(TB_OUT_RDS, "02_target_grid_info.rds"))
tr_proj <- grid_info$tr_polygon_proj

static_names <- c(TB_NAMES_TOP, TB_NAMES_HUM)
static_files <- file.path(TB_PRED_ENMTML_PRESENT, sprintf("%s.tif", static_names))
miss_static  <- static_names[!file.exists(static_files)]
if (length(miss_static))
  stop("static TIFs missing in TB_PRED_ENMTML_PRESENT (run 02 first): ",
       paste(miss_static, collapse=", "))
tb_log(sprintf("static layers to copy per scenario: %d", length(static_names)))

## scenario list: 18 combinations
SCENARIOS <- character(0)
for (p in TB_PERIODS) for (g in TB_GCMS) for (s in TB_SSPS) {
  SCENARIOS <- c(SCENARIOS, sprintf("%s_%s_%s", p, g, s))
}
tb_log(sprintf("%d scenarios to build", length(SCENARIOS)))

## ============================================================================
## 2. LOOP OVER SCENARIOS
## ============================================================================
tb_log_section("2. BUILD SCENARIO SUB-FOLDERS")
tb_tic()

for (scn in SCENARIOS) {
  tb_log_section(sprintf("scenario: %s", scn))
  out_dir <- file.path(TB_PRED_ENMTML_FUTURE, scn)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  src_dir <- file.path(TB_PREDICTORS_TIF_FUTURE, scn)
  if (!dir.exists(src_dir)) {
    tb_log(sprintf("  SKIP: source missing %s", src_dir), "WARN")
    next
  }

  ## --- dynamic: 19 Bio layers ---
  dyn_files <- file.path(src_dir, sprintf("%s.tif", TB_NAMES_BIO))
  miss_dyn  <- TB_NAMES_BIO[!file.exists(dyn_files)]
  if (length(miss_dyn)) {
    tb_log(sprintf("  SKIP: missing %d Bio TIFs: %s", length(miss_dyn),
                   paste(miss_dyn, collapse=",")), "WARN")
    next
  }
  for (nm in TB_NAMES_BIO) {
    src <- file.path(src_dir, sprintf("%s.tif", nm))
    r   <- terra::rast(src)
    if (terra::crs(r) == "") terra::crs(r) <- TB_CRS_WGS
    r   <- terra::project(r, target, method = "bilinear", threads = TRUE)
    r   <- terra::mask(r, tr_mask)
    names(r) <- nm
    fn <- file.path(out_dir, sprintf("%s.tif", nm))
    terra::writeRaster(r, fn, overwrite = TRUE, datatype = "FLT4S",
                       gdal = c("COMPRESS=DEFLATE","PREDICTOR=2","TILED=YES"))
  }
  tb_log(sprintf("  dynamic: %d Bio TIFs written", length(TB_NAMES_BIO)))

  ## --- static: copy 11 from present ENMTML dir ---
  for (nm in static_names) {
    src <- file.path(TB_PRED_ENMTML_PRESENT, sprintf("%s.tif", nm))
    dst <- file.path(out_dir, sprintf("%s.tif", nm))
    file.copy(src, dst, overwrite = TRUE)
  }
  tb_log(sprintf("  static: %d TIFs copied from present", length(static_names)))

  ## --- overview figure ---
  stk <- terra::rast(list.files(out_dir, pattern = "\\.tif$", full.names = TRUE))
  bbox_proj <- sf::st_bbox(tr_proj)
  mk_thumb <- function(r) {
    ggplot() +
      geom_spatraster(data = r, maxcell = 1.5e5) +
      scale_fill_viridis_c(option = "viridis", na.value = "transparent", guide = "none") +
      geom_sf(data = tr_proj, fill = NA, color = TB_COLOR_FRAME, linewidth = 0.25) +
      coord_sf(xlim = c(bbox_proj$xmin, bbox_proj$xmax),
               ylim = c(bbox_proj$ymin, bbox_proj$ymax), expand = FALSE) +
      theme_void() +
      theme(plot.title = element_text(size = 8, face = "bold",
                                      color = TB_COLOR_FRAME, hjust = 0.5)) +
      labs(title = names(r))
  }
  thumbs <- lapply(names(stk), function(nm) mk_thumb(stk[[nm]]))

  parts <- strsplit(scn, "_")[[1]]
  period_key <- sprintf("%s_%s", parts[1], parts[2])
  gcm <- paste(parts[3:(length(parts)-1)], collapse="_")
  ssp <- tail(parts, 1)
  period_label <- TB_PERIOD_LABELS[period_key] %||% period_key
  ssp_label    <- TB_SSP_LABELS[ssp] %||% ssp
  gcm_label    <- TB_GCM_LABELS[gcm] %||% gsub("_", "-", gcm)

  overview <- patchwork::wrap_plots(thumbs, ncol = 6) +
    patchwork::plot_annotation(
      title = sprintf("Future predictors — %s | %s | %s", period_label, gcm_label, ssp_label),
      theme = theme(
        plot.title = element_text(face = "bold", size = 14, color = TB_COLOR_FRAME))
    )
  tb_save_fig(overview, sprintf("03_future_%s_overview", scn),
              w = 18, h = 13, dpi = 200, subdir = "03_future")
}

tb_toc("all scenarios")

## ============================================================================
## 3. VERIFY
## ============================================================================
tb_log_section("3. VERIFY")
ok <- 0; bad <- 0
for (scn in SCENARIOS) {
  d <- file.path(TB_PRED_ENMTML_FUTURE, scn)
  if (!dir.exists(d)) { bad <- bad + 1; next }
  n <- length(list.files(d, pattern = "\\.tif$"))
  if (n == length(TB_NAMES_ALL)) ok <- ok + 1 else {
    bad <- bad + 1
    tb_log(sprintf("scenario %s has %d TIFs (expected %d)",
                   scn, n, length(TB_NAMES_ALL)), "WARN")
  }
}
tb_log(sprintf("verification: %d/%d scenarios OK", ok, length(SCENARIOS)))

tb_log_session()
tb_log("03_predictors_future_to_tif DONE")
