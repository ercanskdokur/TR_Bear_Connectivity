## ============================================================================
## 00_helpers.R
## Project: TR_Bear_Connectivity
## Purpose: Logging, theme, palettes, IO helpers. Source after 00_paths.R.
## ============================================================================

suppressPackageStartupMessages({
  for (.pkg in c("glue","ggplot2","ggspatial","scales","sf","terra")) {
    if (!requireNamespace(.pkg, quietly = TRUE)) {
      install.packages(.pkg, repos = "https://cloud.r-project.org")
    }
  }
  library(glue); library(ggplot2); library(ggspatial); library(scales)
})

## Null-coalescing operator (base R has it only since 4.4; we run on 4.3.2)
if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b
  assign("%||%", `%||%`, envir = .GlobalEnv)
}

## ---- Logging ---------------------------------------------------------------
.tb_log_file <- NULL

tb_log_init <- function(script_name) {
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  fn <- file.path(TB_OUT_LOGS, sprintf("%s_%s.log", script_name, ts))
  .tb_log_file <<- fn
  cat(sprintf("# log file: %s\n", fn), file = fn, append = FALSE)
  tb_log(sprintf("script: %s", script_name))
  tb_log(sprintf("R: %s | platform: %s | host: %s",
                 R.version.string, R.version$platform, Sys.info()[["nodename"]]))
  tb_log(sprintf("env: TB_ENV=%s | data=%s | out=%s", TB_ENV, TB_DATA_ROOT, TB_OUT_ROOT))
  invisible(fn)
}

tb_log <- function(msg, level = "INFO") {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] [%s] %s", ts, level, msg)
  message(line)
  if (!is.null(.tb_log_file)) cat(line, "\n", file = .tb_log_file, append = TRUE, sep = "")
  invisible(line)
}

tb_log_section <- function(title) {
  bar <- paste(rep("=", 78), collapse = "")
  tb_log(bar); tb_log(title); tb_log(bar)
}

tb_log_session <- function() {
  if (is.null(.tb_log_file)) return(invisible())
  cat("\n## sessionInfo()\n", file = .tb_log_file, append = TRUE)
  capture.output(sessionInfo(), file = .tb_log_file, append = TRUE)
}

tb_tic <- function() assign(".tb_tic_t", Sys.time(), envir = .GlobalEnv)
tb_toc <- function(msg = "elapsed") {
  if (!exists(".tb_tic_t", envir = .GlobalEnv)) return(invisible())
  dt <- as.numeric(difftime(Sys.time(), get(".tb_tic_t", envir = .GlobalEnv), units = "secs"))
  tb_log(sprintf("%s: %.1f s (%.2f min)", msg, dt, dt/60))
}

## ---- Safe IO ---------------------------------------------------------------
tb_save_rds <- function(obj, name) {
  fn <- file.path(TB_OUT_RDS, sprintf("%s.rds", name))
  saveRDS(obj, fn)
  tb_log(sprintf("saved RDS: %s (%.1f MB)", fn, file.info(fn)$size / 1e6))
  invisible(fn)
}

tb_read_rds <- function(name) {
  fn <- file.path(TB_OUT_RDS, sprintf("%s.rds", name))
  if (!file.exists(fn)) stop("RDS missing: ", fn)
  readRDS(fn)
}

tb_save_raster <- function(r, name, subdir = NULL, datatype = "FLT4S") {
  if (!requireNamespace("terra", quietly = TRUE)) stop("terra needed")
  d <- if (is.null(subdir)) TB_OUT_RASTERS else file.path(TB_OUT_RASTERS, subdir)
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  fn <- file.path(d, sprintf("%s.tif", name))
  terra::writeRaster(r, fn, overwrite = TRUE, datatype = datatype,
                     gdal = c("COMPRESS=DEFLATE","PREDICTOR=2","TILED=YES"))
  tb_log(sprintf("saved raster: %s", fn))
  invisible(fn)
}

tb_save_vector <- function(v, name, subdir = NULL) {
  if (!requireNamespace("sf", quietly = TRUE)) stop("sf needed")
  d <- if (is.null(subdir)) TB_OUT_VECTORS else file.path(TB_OUT_VECTORS, subdir)
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  fn <- file.path(d, sprintf("%s.gpkg", name))
  sf::st_write(v, fn, delete_dsn = TRUE, quiet = TRUE)
  tb_log(sprintf("saved vector: %s", fn))
  invisible(fn)
}

tb_save_table <- function(df, name, subdir = NULL) {
  d <- if (is.null(subdir)) TB_OUT_TABLES else file.path(TB_OUT_TABLES, subdir)
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  fn <- file.path(d, sprintf("%s.csv", name))
  write.csv(df, fn, row.names = FALSE)
  tb_log(sprintf("saved table: %s (%d rows)", fn, nrow(df)))
  invisible(fn)
}

tb_save_fig <- function(p, name, w = 14, h = 9, dpi = 600, subdir = NULL) {
  d <- if (is.null(subdir)) TB_OUT_FIGURES else file.path(TB_OUT_FIGURES, subdir)
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  fn <- file.path(d, sprintf("%s.png", name))
  ggsave(fn, plot = p, width = w, height = h, dpi = dpi, bg = "white")
  tb_log(sprintf("saved figure: %s", fn))
  invisible(fn)
}

## ---- Colorblind-friendly palettes ------------------------------------------
TB_PAL_OKABE_ITO <- c(
  "#000000","#E69F00","#56B4E9","#009E73","#F0E442",
  "#0072B2","#D55E00","#CC79A7"
)

TB_PAL_GAINLOSS <- c(
  "Loss"          = "#D55E00",
  "Stable unsuit" = "#BDBDBD",
  "Stable suit"   = "#56B4E9",
  "Gain"          = "#009E73"
)

TB_PAL_BINARY <- c("Unsuitable" = "#BDBDBD", "Suitable" = "#009E73")

TB_PAL_BIOREG <- c(
  "Euro-Siberian"  = "#0072B2",
  "Irano-Turanian" = "#E69F00",
  "Mediterranean"  = "#009E73"
)

TB_PAL_ACTIVITY <- c(
  "Beekeeping"                 = "#E69F00",
  "Livestock"                  = "#D55E00",
  "Farming Activity"           = "#F0E442",
  "Human Activity in Forest"   = "#009E73",
  "Road accident"              = "#CC79A7",
  "Foraging around settlement" = "#56B4E9",
  "Poaching"                   = "#000000"
)

TB_FILL_LAND    <- "#F2F2F2"
TB_COLOR_LAND   <- "#BDBDBD"
TB_FILL_SEA     <- "#E3F2FD"
TB_COLOR_FRAME  <- "#263238"
TB_COLOR_AXIS   <- "#455A64"

## ---- Theme -----------------------------------------------------------------
theme_trbear <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.background = element_rect(fill = TB_FILL_SEA, color = NA),
      panel.grid.major = element_line(color = "white", linewidth = 0.2),
      panel.grid.minor = element_blank(),
      axis.text  = element_text(color = TB_COLOR_AXIS, size = base_size - 2),
      axis.title = element_blank(),
      plot.title    = element_text(face = "bold", size = base_size + 4, color = TB_COLOR_FRAME),
      plot.subtitle = element_text(size = base_size, color = TB_COLOR_AXIS),
      plot.caption  = element_text(size = base_size - 3, color = TB_COLOR_AXIS, hjust = 1),
      panel.border  = element_rect(color = TB_COLOR_FRAME, fill = NA, linewidth = 1.2),
      legend.background = element_rect(fill = "white", color = NA),
      legend.title  = element_text(face = "bold", size = base_size - 1),
      plot.margin   = ggplot2::margin(t = 20, r = 20, b = 20, l = 20)
    )
}

theme_trbear_bar <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major.y = element_line(color = "gray85", linewidth = 0.3),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text  = element_text(color = TB_COLOR_AXIS),
      axis.title = element_text(color = TB_COLOR_FRAME, face = "bold"),
      plot.title    = element_text(face = "bold", size = base_size + 4, color = TB_COLOR_FRAME),
      plot.subtitle = element_text(size = base_size, color = TB_COLOR_AXIS),
      panel.border  = element_rect(color = TB_COLOR_FRAME, fill = NA, linewidth = 1.0),
      legend.position = "right",
      legend.title    = element_text(face = "bold", size = base_size - 1)
    )
}

tb_map_decorations <- function(north = "tr", scale = "bl") {
  list(
    annotation_north_arrow(
      location = north, which_north = "true",
      pad_x = unit(0.3, "in"), pad_y = unit(0.3, "in"),
      style = north_arrow_minimal(line_col = TB_COLOR_FRAME, text_col = TB_COLOR_FRAME)
    ),
    annotation_scale(
      location = scale, width_hint = 0.18, style = "ticks",
      text_col = TB_COLOR_FRAME, line_col = TB_COLOR_FRAME
    )
  )
}

tb_basemap_world <- function(crs_target = TB_CRS_PROJ, scale_res = "medium") {
  if (!requireNamespace("rnaturalearth", quietly = TRUE)) {
    install.packages("rnaturalearth", repos = "https://cloud.r-project.org")
  }
  w <- rnaturalearth::ne_countries(scale = scale_res, returnclass = "sf")
  sf::st_transform(w, crs_target)
}

## ---- Utility ---------------------------------------------------------------
tb_to_proj <- function(x) {
  if (inherits(x, "sf"))         return(sf::st_transform(x, TB_CRS_PROJ))
  if (inherits(x, "SpatRaster")) return(terra::project(x, TB_CRS_PROJ))
  if (inherits(x, "SpatVector")) return(terra::project(x, TB_CRS_PROJ))
  stop("tb_to_proj: unsupported class")
}

tb_pkg_versions <- function(pkgs = c("terra","sf","landscapemetrics")) {
  v <- sapply(pkgs, function(p)
    tryCatch(as.character(packageVersion(p)), error = function(e) NA))
  tb_log(paste(names(v), v, sep = "=", collapse = " | "))
}

## ---- Spatial thinning -------------------------------------------------------
## Drop-in replacement for spThin: keep one point per grid cell of a reference
## raster (or one per radius using distance-based filtering).
##
## METHOD = "grid": rasterize-based, deterministic given seed, O(n). Recommended
##                  for 1-km thinning on a known prediction grid.
## METHOD = "dist": iterative pairwise distance filter (spThin-style). O(n^2)
##                  for n points; fine for n < 10k.

tb_thin_grid <- function(points_sf, ref_raster, seed = 42L) {
  stopifnot(inherits(points_sf, "sf"), inherits(ref_raster, "SpatRaster"))
  if (sf::st_crs(points_sf) != sf::st_crs(ref_raster)) {
    points_sf <- sf::st_transform(points_sf, terra::crs(ref_raster))
  }
  pv <- terra::vect(points_sf)
  cells <- terra::cells(ref_raster[[1]], pv)[, "cell"]
  set.seed(seed)
  ord  <- sample(seq_along(cells))
  keep <- !duplicated(cells[ord])
  idx  <- ord[keep]
  out  <- points_sf[idx, , drop = FALSE]
  attr(out, "tb_thin_method") <- "grid"
  attr(out, "tb_thin_n_in")   <- nrow(points_sf)
  attr(out, "tb_thin_n_out")  <- nrow(out)
  tb_log(sprintf("tb_thin_grid: %d -> %d points (cell-based, seed=%d)",
                 nrow(points_sf), nrow(out), seed))
  out
}

tb_thin_dist <- function(points_sf, min_dist_m, n_reps = 10L, seed = 42L) {
  stopifnot(inherits(points_sf, "sf"))
  set.seed(seed)
  best <- NULL
  for (r in seq_len(n_reps)) {
    ord <- sample(seq_len(nrow(points_sf)))
    p   <- points_sf[ord, , drop = FALSE]
    coords <- sf::st_coordinates(p)
    keep <- logical(nrow(p))
    keep[1] <- TRUE
    kept_coords <- coords[1, , drop = FALSE]
    for (i in seq.int(2, nrow(p))) {
      d <- sqrt(rowSums((sweep(kept_coords, 2, coords[i, ]))^2))
      if (all(d >= min_dist_m)) {
        keep[i] <- TRUE
        kept_coords <- rbind(kept_coords, coords[i, , drop = FALSE])
      }
    }
    out_r <- p[keep, , drop = FALSE]
    if (is.null(best) || nrow(out_r) > nrow(best)) best <- out_r
  }
  attr(best, "tb_thin_method")    <- "dist"
  attr(best, "tb_thin_min_dist")  <- min_dist_m
  attr(best, "tb_thin_n_in")      <- nrow(points_sf)
  attr(best, "tb_thin_n_out")     <- nrow(best)
  tb_log(sprintf("tb_thin_dist: %d -> %d points (min_dist=%.0f m, reps=%d)",
                 nrow(points_sf), nrow(best), min_dist_m, n_reps))
  best
}
