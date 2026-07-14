## ============================================================================
## 13_unicor_prep.R
## Project: TR_Bear_Connectivity
## Purpose: Prepare UNICOR inputs for present + 6 future scenarios.
##   - Generate source points from PRESENT binary suitable raster:
##       connected patches ≥ TB_PATCH_MIN_KM2 → centroid
##       (shared across all scenarios to make connectivity comparable)
##   - Convert each resistance .tif → .rsg (Arc/Info ASCII grid)
##   - Write per-scenario .rip config (UNICOR parameter file)
##
## Outputs (derived/unicor/):
##   sources.xy                    (shared source points, X,Y header)
##   sources.gpkg                  (sf version for plotting / QC)
##   <scenario>/<scenario>.rsg     (resistance grid in ASCII)
##   <scenario>/<scenario>.rip     (UNICOR config)
##   <scenario>/sources.xy         (copy of shared sources for UNICOR)
##   tables/13_source_points.csv   (patch ID, area, centroid X/Y)
##   figures/13_unicor/fig13a_source_points.png
##                     fig13b_source_patches.png
## ============================================================================

suppressPackageStartupMessages({
  library(terra); library(sf); library(ggplot2); library(dplyr)
  library(tidyterra); library(rnaturalearth)
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
tb_log_init("13_unicor_prep")

FIG_SUBDIR <- "13_unicor"

UNICOR_PROCS        <- 16L
UNICOR_EDGE_DIST    <- 500000L     ## permissive — refine after first run
UNICOR_KDE_GRIDSIZE <- 5L

## ----------------------------------------------------------------------------
## 1) Generate shared source points from PRESENT binary suitable raster
## ----------------------------------------------------------------------------
tb_log_section("Source points from present binary")

pres_bin <- terra::rast(file.path(TB_OUT_HS_BINARY, "present_wmean.tif"))
names(pres_bin) <- "suit"

cell_area_km2 <- prod(terra::res(pres_bin)) / 1e6
tb_log(sprintf("cell area = %.4f km²", cell_area_km2))

## terra::patches: connected components (8-direction adjacency)
patches <- terra::patches(pres_bin, directions = 8, zeroAsNA = TRUE)
names(patches) <- "patch"

## Compute area per patch
patch_stats <- terra::freq(patches) |>
  as.data.frame() |>
  rename(patch_id = value, n_cells = count) |>
  mutate(area_km2 = n_cells * cell_area_km2) |>
  filter(!is.na(patch_id)) |>
  arrange(desc(area_km2))

tb_log(sprintf("total patches: %d | ≥%d km²: %d",
               nrow(patch_stats),
               TB_PATCH_MIN_KM2,
               sum(patch_stats$area_km2 >= TB_PATCH_MIN_KM2)))

## Keep patches ≥ threshold
keep_patches <- patch_stats$patch_id[patch_stats$area_km2 >= TB_PATCH_MIN_KM2]
patches_keep <- terra::ifel(patches %in% keep_patches, patches, NA)
names(patches_keep) <- "patch"

## Centroid of each kept patch (mean of cell coordinates)
get_centroid <- function(patches_r, ids) {
  out <- lapply(ids, function(i) {
    p <- terra::ifel(patches_r == i, 1, NA)
    pts <- terra::as.points(p)
    if (nrow(pts) == 0) return(NULL)
    crd <- terra::crds(pts)
    data.frame(patch_id = i,
               n_cells  = nrow(crd),
               x        = mean(crd[, 1]),
               y        = mean(crd[, 2]))
  })
  do.call(rbind, out)
}

source_df <- get_centroid(patches_keep, keep_patches)
source_df$area_km2 <- source_df$n_cells * cell_area_km2
source_df <- source_df |> arrange(desc(area_km2))

tb_log(sprintf("source points: %d (largest %.0f km², smallest %.0f km²)",
               nrow(source_df), max(source_df$area_km2), min(source_df$area_km2)))

tb_save_table(source_df, "13_source_points")

## Save XY file (UNICOR format: "X,Y" header, no patch ID)
xy_master <- file.path(TB_OUT_UNICOR_DIR, "sources.xy")
write.table(
  data.frame(X = source_df$x, Y = source_df$y),
  xy_master,
  sep = ",", row.names = FALSE, quote = FALSE)
tb_log(sprintf("wrote %s (%d points)", xy_master, nrow(source_df)))

## Sf version for QC
sources_sf <- sf::st_as_sf(source_df, coords = c("x", "y"), crs = TB_CRS_PROJ)
tb_save_vector(sources_sf, "sources", subdir = NULL)

## ----------------------------------------------------------------------------
## 2) Per-scenario UNICOR input files
## ----------------------------------------------------------------------------
tb_log_section("Per-scenario UNICOR inputs")

scenarios <- c("present",
               sprintf("%s_%s", rep(TB_PERIODS, each = length(TB_SSPS)),
                                rep(TB_SSPS,    times = length(TB_PERIODS))))

write_rip <- function(rip_path, session_label, rsg_basename, xy_basename) {
  txt <- c(
    sprintf("Session_label\t%s", session_label),
    sprintf("Grid_Filename\t%s", rsg_basename),
    sprintf("XY_Filename\t%s", xy_basename),
    "Use_Direction\tFALSE",
    "Type_Direction\tFlowAcc",
    "Use_Resistance\tTRUE",
    "Barrier_or_U_Filename\tNA",
    "Direction_or_V_Filename\tNA",
    "Speed_To_Resistance_Scale\t0;10",
    "Use_ED_threshold\tFalse",
    sprintf("ED_Distance\t%d", UNICOR_EDGE_DIST),
    "Edge_Type\tnormal",
    "Transform_function\tlinear",
    "Const_kernal_vol\tFalse",
    "Kernel_volume\t10000",
    sprintf("Edge_Distance\t%d", UNICOR_EDGE_DIST),
    sprintf("Number_of_Processes\t%d", UNICOR_PROCS),
    "KDE_Function\tGaussian",
    sprintf("KDE_GridSize\t%d", UNICOR_KDE_GRIDSIZE),
    "Number_of_Categories\t5",
    "Save_Path_Output\tTRUE",
    "Save_IndividualPaths_Output\tFALSE",
    "Save_GraphMetrics_Output\tTRUE",
    "Save_KDE_Output\tTRUE",
    "Save_Category_Output\tTRUE",
    "Save_CDmatrix_Output\tTRUE"
  )
  writeLines(txt, rip_path)
}

for (scen in scenarios) {
  out_dir <- file.path(TB_OUT_UNICOR_DIR, scen)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  src_tif <- file.path(TB_OUT_RESISTANCE, sprintf("%s.tif", scen))
  if (!file.exists(src_tif)) {
    tb_log(sprintf("MISSING resistance %s", src_tif), "WARN"); next
  }
  R <- terra::rast(src_tif)
  if (terra::nlyr(R) > 1) R <- R[[1]]

  ## Replace NA with NODATA sentinel that UNICOR expects (-9999), via AAIGrid
  rsg_path <- file.path(out_dir, sprintf("%s.rsg", scen))
  terra::writeRaster(R, rsg_path,
                     overwrite = TRUE,
                     filetype  = "AAIGrid",
                     NAflag    = -9999,
                     datatype  = "FLT4S")
  tb_log(sprintf("wrote %s", rsg_path))

  ## Copy XY into scenario folder (UNICOR uses basename relative to .rip)
  file.copy(xy_master, file.path(out_dir, "sources.xy"), overwrite = TRUE)

  ## Write RIP
  write_rip(rip_path  = file.path(out_dir, sprintf("%s.rip", scen)),
            session_label = scen,
            rsg_basename  = sprintf("%s.rsg", scen),
            xy_basename   = "sources.xy")
  tb_log(sprintf("wrote %s.rip", scen))
}

## ----------------------------------------------------------------------------
## FIGURES
## ----------------------------------------------------------------------------
tb_log_section("Figures")

world_sf <- tryCatch(
  rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") |>
    sf::st_transform(TB_CRS_PROJ),
  error = function(e) NULL)
tr_mask_shp <- file.path(TB_DATA_ROOT, "TR_mask", "TR_mask.shp")
tr_mask_sf <- if (file.exists(tr_mask_shp))
  sf::st_read(tr_mask_shp, quiet = TRUE) |> sf::st_transform(TB_CRS_PROJ) else NULL

e <- terra::ext(pres_bin); ext_pad <- 30000
xlim_p <- c(e$xmin - ext_pad, e$xmax + ext_pad)
ylim_p <- c(e$ymin - ext_pad, e$ymax + ext_pad)

.clip_to_tr <- function(r) {
  if (is.null(tr_mask_sf)) return(r)
  terra::mask(r, terra::vect(tr_mask_sf))
}

## ---- fig13a: source points on present binary HS -----------------------------
pres_bin_factor <- .clip_to_tr(pres_bin)
levels(pres_bin_factor) <- data.frame(value = c(0, 1),
                                       class = c("Unsuitable", "Suitable"))

p13a <- ggplot()
if (!is.null(world_sf)) p13a <- p13a +
  geom_sf(data = world_sf, fill = "#E8E8E8",
          color = "#7C8A93", linewidth = 0.4)
p13a <- p13a +
  tidyterra::geom_spatraster(data = pres_bin_factor, na.rm = TRUE) +
  scale_fill_manual(values = TB_PAL_BINARY, na.translate = FALSE,
                    name = "Habitat")
if (!is.null(tr_mask_sf)) p13a <- p13a +
  geom_sf(data = tr_mask_sf, fill = NA,
          color = TB_COLOR_FRAME, linewidth = 0.5)
p13a <- p13a +
  geom_sf(data = sources_sf,
          aes(size = area_km2),
          fill = "#D55E00", color = "#000000",
          shape = 21, stroke = 0.5, alpha = 0.85) +
  scale_size_continuous(range = c(2, 8),
                        breaks = c(100, 500, 1000, 5000, 10000),
                        labels = scales::label_comma(),
                        name = "Patch area\n(km²)") +
  coord_sf(xlim = xlim_p, ylim = ylim_p, datum = sf::st_crs(4326),
           expand = FALSE) +
  tb_map_decorations() +
  labs(title    = sprintf("Habitat patches (source nodes) for connectivity analysis (n=%d)",
                          nrow(source_df)),
       subtitle = sprintf(
         "Centroids of suitable habitat patches ≥ %d km² in the present-day binary HS map (W_MEAN MAX_TSS).",
         TB_PATCH_MIN_KM2)) +
  theme_trbear()

tb_save_fig(p13a, "fig13a_source_points", w = 14, h = 9, subdir = FIG_SUBDIR)

## ---- fig13b: kept patches coloured by AREA (km²) ----------------------------
## Earlier this panel coloured polygons by raw patch ID, which carried no
## meaning (assignment order from terra::patches()).  Now we recolour by
## patch area — a meaningful conservation metric — and expose the legend.
patches_keep_clip <- .clip_to_tr(patches_keep)
area_lookup <- setNames(patch_stats$area_km2[match(keep_patches,
                                                     patch_stats$patch_id)],
                         as.character(keep_patches))
patches_area <- terra::classify(patches_keep_clip,
                                rcl = cbind(as.numeric(names(area_lookup)),
                                             as.numeric(area_lookup)))
names(patches_area) <- "area_km2"

p13b <- ggplot()
if (!is.null(world_sf)) p13b <- p13b +
  geom_sf(data = world_sf, fill = "#E8E8E8",
          color = "#7C8A93", linewidth = 0.4)
p13b <- p13b +
  tidyterra::geom_spatraster(data = patches_area, na.rm = TRUE) +
  scale_fill_viridis_c(option = "viridis", trans = "log10",
                       na.value = "transparent",
                       name = "Patch area\n(km², log10)",
                       breaks = c(100, 500, 1000, 5000, 10000, 40000),
                       labels = scales::label_comma())
if (!is.null(tr_mask_sf)) p13b <- p13b +
  geom_sf(data = tr_mask_sf, fill = NA,
          color = TB_COLOR_FRAME, linewidth = 0.5)
p13b <- p13b +
  geom_sf(data = sources_sf, fill = "white", color = "#000000",
          shape = 21, size = 1.5, stroke = 0.4) +
  coord_sf(xlim = xlim_p, ylim = ylim_p, datum = sf::st_crs(4326),
           expand = FALSE) +
  tb_map_decorations() +
  labs(title    = sprintf("Habitat patches (source nodes), n=%d",
                          nrow(source_df)),
       subtitle = sprintf(
         "Coloured polygons: connected suitable patches ≥ %d km², shaded by area on a log scale.  White dots: patch centroids.",
         TB_PATCH_MIN_KM2)) +
  theme_trbear()

tb_save_fig(p13b, "fig13b_source_patches", w = 14, h = 9, subdir = FIG_SUBDIR)

tb_log_session()
tb_log("13_unicor_prep DONE")
