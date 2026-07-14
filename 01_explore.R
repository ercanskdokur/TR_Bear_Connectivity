## ============================================================================
## 01_explore_data.R
## Project: TR_Bear_Connectivity
## Purpose: Sanity-check all input data: points, predictors, PA .gdb, roads.
##          Produces a structured summary log + a QC table + one quick TR map.
## Inputs:
##   data/points/PresencePoints.txt
##   data/points/ConflictPoints.txt
##   data/points/conflict_table_df.xlsx
##   data/Predictors_TIF/present/*.tif         (30 present predictor layers)
##   data/ProtectedAreas/PAs/PAs.gpkg
##   data/Roads/gis_osm_roads_free_1.shp
## Outputs:
##   outputs/tables/01_data_summary.csv
##   outputs/tables/01_predictor_stats_present.csv
##   outputs/tables/01_conflict_typology.csv
##   outputs/rds_files/01_explore_summary.rds
##   outputs/figures/01_qc/01_qc_points_overview.png
##   outputs/logs/01_explore_data_<ts>.log
## ============================================================================

suppressPackageStartupMessages({
  library(terra); library(sf); library(dplyr); library(readr); library(readxl)
  library(ggplot2); library(rnaturalearth)
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
  stop("Cannot find 00_paths.R in any known location")
}
setwd(.tb_find_paths_R())
cat(sprintf("[bootstrap] wd = %s\n", getwd()))
source("00_paths.R")
source("00_helpers.R")
tb_log_init("01_explore_data")
tb_pkg_versions(c("terra","sf","blockCV","landscapemetrics","readxl"))

summary_rows <- list()
add_row <- function(item, value, status = "OK", note = "") {
  v <- tryCatch(as.character(value), error = function(e) NA_character_)
  if (length(v) == 0L) v <- NA_character_
  if (length(v) > 1L)  v <- paste(v, collapse = ";")
  summary_rows[[length(summary_rows) + 1]] <<-
    data.frame(item = item, value = v,
               status = status, note = note, stringsAsFactors = FALSE)
}

## ============================================================================
## 1. POINTS — Presence
## ============================================================================
tb_log_section("1. PRESENCE POINTS")
tb_tic()
pres <- read.table(TB_PRESENCE_TXT, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, quote = "\"")
tb_log(sprintf("rows = %d, cols = %s", nrow(pres), paste(names(pres), collapse=", ")))
tb_log(sprintf("species values: %s",
               paste(sort(unique(pres$species)), collapse=", ")))
tb_log(sprintf("x range: [%.4f, %.4f]", min(pres$x), max(pres$x)))
tb_log(sprintf("y range: [%.4f, %.4f]", min(pres$y), max(pres$y)))
tb_log(sprintf("NA: x=%d, y=%d", sum(is.na(pres$x)), sum(is.na(pres$y))))
dup_xy <- sum(duplicated(pres[, c("x","y")]))
tb_log(sprintf("duplicate (x,y) pairs: %d", dup_xy))
add_row("presence_n", nrow(pres))
add_row("presence_duplicates_xy", dup_xy,
        if (dup_xy > 0) "WARN" else "OK",
        "to be thinned in 06_points_prep")
tb_toc("presence read")

## ============================================================================
## 2. POINTS — Conflict (.txt)
## ============================================================================
tb_log_section("2. CONFLICT POINTS (.txt)")
conf_txt <- read.table(TB_CONFLICT_TXT, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE, quote = "\"")
tb_log(sprintf("rows = %d, cols = %s", nrow(conf_txt), paste(names(conf_txt), collapse=", ")))
tb_log(sprintf("x range: [%.4f, %.4f]", min(conf_txt$x), max(conf_txt$x)))
tb_log(sprintf("y range: [%.4f, %.4f]", min(conf_txt$y), max(conf_txt$y)))
add_row("conflict_txt_n", nrow(conf_txt))

## ============================================================================
## 3. POINTS — Conflict typology (.xlsx)
## ============================================================================
tb_log_section("3. CONFLICT TABLE (xlsx)")
conf_xlsx <- readxl::read_excel(TB_CONFLICT_XLSX)
tb_log(sprintf("rows = %d, cols = %d", nrow(conf_xlsx), ncol(conf_xlsx)))
tb_log(sprintf("cols: %s", paste(names(conf_xlsx), collapse=", ")))

## normalize Biogeographic_Region spelling variants
br_raw <- table(conf_xlsx$Biogeographic_Region, useNA = "ifany")
tb_log("Biogeographic_Region raw values:")
for (k in names(br_raw)) tb_log(sprintf("  %s -> %d", k, br_raw[[k]]))

br_norm <- conf_xlsx$Biogeographic_Region
br_norm <- gsub("^Eurosiberian$", "Euro-Siberian", br_norm)
br_norm <- gsub("^IranoTuranian$", "Irano-Turanian", br_norm)
conf_xlsx$Biogeographic_Region_norm <- br_norm
tb_log("Biogeographic_Region normalized values:")
for (k in names(table(br_norm))) tb_log(sprintf("  %s -> %d", k, table(br_norm)[[k]]))

## Activity / Sub_reason / Behaviour
act_tab <- table(conf_xlsx$Activity, useNA = "ifany")
tb_log("Activity counts:")
for (k in names(sort(act_tab, decreasing = TRUE))) tb_log(sprintf("  %s -> %d", k, act_tab[[k]]))
sub_tab <- table(conf_xlsx$Sub_reason, useNA = "ifany")
tb_log(sprintf("Sub_reason: %d unique values", length(sub_tab)))
beh_tab <- table(conf_xlsx$Behaviour, useNA = "ifany")
tb_log(sprintf("Behaviour: Foraging=%d, Non-Foraging=%d",
               beh_tab["Foraging"] %||% 0, beh_tab["Non-Foraging"] %||% 0))

## Damage / casualties
sum_fin  <- sum(conf_xlsx$Financial_ == 1, na.rm = TRUE)
sum_hdth <- sum(conf_xlsx$Human_deat == 1, na.rm = TRUE)
sum_bdth <- sum(conf_xlsx$Bear_death >= 1, na.rm = TRUE)
sum_hinj <- sum(conf_xlsx$Human_inju >= 1, na.rm = TRUE)
sum_binj <- sum(conf_xlsx$Bear_injur >= 1, na.rm = TRUE)
tb_log(sprintf("Damages: financial=%d, human_death=%d, bear_death=%d, human_injury=%d, bear_injury=%d",
               sum_fin, sum_hdth, sum_bdth, sum_hinj, sum_binj))
add_row("conflict_xlsx_n", nrow(conf_xlsx))
add_row("conflict_human_deaths", sum_hdth)
add_row("conflict_bear_deaths", sum_bdth)

## save conflict typology summary
typology <- data.frame(
  Activity = names(act_tab),
  n = as.integer(act_tab),
  stringsAsFactors = FALSE
)
typology <- typology[order(-typology$n), ]
tb_save_table(typology, "01_conflict_typology")

## sanity: do xlsx Long/Lat match conf_txt x/y? (478 vs 478)
add_row("conflict_xlsx_matches_txt", nrow(conf_xlsx) == nrow(conf_txt),
        if (nrow(conf_xlsx) == nrow(conf_txt)) "OK" else "WARN")

## ============================================================================
## 4. PREDICTORS — Present (TIF source)
## ============================================================================
tb_log_section("4. PREDICTORS — Present TIFs")
tb_log(sprintf("source dir: %s", TB_PREDICTORS_TIF_PRESENT))

tif_files <- list.files(TB_PREDICTORS_TIF_PRESENT, pattern = "\\.tif$", full.names = TRUE)
tb_log(sprintf("found %d TIFs", length(tif_files)))
expected_names <- TB_NAMES_ALL
found_names <- tools::file_path_sans_ext(basename(tif_files))
missing_tifs <- setdiff(expected_names, found_names)
extra_tifs   <- setdiff(found_names, expected_names)
if (length(missing_tifs)) tb_log(sprintf("MISSING TIFs: %s", paste(missing_tifs, collapse=", ")), "WARN")
if (length(extra_tifs))   tb_log(sprintf("UNEXPECTED TIFs: %s", paste(extra_tifs, collapse=", ")), "WARN")
add_row("present_tifs_found",   length(tif_files))
add_row("present_tifs_missing", length(missing_tifs),
        if (length(missing_tifs) == 0) "OK" else "WARN",
        paste(missing_tifs, collapse=";"))

predictor_stats <- list()
report_layer <- function(group, fn) {
  r  <- terra::rast(fn)
  nm <- names(r)
  mm <- terra::minmax(r)
  n_nonna <- terra::global(r, fun = "notNA")[[1]]
  n_na    <- terra::global(r, fun = "isNA")[[1]]
  tb_log(sprintf("  [%s] %s: [%.3f, %.3f], notNA=%s, NA=%s, dim=%s, res=%s",
                 group, nm, mm[1], mm[2], n_nonna, n_na,
                 paste(dim(r)[1:2], collapse="x"),
                 paste(round(res(r), 6), collapse="x")))
  predictor_stats[[length(predictor_stats) + 1]] <<- data.frame(
    stack = group, layer = nm,
    min = mm[1], max = mm[2],
    n_nonNA = n_nonna, n_NA = n_na,
    stringsAsFactors = FALSE
  )
}

group_of <- function(nm) {
  if (nm %in% TB_NAMES_TOP) "rTop"
  else if (nm %in% TB_NAMES_BIO) "rBio"
  else if (nm %in% TB_NAMES_HUM) "rHum"
  else "other"
}
for (nm in c(TB_NAMES_TOP, TB_NAMES_BIO, TB_NAMES_HUM)) {
  fn <- file.path(TB_PREDICTORS_TIF_PRESENT, sprintf("%s.tif", nm))
  if (file.exists(fn)) report_layer(group_of(nm), fn)
  else tb_log(sprintf("  [skip] %s: file missing", nm), "WARN")
}

## quick CRS/extent consistency check (use Bio01 + Elevation as anchors)
anchor <- terra::rast(file.path(TB_PREDICTORS_TIF_PRESENT, "Bio01.tif"))
tb_log(sprintf("anchor CRS: %s", substr(terra::crs(anchor, proj=TRUE), 1, 80)))
tb_log(sprintf("anchor extent: %s",
               paste(round(as.vector(terra::ext(anchor)), 4), collapse=", ")))

pstats_df <- do.call(rbind, predictor_stats)
tb_save_table(pstats_df, "01_predictor_stats_present")
add_row("present_predictors_total", nrow(pstats_df))

## ============================================================================
## 5. PROTECTED AREAS .gdb — list layers (don't merge here)
## ============================================================================
tb_log_section("5. PROTECTED AREAS .gdb LAYERS")
pa_layers <- tryCatch(sf::st_layers(TB_PA_GDB), error = function(e) NULL)
if (is.null(pa_layers)) {
  tb_log("FAILED to read PA .gdb", "ERROR")
  add_row("pa_gdb_readable", FALSE, "FAIL")
} else {
  tb_log(sprintf("layers in .gdb (%d total):", length(pa_layers$name)))
  for (i in seq_along(pa_layers$name)) {
    tb_log(sprintf("  [%2d] %s | geom=%s | features=%s",
                   i, pa_layers$name[i],
                   pa_layers$geomtype[[i]] %||% "?",
                   pa_layers$features[i] %||% "?"))
  }
  add_row("pa_gdb_n_layers", length(pa_layers$name))
}

## ============================================================================
## 6. ROADS shapefile
## ============================================================================
tb_log_section("6. ROADS shapefile")
roads_info <- tryCatch(sf::st_layers(TB_ROADS_SHP), error = function(e) NULL)
if (!is.null(roads_info)) {
  tb_log(sprintf("layer: %s | geom=%s | features=%s",
                 roads_info$name[1], roads_info$geomtype[[1]] %||% "?",
                 roads_info$features[1] %||% "?"))
  add_row("roads_features", roads_info$features[1] %||% NA_character_)
} else {
  tb_log("FAILED to read roads shp", "ERROR")
  add_row("roads_readable", FALSE, "FAIL")
}

## peek at attribute names (first 5 features only, to avoid loading whole file)
roads_peek <- tryCatch(
  sf::st_read(TB_ROADS_SHP, query = sprintf("SELECT * FROM \"%s\" LIMIT 5",
                                            tools::file_path_sans_ext(basename(TB_ROADS_SHP))),
              quiet = TRUE),
  error = function(e) NULL
)
if (!is.null(roads_peek)) {
  tb_log(sprintf("roads attribute names: %s", paste(names(roads_peek), collapse=", ")))
  if ("fclass" %in% names(roads_peek)) {
    tb_log(sprintf("roads fclass sample: %s",
                   paste(unique(roads_peek$fclass), collapse=", ")))
  }
}

## ============================================================================
## 7. QUICK QC FIGURE — TR + points (WGS84, no reproject yet)
## ============================================================================
tb_log_section("7. QC FIGURE: points on TR (WGS84)")
tb_tic()
tr <- tryCatch(
  rnaturalearth::ne_countries(country = "Turkey", scale = "medium", returnclass = "sf"),
  error = function(e) NULL
)

pres_sf <- sf::st_as_sf(pres, coords = c("x","y"), crs = TB_CRS_WGS)
conf_sf <- sf::st_as_sf(conf_txt, coords = c("x","y"), crs = TB_CRS_WGS)

p <- ggplot()
if (!is.null(tr)) {
  p <- p + geom_sf(data = tr, fill = TB_FILL_LAND, color = TB_COLOR_LAND, linewidth = 0.4)
}
p <- p +
  geom_sf(data = pres_sf, color = "#0072B2", alpha = 0.55, size = 0.9) +
  geom_sf(data = conf_sf, color = "#D55E00", alpha = 0.85, size = 1.2, shape = 17) +
  theme_trbear() +
  tb_map_decorations() +
  labs(
    title = "Bear presence (blue) & conflict (orange) — QC overview",
    subtitle = sprintf("Presence n=%d  |  Conflict n=%d  |  CRS: WGS84 (raw)",
                       nrow(pres_sf), nrow(conf_sf)),
    caption = "Source: PresencePoints.txt + ConflictPoints.txt"
  )
tb_save_fig(p, "01_qc_points_overview", w = 13, h = 8, subdir = "01_qc")
tb_toc("qc figure")

## ============================================================================
## SAVE SUMMARY
## ============================================================================
tb_log_section("WRITING SUMMARY")
summary_df <- do.call(rbind, summary_rows)
tb_save_table(summary_df, "01_data_summary")
tb_save_rds(list(
  presence_n        = nrow(pres),
  conflict_n        = nrow(conf_txt),
  conflict_xlsx_cols= names(conf_xlsx),
  current_layers    = list(rTop = TB_NAMES_TOP, rBio = TB_NAMES_BIO, rHum = TB_NAMES_HUM),
  pa_layers         = pa_layers$name,
  conflict_typology = typology,
  predictor_stats   = pstats_df
), "01_explore_summary")

tb_log_session()
tb_log("01_explore_data DONE")
