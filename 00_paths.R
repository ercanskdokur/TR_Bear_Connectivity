## ============================================================================
## 00_paths.R
## Project: TR_Bear_Connectivity (Turkey-wide brown bear SDM + UNICOR connectivity)
## Pipeline: ENMTML-based
## Purpose: Single source of truth for paths and parameters. Auto-detects
##          local vs cluster environment.
## ============================================================================

suppressPackageStartupMessages({
  if (!requireNamespace("fs", quietly = TRUE)) install.packages("fs", repos = "https://cloud.r-project.org")
})

.tb_detect_env <- function() {
  host <- tolower(Sys.info()[["nodename"]])
  if (nzchar(Sys.getenv("SLURM_JOB_ID"))) return("cluster")
  if (grepl("kuacc|login|node|compute|hpc", host)) return("cluster")
  "local"
}

TB_ENV <- .tb_detect_env()

## ---- Project location ------------------------------------------------------
## All paths resolve relative to the project root: the directory that holds
## these scripts and, by default, the data/ and outputs/ sub-folders. No
## absolute paths are hard-coded, so the archive runs from wherever it is
## unpacked. The runner scripts setwd() to their own folder before sourcing
## this file, so getwd() below is that folder.
##
## Every location can be overridden with an environment variable, which is the
## recommended way to reproduce an HPC layout where data, outputs and programs
## live on separate volumes:
##   TB_ROOT       project root       (default: the scripts' own directory)
##   TB_DATA_ROOT  input data root    (default: <TB_ROOT>/data)
##   TB_OUT_ROOT   output root        (default: <TB_ROOT>/outputs)
##   TB_PROGRAMS   external programs  (default: <TB_ROOT>/programs; UNICOR, Conefor)
TB_HPC_ROOT  <- Sys.getenv("TB_ROOT",      unset = getwd())
TB_DATA_ROOT <- Sys.getenv("TB_DATA_ROOT", unset = file.path(TB_HPC_ROOT, "data"))
TB_OUT_ROOT  <- Sys.getenv("TB_OUT_ROOT",  unset = file.path(TB_HPC_ROOT, "outputs"))
TB_SCRIPTS   <- Sys.getenv("TB_SCRIPTS",   unset = TB_HPC_ROOT)
TB_PROGRAMS  <- Sys.getenv("TB_PROGRAMS",  unset = file.path(TB_HPC_ROOT, "programs"))

## ---- Raw input paths (read-only) -------------------------------------------
TB_POINTS_DIR     <- file.path(TB_DATA_ROOT, "points")
TB_PRESENCE_TXT   <- file.path(TB_POINTS_DIR, "PresencePoints.txt")
TB_CONFLICT_TXT   <- file.path(TB_POINTS_DIR, "ConflictPoints.txt")
TB_CONFLICT_XLSX  <- file.path(TB_POINTS_DIR, "conflict_table_df.xlsx")

## TIF-based predictor source directories (pre-converted from RData, read-only inputs)
TB_PREDICTORS_TIF_DIR     <- file.path(TB_DATA_ROOT, "Predictors_TIF")
TB_PREDICTORS_TIF_PRESENT <- file.path(TB_PREDICTORS_TIF_DIR, "present")
TB_PREDICTORS_TIF_FUTURE  <- file.path(TB_PREDICTORS_TIF_DIR, "future")

TB_PA_GDB         <- file.path(TB_DATA_ROOT, "ProtectedAreas", "PAs", "PAs.gpkg")
TB_ROADS_SHP      <- file.path(TB_DATA_ROOT, "Roads", "gis_osm_roads_free_1.shp")

## ---- ENMTML inputs (script 02/03/04/05 build these) ------------------------
TB_PRED_ENMTML_DIR        <- file.path(TB_DATA_ROOT, "predictors_enmtml")
TB_PRED_ENMTML_PRESENT    <- file.path(TB_PRED_ENMTML_DIR, "present")
TB_PRED_ENMTML_FUTURE     <- file.path(TB_PRED_ENMTML_DIR, "future")

TB_OCC_ENMTML             <- file.path(TB_POINTS_DIR, "occ_enmtml.txt")
TB_OCC_CONFLICT_ENMTML    <- file.path(TB_POINTS_DIR, "occ_conflict_enmtml.txt")
TB_TR_MASK_SHP            <- file.path(TB_DATA_ROOT, "TR_mask", "TR_mask.shp")

## ---- Output paths ----------------------------------------------------------
TB_OUT_LOGS       <- file.path(TB_OUT_ROOT, "logs")
TB_OUT_FIGURES    <- file.path(TB_OUT_ROOT, "figures")
TB_OUT_TABLES     <- file.path(TB_OUT_ROOT, "tables")
TB_OUT_RDS        <- file.path(TB_OUT_ROOT, "rds_files")
TB_OUT_RASTERS    <- file.path(TB_OUT_ROOT, "rasters")
TB_OUT_VECTORS    <- file.path(TB_OUT_ROOT, "vectors")

## ENMTML's own result directory (Algorithm/, Ensemble/, Projection/, *.txt)
TB_OUT_ENMTML            <- file.path(TB_OUT_ROOT, "enmtml_result")
TB_OUT_ENMTML_CONFLICT   <- file.path(TB_OUT_ROOT, "enmtml_conflict_result")

## Post-processing derived outputs
TB_OUT_DERIVED       <- file.path(TB_OUT_ROOT, "derived")
TB_OUT_HS_PRESENT    <- file.path(TB_OUT_DERIVED, "hs_present")
TB_OUT_HS_FUTURE     <- file.path(TB_OUT_DERIVED, "hs_future")        # 18 GCM rasters
TB_OUT_HS_AVG        <- file.path(TB_OUT_DERIVED, "hs_gcm_avg")       # 6 final scenarios
TB_OUT_HS_BINARY     <- file.path(TB_OUT_DERIVED, "hs_binary")
TB_OUT_GAINLOSS      <- file.path(TB_OUT_DERIVED, "gain_loss")
TB_OUT_RESISTANCE    <- file.path(TB_OUT_DERIVED, "resistance")
TB_OUT_UNICOR_DIR    <- file.path(TB_OUT_DERIVED, "unicor")
TB_OUT_CONEFOR_DIR   <- file.path(TB_OUT_DERIVED, "conefor")

.tb_dirs <- c(
  TB_OUT_ROOT, TB_OUT_LOGS, TB_OUT_FIGURES, TB_OUT_TABLES, TB_OUT_RDS,
  TB_OUT_RASTERS, TB_OUT_VECTORS,
  TB_OUT_ENMTML, TB_OUT_ENMTML_CONFLICT,
  TB_OUT_DERIVED, TB_OUT_HS_PRESENT, TB_OUT_HS_FUTURE, TB_OUT_HS_AVG,
  TB_OUT_HS_BINARY, TB_OUT_GAINLOSS, TB_OUT_RESISTANCE,
  TB_OUT_UNICOR_DIR, TB_OUT_CONEFOR_DIR,
  TB_PRED_ENMTML_DIR, TB_PRED_ENMTML_PRESENT, TB_PRED_ENMTML_FUTURE
)
invisible(lapply(.tb_dirs, function(d) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}))

## ---- CRS and resolution ----------------------------------------------------
TB_CRS_PROJ <- "+proj=aea +lat_1=37 +lat_2=41 +lat_0=39 +lon_0=35 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
TB_CRS_WGS  <- "EPSG:4326"
TB_RES_M    <- 1000

## ---- Climate scenarios ----------------------------------------------------
TB_GCMS     <- c("GFDL_ESM4", "IPSL_CM6A_LR", "MPI_ESM1_2_HR")
TB_SSPS     <- c("ssp126", "ssp370", "ssp585")
TB_PERIODS  <- c("2041_2070", "2071_2100")

TB_PERIOD_LABELS <- c(
  "2041_2070" = "Near future (2070s)",
  "2071_2100" = "Far future (2100s)"
)
TB_SSP_LABELS <- c(
  "ssp126" = "SSP126",
  "ssp370" = "SSP370",
  "ssp585" = "SSP585"
)
## Pretty GCM names (use "-" instead of "_") for figure titles
TB_GCM_LABELS <- c(
  "GFDL_ESM4"     = "GFDL-ESM4",
  "IPSL_CM6A_LR"  = "IPSL-CM6A-LR",
  "MPI_ESM1_2_HR" = "MPI-ESM1-2-HR"
)
## Pretty scenario id (used by many downstream figures)
tb_pretty_scen <- function(s) {
  out <- s
  for (g in names(TB_GCM_LABELS)) out <- gsub(g, TB_GCM_LABELS[[g]], out, fixed = TRUE)
  for (p in names(TB_SSP_LABELS)) out <- gsub(p, TB_SSP_LABELS[[p]], out, fixed = TRUE)
  out <- gsub("2041_2070", "2070s",  out, fixed = TRUE)
  out <- gsub("2071_2100", "2100s",  out, fixed = TRUE)
  out <- gsub("_", " ", out, fixed = TRUE)
  out
}

## ---- Predictor name groups -------------------------------------------------
TB_NAMES_TOP <- c("Elevation","Slope","Aspect","Roughness")
TB_NAMES_BIO <- sprintf("Bio%02d", 1:19)
TB_NAMES_HUM <- c("d2Forest","d2Agricultural","d2ArtificialSurfaces",
                  "d2Water","d2Roads","PopDen","GHMTS")
TB_NAMES_ALL <- c(TB_NAMES_TOP, TB_NAMES_BIO, TB_NAMES_HUM)

## Future-varying predictors (climate only). Others held constant.
TB_NAMES_FUTURE_VARY <- TB_NAMES_BIO

## ---- ENMTML configuration constants ----------------------------------------
TB_ENM_ALGORITHMS   <- c("BIO","GLM","GAM","SVM","RDF","BRT","MXD","MAH")
## NOTE: GEO_CONST → ENV_CONST. ENMTML's inv_geo() always calls
## dismo::circles(lonlat=TRUE), which corrupts buffer math when data is in
## projected AEA-meters. ENV_CONST uses inv_bio() → dismo::bioclim, which is
## CRS-agnostic. Methodologically valid: pseudo-absences are placed in
## environmentally dissimilar cells from presences.
TB_ENM_PA_METHOD    <- c(method = "ENV_CONST")
TB_ENM_PA_RATIO     <- 1
TB_ENM_COLIN_VAR    <- c(method = "VIF")
TB_ENM_THIN_OCC     <- c(method = "CELLSIZE")
## NOTE 2026-05-23: BLOCK → BOOT switch. BLOCK + ENV_CONST combination triggered
## "variable lengths differ" in coef() during GLM/GAM parallel fit (likely uneven
## PA distribution across spatial folds). BOOT uses simpler non-spatial CV +
## OptimRandomPoints PA — robust and well-tested. Trade-off: it relaxes spatial
## cross-validation rigour in exchange for a reliable, complete pipeline.
TB_ENM_PART         <- c(method = "BOOT", replicates = "10", proportion = "0.7")
TB_ENM_THR          <- c(type = "MAX_TSS")
## NOTE 2026-05-25: ENMTML expects ensemble as a NAMED VECTOR, not list-of-vectors.
## Previously used list() form caused silent skip of Ensemble step
## (grep("method", names(list())) returns integer(0)).
## Correct syntax: c(method = c("W_MEAN","MEAN"), metric = "TSS")
TB_ENM_ENSEMBLE     <- c(method = c("W_MEAN", "MEAN"), metric = "TSS")
TB_ENM_MSDM         <- NULL
TB_ENM_EXTRAPOLATION<- TRUE
TB_ENM_MIN_OCC      <- 10
TB_ENM_IMP_VAR      <- TRUE
TB_ENM_SAVE_PART    <- TRUE
TB_ENM_SAVE_FINAL   <- TRUE
TB_ENM_CORES        <- if (TB_ENV == "cluster") 16L else 4L

## ---- Downstream parameters -------------------------------------------------
TB_PATCH_MIN_KM2 <- 83              # source patch threshold (user-set 2026-05-25)

## Resistance: Trainor et al. (2013) negative-exponential form, parameterised
## per Shokri et al. (2021):
##
##   R(h) = 100 − 99 · ((1 − exp(−c·h)) / (1 − exp(−c)))
##
## where h ∈ [0, 1] is habitat suitability and c is a shape constant.
## c = 4 reflects bears' ability to move through less-suitable habitat
## (smaller c = more linear; larger c = stronger drop-off as h increases).
##   h = 1 → R = 1   (best habitat, lowest resistance)
##   h = 0 → R = 100 (worst habitat, maximum barrier)
TB_RESIST_C      <- 4
TB_RESIST_MIN    <- 1
TB_RESIST_MAX    <- 100

cat(sprintf("[paths] TB_ENV = %s\n", TB_ENV))
cat(sprintf("[paths] DATA  = %s\n", TB_DATA_ROOT))
cat(sprintf("[paths] OUT   = %s\n", TB_OUT_ROOT))
cat(sprintf("[paths] ENMTML = %s\n", TB_OUT_ENMTML))
