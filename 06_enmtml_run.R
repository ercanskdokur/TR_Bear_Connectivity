## ============================================================================
## 06_enmtml_run.R
## Project: TR_Bear_Connectivity — ENMTML pipeline
## Purpose: SINGLE master ENMTML() call covering the entire SDM workflow:
##   - 8 algorithms (BIO, GLM, GAM, SVM, RDF, BRT, MXD, MAH)
##   - VIF collinearity reduction
##   - ENV_CONST pseudo-absence (environmental dissimilarity), balanced 1:1
##     (ENMTML's GEO_CONST assumes lon/lat and breaks on projected AEA-m data)
##   - BOOT cross-validation (bootstrap, 10 replicates, 70/30 split)
##   - TR_mask.shp accessible area
##   - 8 thresholds-driven binary maps (MAX_TSS)
##   - W_MEAN (TSS-weighted) + MEAN ensembles
##   - 18 future projections (2 periods × 3 SSPs × 3 GCMs)
##   - MOP extrapolation analysis
## Inputs:
##   data/predictors_enmtml/present/*.tif        (30 layers)
##   data/predictors_enmtml/future/<scenario>/*.tif × 18
##   data/points/occ_enmtml.txt                  (species, x, y)
##   data/TR_mask/TR_mask.shp
## Outputs (ENMTML auto-builds under result_dir):
##   outputs/enmtml_result/
##     Algorithm/{BIO,GLM,GAM,SVM,RDF,BRT,MXD,MAH}/
##     Ensemble/{W_MEAN,MEAN}/
##     Projection/<scenario>/{Algorithm,Ensemble}/
##     Evaluation_Table.txt
##     Thresholds_Algorithm.txt   Thresholds_Ensemble.txt
##     InfoModelling.txt
##     Occurrences_Cleaned.txt    Occurrences_Filtered.txt
##     Moran_*  Mess_*
##     Extent_Masks/, BLOCK/, Pseudoabsence/
## ============================================================================

## raster must load BEFORE terra so terra's functions take precedence
suppressPackageStartupMessages({
  library(raster); library(ENMTML); library(sf); library(terra)
})

## ============================================================================
## RUNTIME PATCH for ENMTML::FitENM_TMLA_Parallel (2026-05-24)
## Bug: lapply(InfoModeling, write, ...) crashes with
##   "argument 1 (type 'list') cannot be handled by 'cat'"
## because one element of InfoModeling list is itself a list/complex object
## that base::write -> base::cat cannot handle. This happens AFTER
## "Models fitted!" but blocks pipeline from proceeding to replicate 2 / ensemble.
## Fix: replace the lapply call to coerce each element via as.character(unlist())
## before writing.
## ============================================================================
local({
  ns <- asNamespace("ENMTML")
  orig <- get("FitENM_TMLA_Parallel", envir = ns)
  bd_text <- deparse(body(orig), width.cutoff = 500L)
  hit <- grep("lapply\\(InfoModeling,\\s*write\\b", bd_text)
  if (length(hit) == 0) {
    warning("PATCH FAIL: could not locate 'lapply(InfoModeling, write,' in FitENM_TMLA_Parallel body")
  } else {
    ## Original:  lapply(InfoModeling, write, <path>, append=T, ncolumns=20, sep="\t")
    ## Patched:   lapply(InfoModeling, function(.elem) base::write(as.character(unlist(.elem)), <path>, append=T, ncolumns=20, sep="\t"))
    ## Net +1 unmatched open paren; append ) at end of line.
    bd_text[hit] <- sub(
      "lapply\\(InfoModeling,\\s*write,\\s*",
      "lapply(InfoModeling, function(.elem) base::write(as.character(unlist(.elem)), ",
      bd_text[hit]
    )
    bd_text[hit] <- paste0(bd_text[hit], ")")
    new_body <- parse(text = paste(bd_text, collapse = "\n"))[[1]]
    body(orig) <- new_body
    assignInNamespace("FitENM_TMLA_Parallel", orig, ns = "ENMTML")
    cat(sprintf("[patch] ENMTML::FitENM_TMLA_Parallel InfoModeling-write patched (line %d)\n", hit[1]))
  }
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
tb_log_init("06_enmtml_run")
tb_pkg_versions(c("ENMTML","terra","sf","dismo","randomForest","gbm","maxnet","kernlab"))

## ============================================================================
## 1. SANITY CHECKS
## ============================================================================
tb_log_section("1. SANITY CHECKS")

stopifnot(dir.exists(TB_PRED_ENMTML_PRESENT))
stopifnot(dir.exists(TB_PRED_ENMTML_FUTURE))
stopifnot(file.exists(TB_OCC_ENMTML))
stopifnot(file.exists(TB_TR_MASK_SHP))

n_present  <- length(list.files(TB_PRED_ENMTML_PRESENT, pattern = "\\.tif$"))
n_scenarios<- length(list.dirs(TB_PRED_ENMTML_FUTURE, recursive = FALSE))
occ_df     <- read.table(TB_OCC_ENMTML, header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE)

tb_log(sprintf("present predictors  : %d TIFs in %s", n_present, TB_PRED_ENMTML_PRESENT))
tb_log(sprintf("future scenarios    : %d sub-folders in %s", n_scenarios, TB_PRED_ENMTML_FUTURE))
tb_log(sprintf("occurrence records  : %d rows in %s", nrow(occ_df), TB_OCC_ENMTML))
tb_log(sprintf("accessible-area shp : %s", TB_TR_MASK_SHP))
tb_log(sprintf("result directory    : %s", TB_OUT_ENMTML))

if (n_present < 30)   tb_log("WARN: <30 present TIFs (expected 30)", "WARN")
if (n_scenarios < 18) tb_log("WARN: <18 future sub-folders (expected 18)", "WARN")

## verify each scenario sub-folder has the same TIFs as present
present_names <- sort(tools::file_path_sans_ext(
  list.files(TB_PRED_ENMTML_PRESENT, pattern = "\\.tif$")))
scenarios <- list.dirs(TB_PRED_ENMTML_FUTURE, recursive = FALSE)
for (s in scenarios) {
  scn_names <- sort(tools::file_path_sans_ext(
    list.files(s, pattern = "\\.tif$")))
  if (!identical(scn_names, present_names))
    tb_log(sprintf("WARN: scenario %s predictors differ from present", basename(s)), "WARN")
}

## ============================================================================
## 2. CLEAN UP RESULT DIR (ENMTML refuses to overwrite some outputs)
## ============================================================================
tb_log_section("2. PREPARE RESULT DIR")
if (dir.exists(TB_OUT_ENMTML)) {
  tb_log(sprintf("clearing existing result dir: %s", TB_OUT_ENMTML))
  unlink(TB_OUT_ENMTML, recursive = TRUE, force = TRUE)
}
dir.create(TB_OUT_ENMTML, recursive = TRUE)

## ============================================================================
## 2b. REPAIR ACCESSIBLE AREA GEOMETRY
## ============================================================================
tb_log_section("2b. REPAIR ACCESSIBLE AREA MASK")

mask_orig <- sf::st_read(TB_TR_MASK_SHP, quiet = TRUE)
tb_log(sprintf("mask original: features=%d  CRS=%s  valid=%s",
               nrow(mask_orig),
               sf::st_crs(mask_orig)$input %||% "unknown",
               all(sf::st_is_valid(mask_orig))))

## Fix invalid geometries, extract only polygons (drops lines/points from
## GeometryCollections), then dissolve to one geometry
mask_fixed <- sf::st_make_valid(mask_orig)
mask_fixed <- sf::st_collection_extract(mask_fixed, "POLYGON")
mask_fixed <- sf::st_union(mask_fixed)            # sfc of one MULTIPOLYGON
mask_fixed <- sf::st_make_valid(mask_fixed)        # second pass after union
mask_fixed <- sf::st_sf(geometry = mask_fixed)     # proper sf object

tb_log(sprintf("mask repaired: features=%d  CRS=%s  valid=%s",
               nrow(mask_fixed),
               sf::st_crs(mask_fixed)$input %||% "unknown",
               all(sf::st_is_valid(mask_fixed))))

mask_fixed_dir <- file.path(TB_OUT_ROOT, "tmp_mask")
dir.create(mask_fixed_dir, recursive = TRUE, showWarnings = FALSE)
mask_fixed_shp <- file.path(mask_fixed_dir, "TR_mask_valid.shp")
sf::st_write(mask_fixed, mask_fixed_shp, delete_dsn = TRUE, quiet = TRUE)
tb_log(sprintf("mask written -> %s", mask_fixed_shp))

## ============================================================================
## 3. CALL ENMTML
## ============================================================================
tb_log_section("3. ENMTML CALL")
tb_log(sprintf("algorithms     : %s", paste(TB_ENM_ALGORITHMS, collapse = ", ")))
tb_log(sprintf("PA method      : %s", paste(TB_ENM_PA_METHOD,  collapse = "=")))
tb_log(sprintf("colinearity    : %s", paste(TB_ENM_COLIN_VAR,  collapse = "=")))
tb_log(sprintf("partition      : %s", paste(TB_ENM_PART,       collapse = "=")))
tb_log(sprintf("threshold      : %s", paste(TB_ENM_THR,        collapse = "=")))
tb_log(sprintf("ensemble #1    : %s", paste(TB_ENM_ENSEMBLE[[1]], collapse = "=")))
tb_log(sprintf("ensemble #2    : %s", paste(TB_ENM_ENSEMBLE[[2]], collapse = "=")))
tb_log(sprintf("PA ratio       : %s", TB_ENM_PA_RATIO))
tb_log(sprintf("cores          : %s", TB_ENM_CORES))
tb_log(sprintf("extrapolation  : %s", TB_ENM_EXTRAPOLATION))

tb_tic()
ENMTML::ENMTML(
  pred_dir            = TB_PRED_ENMTML_PRESENT,
  proj_dir            = TB_PRED_ENMTML_FUTURE,
  result_dir          = TB_OUT_ENMTML,
  occ_file            = TB_OCC_ENMTML,
  sp                  = "species",
  x                   = "x",
  y                   = "y",

  min_occ             = TB_ENM_MIN_OCC,
  thin_occ            = TB_ENM_THIN_OCC,
  eval_occ            = NULL,

  colin_var           = TB_ENM_COLIN_VAR,
  imp_var             = TB_ENM_IMP_VAR,

  sp_accessible_area  = c(method = "MASK", filepath = mask_fixed_shp),
  pseudoabs_method    = TB_ENM_PA_METHOD,
  pres_abs_ratio      = TB_ENM_PA_RATIO,

  part                = TB_ENM_PART,
  save_part           = TB_ENM_SAVE_PART,
  save_final          = TB_ENM_SAVE_FINAL,

  algorithm           = TB_ENM_ALGORITHMS,
  thr                 = TB_ENM_THR,
  msdm                = TB_ENM_MSDM,
  ensemble            = TB_ENM_ENSEMBLE,
  extrapolation       = TB_ENM_EXTRAPOLATION,
  cores               = TB_ENM_CORES
)
tb_toc("ENMTML total")

## ============================================================================
## 4. POST-CHECK + INDEX OUTPUTS
## ============================================================================
tb_log_section("4. POST-CHECK")

list_outputs <- function(dir, label) {
  if (!dir.exists(dir)) { tb_log(sprintf("MISSING: %s (%s)", label, dir), "WARN"); return() }
  f <- list.files(dir, recursive = TRUE)
  tb_log(sprintf("%-30s n_files=%d", label, length(f)))
}

list_outputs(file.path(TB_OUT_ENMTML, "Algorithm"),  "Algorithm/")
list_outputs(file.path(TB_OUT_ENMTML, "Ensemble"),   "Ensemble/")
list_outputs(file.path(TB_OUT_ENMTML, "Projection"), "Projection/")

eval_table <- file.path(TB_OUT_ENMTML, "Evaluation_Table.txt")
if (file.exists(eval_table)) {
  ev <- read.table(eval_table, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  tb_log(sprintf("Evaluation_Table.txt rows=%d cols=%d", nrow(ev), ncol(ev)))
  tb_log(sprintf("metrics: %s", paste(setdiff(names(ev), c("Sp","Alg","Part","Thr")),
                                      collapse = ", ")))
  tb_save_table(ev, "06_enmtml_eval_raw")
} else {
  tb_log("Evaluation_Table.txt missing — ENMTML failed?", "WARN")
}

tb_log_session()
tb_log("06_enmtml_run DONE")
