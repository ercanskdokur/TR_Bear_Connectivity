## ============================================================================
## 20_conflict_enmtml.R
## Project: TR_Bear_Connectivity
## Purpose: Run ENMTML on the 478 bear–human conflict points using the SAME
##   predictor stack and configuration as the bear-occurrence SDM (06_enmtml_run.R).
##   This produces a "conflict suitability" surface comparable to the bear-HS
##   surface; the difference between them is the conflict-mismatch signal
##   used in 21_conflict_overlay.
##
## Configuration mirrors 06_enmtml_run.R exactly:
##   8 algorithms × ENV_CONST PA × VIF colinearity × BOOT(rep=10, p=0.7) ×
##   MAX_TSS threshold × W_MEAN + MEAN ensemble × 18 future projections.
##
## Inputs:
##   data/predictors_enmtml/present/*.tif        (30 layers, same as bear SDM)
##   data/predictors_enmtml/future/<scenario>/*.tif × 18
##   data/points/occ_conflict_enmtml.txt         (Ursus_arctos_conflict, 478)
##   data/TR_mask/TR_mask.shp                    (accessible area)
## Outputs (ENMTML auto-builds under TB_OUT_ENMTML_CONFLICT):
##   Algorithm/, Ensemble/, Projection/, Evaluation_Table.txt, …
## ============================================================================

suppressPackageStartupMessages({
  library(raster); library(ENMTML); library(sf); library(terra)
})

## ============================================================================
## RUNTIME PATCH (same as 06): ENMTML::FitENM_TMLA_Parallel write fix
## ============================================================================
local({
  ns <- asNamespace("ENMTML")
  orig <- get("FitENM_TMLA_Parallel", envir = ns)
  bd_text <- deparse(body(orig), width.cutoff = 500L)
  hit <- grep("lapply\\(InfoModeling,\\s*write\\b", bd_text)
  if (length(hit) == 0) {
    warning("PATCH FAIL: could not locate 'lapply(InfoModeling, write,' in FitENM_TMLA_Parallel body")
  } else {
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
tb_log_init("20_conflict_enmtml")
tb_pkg_versions(c("ENMTML","terra","sf","dismo","randomForest","gbm","maxnet","kernlab"))

## ============================================================================
## 1. SANITY CHECKS
## ============================================================================
tb_log_section("1. SANITY CHECKS")

stopifnot(dir.exists(TB_PRED_ENMTML_PRESENT))
stopifnot(dir.exists(TB_PRED_ENMTML_FUTURE))
stopifnot(file.exists(TB_OCC_CONFLICT_ENMTML))
stopifnot(file.exists(TB_TR_MASK_SHP))

n_present  <- length(list.files(TB_PRED_ENMTML_PRESENT, pattern = "\\.tif$"))
n_scenarios<- length(list.dirs(TB_PRED_ENMTML_FUTURE, recursive = FALSE))
occ_df     <- read.table(TB_OCC_CONFLICT_ENMTML, header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE)

tb_log(sprintf("present predictors  : %d TIFs in %s", n_present, TB_PRED_ENMTML_PRESENT))
tb_log(sprintf("future scenarios    : %d sub-folders in %s", n_scenarios, TB_PRED_ENMTML_FUTURE))
tb_log(sprintf("conflict records    : %d rows in %s", nrow(occ_df), TB_OCC_CONFLICT_ENMTML))
tb_log(sprintf("accessible-area shp : %s", TB_TR_MASK_SHP))
tb_log(sprintf("result directory    : %s", TB_OUT_ENMTML_CONFLICT))

## ============================================================================
## 2. CLEAN UP RESULT DIR
## ============================================================================
tb_log_section("2. PREPARE RESULT DIR")
if (dir.exists(TB_OUT_ENMTML_CONFLICT)) {
  tb_log(sprintf("clearing existing result dir: %s", TB_OUT_ENMTML_CONFLICT))
  unlink(TB_OUT_ENMTML_CONFLICT, recursive = TRUE, force = TRUE)
}
dir.create(TB_OUT_ENMTML_CONFLICT, recursive = TRUE)

## ============================================================================
## 2b. REPAIR ACCESSIBLE AREA GEOMETRY (re-use validated mask if present)
## ============================================================================
tb_log_section("2b. ACCESSIBLE AREA MASK")

mask_fixed_dir <- file.path(TB_OUT_ROOT, "tmp_mask")
mask_fixed_shp <- file.path(mask_fixed_dir, "TR_mask_valid.shp")

if (!file.exists(mask_fixed_shp)) {
  mask_orig <- sf::st_read(TB_TR_MASK_SHP, quiet = TRUE)
  tb_log(sprintf("mask original: features=%d  valid=%s",
                 nrow(mask_orig), all(sf::st_is_valid(mask_orig))))
  mask_fixed <- sf::st_make_valid(mask_orig)
  mask_fixed <- sf::st_collection_extract(mask_fixed, "POLYGON")
  mask_fixed <- sf::st_union(mask_fixed)
  mask_fixed <- sf::st_make_valid(mask_fixed)
  mask_fixed <- sf::st_sf(geometry = mask_fixed)
  dir.create(mask_fixed_dir, recursive = TRUE, showWarnings = FALSE)
  sf::st_write(mask_fixed, mask_fixed_shp, delete_dsn = TRUE, quiet = TRUE)
  tb_log(sprintf("mask repaired and written -> %s", mask_fixed_shp))
} else {
  tb_log(sprintf("reusing existing validated mask: %s", mask_fixed_shp))
}

## ============================================================================
## 3. CALL ENMTML
## ============================================================================
tb_log_section("3. ENMTML CALL")
tb_log(sprintf("algorithms     : %s", paste(TB_ENM_ALGORITHMS, collapse = ", ")))
tb_log(sprintf("PA method      : %s", paste(TB_ENM_PA_METHOD,  collapse = "=")))
tb_log(sprintf("colinearity    : %s", paste(TB_ENM_COLIN_VAR,  collapse = "=")))
tb_log(sprintf("partition      : %s", paste(TB_ENM_PART,       collapse = "=")))
tb_log(sprintf("threshold      : %s", paste(TB_ENM_THR,        collapse = "=")))
tb_log(sprintf("PA ratio       : %s", TB_ENM_PA_RATIO))
tb_log(sprintf("cores          : %s", TB_ENM_CORES))
tb_log(sprintf("extrapolation  : %s", TB_ENM_EXTRAPOLATION))

tb_tic()
ENMTML::ENMTML(
  pred_dir            = TB_PRED_ENMTML_PRESENT,
  proj_dir            = TB_PRED_ENMTML_FUTURE,
  result_dir          = TB_OUT_ENMTML_CONFLICT,
  occ_file            = TB_OCC_CONFLICT_ENMTML,
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
tb_toc("ENMTML conflict total")

## ============================================================================
## 4. POST-CHECK
## ============================================================================
tb_log_section("4. POST-CHECK")

list_outputs <- function(dir, label) {
  if (!dir.exists(dir)) { tb_log(sprintf("MISSING: %s (%s)", label, dir), "WARN"); return() }
  f <- list.files(dir, recursive = TRUE)
  tb_log(sprintf("%-30s n_files=%d", label, length(f)))
}

list_outputs(file.path(TB_OUT_ENMTML_CONFLICT, "Algorithm"),  "Algorithm/")
list_outputs(file.path(TB_OUT_ENMTML_CONFLICT, "Ensemble"),   "Ensemble/")
list_outputs(file.path(TB_OUT_ENMTML_CONFLICT, "Projection"), "Projection/")

eval_table <- file.path(TB_OUT_ENMTML_CONFLICT, "Evaluation_Table.txt")
if (file.exists(eval_table)) {
  ev <- read.table(eval_table, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  tb_log(sprintf("Evaluation_Table.txt rows=%d cols=%d", nrow(ev), ncol(ev)))
  tb_save_table(ev, "20_conflict_enmtml_eval_raw")
} else {
  tb_log("Evaluation_Table.txt missing — ENMTML failed?", "WARN")
}

tb_log_session()
tb_log("20_conflict_enmtml DONE")
