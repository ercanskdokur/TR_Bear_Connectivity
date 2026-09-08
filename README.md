# Severing the Anatolian land bridge: climate change and protection gaps threaten brown bear connectivity

This repository contains the complete analysis code for the brown bear
species-distribution-model (SDM) + landscape-connectivity study for Türkiye.
The pipeline links an ENMTML ensemble SDM to a resistance surface, derives
least-cost corridors (UNICOR) and graph-theoretic connectivity indices
(PC, IIC, dPC, dIIC and the area-comparable ECA, computed in R following the
published Conefor definitions), projects all of these onto present and future
(CMIP6) climate, and prioritises core habitats, corridors, pinch points and
protected-area gaps.

**A note on the ensemble, because two variants exist in the outputs.** The SDM
is fitted with eight algorithms. `06z2_ensemble_threshold_fix.R` derives the
MAX_TSS threshold from, and writes the binary habitat map from, the full
**eight-algorithm** ensemble; that binary is what `13_unicor_prep.R` turns into
the 93 source patches, so every graph result rests on it. `24_ensemble_no_mah.R`
then drops Mahalanobis distance (2.5% of the TSS weight) and overwrites the
*continuous* `Ensemble/W_MEAN/<sp>.tif` in place, backing the original up as
`<sp>_orig8alg.tif`; that **seven-algorithm** surface is what `12_resistance.R`
converts into movement resistance. The two surfaces differ by a mean absolute
0.006 on the 0-1 scale (r = 0.9999), so the split is immaterial, but scripts
must not be re-ordered on the assumption that a single ensemble is in play.
`C1_maxTSS_sensitivity.R` sweeps both variants and reproduces the reported
baseline (93 cores, 141,049 km2) exactly.

All scripts are written in R (analysis) and Bash/SLURM (HPC job wrappers).
Comments and outputs are in English; only proper nouns (e.g. "Türkiye",
Turkish province names) retain their native spelling.

---

## 1. How to read this archive

* Files are numbered in **execution order**. Run them by ascending prefix:
  `00_*` (setup) → `01`…`38`, then the supporting figure script (`S22_paired.R`),
  and finally the robustness strands (`master_submit_all.sh`, §4 Phase 9).
* The numbers `17` and `22`–`23` are intentionally absent. They were exploratory
  steps (a Euclidean-distance connectivity graph and two early figure/table
  assemblers) that the final manuscript does not use; they were removed rather
  than shipped as dead code. Numbering of the remaining steps is unchanged so
  that it still matches the outputs they write.
* Every analysis step exists as a pair: an **`.R`** file (the analysis) and a
  matching **`.slurm`** file (the HPC job that runs it). A few steps are bundled
  in one wrapper (`28_33_extra_analyses.slurm`, `36_37_sensitivity.slurm`).
* Every R script begins with `source("00_paths.R"); source("00_helpers.R")`,
  so those two files must sit in the working directory. No other script sources
  another — steps communicate only through files written to `outputs/`.
* **Paths are relative / environment-driven — no absolute paths are hard-coded.**
  `00_paths.R` resolves everything from the project root (the folder that holds
  the scripts, with `data/` and `outputs/` beside it), so the archive runs from
  wherever it is unpacked. To reproduce the original split-volume HPC layout,
  override any of `TB_ROOT`, `TB_DATA_ROOT`, `TB_OUT_ROOT`, `TB_PROGRAMS` as
  environment variables.
* The `.slurm` files document the exact HPC invocation (SLURM + Apptainer) and
  are **templates**: set `PROJECT_ROOT` (and, if needed, `TB_DATA_ROOT`) in your
  environment, create a `logs/` folder next to the scripts, and adapt the
  `#SBATCH` directives to your own scheduler. The science is fully contained in
  the `.R` files.

## 2. Software environment

**The `Dockerfile` in this repository is the computational environment.** No
prebuilt image is distributed; the environment is rebuilt from this recipe. It
pins the base image (`rocker/geospatial:4.3.2`), installs every R package from a
dated CRAN snapshot (Posit PPM `2023-12-01`; `rgdal`/`rgeos` from `2023-10-01`,
just before their CRAN archival), and exposes the two source-installed tools as
build arguments so they can be pinned to exact commits.

| File | Role |
|---|---|
| `Dockerfile` | **Canonical** environment recipe (R + geospatial stack + ENMTML + UNICOR) |
| `00_singularity_trbear.def` | Equivalent Apptainer/Singularity recipe, for building directly on an HPC |
| `00_build_docker_local.ps1` | Convenience wrapper for building the image on Windows/PowerShell |
| `00_convert_to_sif.slurm` | Convert a locally built Docker image to a `.sif` Apptainer image on a cluster |

Build it:

```bash
docker build -t trbear:latest -f Dockerfile .
```

To reproduce the exact study environment, additionally pin the two upstream
tools to the commits used (this study used **ENMTML 1.0.0**):

```bash
docker build -t trbear:latest -f Dockerfile . \
    --build-arg ENMTML_REF=<commit-sha> --build-arg UNICOR_REF=<commit-sha>
```

To run on an HPC with Apptainer/Singularity, either build from
`00_singularity_trbear.def` directly, or export the Docker image and convert it:

```bash
docker save trbear:latest | gzip > trbear.tar.gz    # transfer to the cluster
apptainer build trbear.sif docker-archive://trbear.tar.gz
```

Every analysis script then runs inside the container, e.g.
`apptainer exec trbear.sif Rscript 08_postprocess_present.R`.

Key external tools: **ENMTML** (ensemble SDM), **UNICOR** (resistant kernel /
least-cost corridors) and the R `terra`/`sf`/`landscapemetrics` stack. A Conefor
command-line binary is also fetched when available, but **no analysis script
uses it** — all connectivity indices are computed in R following the published
Conefor definitions; it is installed only for optional independent cross-checks,
and the build tolerates its (unstable) download failing.

## 3. Input data (deposited on Zenodo)

The scripts expect the following raw inputs (paths resolve under
`<project root>/data/`, or `TB_DATA_ROOT`). These are **not** code; they are
archived as a separate **Zenodo** data deposit (see §6) under CC-BY-4.0:

* **Climate** — CHELSA bioclimatic (Bio1-19) variables (present + CMIP6 future).
* **Topography** — SRTM-derived elevation/slope/aspect/roughness.
* **Anthropogenic** — CORINE land cover, OpenStreetMap roads, Global Human Settlement Layer Population Density data, Global Human Modification of Terrestrial Systems (GHMTS).
* **Occurrence** — brown bear presence points (`PresencePoints.txt`).
* **Conflict** — human–bear conflict records (`ConflictPoints.txt`).
* **Masks / boundaries** — Türkiye land mask, protected-area layer.

### Climate / projection configuration

* **CRS:** Türkiye Albers Equal-Area
  (`+proj=aea +lat_1=37 +lat_2=41 +lat_0=39 +lon_0=35`)
* **Resolution:** 1 km
* **GCMs:** GFDL-ESM4, IPSL-CM6A-LR, MPI-ESM1-2-HR
* **SSPs:** ssp126, ssp370, ssp585
* **Periods:** 2041–2070, 2071–2100 (3 GCM × 3 SSP × 2 periods = 18 future runs;
  averaged to 6 GCM-mean scenarios)
* **SDM algorithms:** BIO, GLM, GAM, SVM, RDF, BRT, MXD, MAH (8), with bootstrap
  cross-validation (BOOT, 10 replicates, 70/30 split); MAH later dropped (random
  skill) → 7-algorithm working ensemble.
* **Dispersal distances:** 50, 100, 150, 200, 300, 400 km (focal d = 100 km).

## 4. Pipeline — script by script

### Phase 0 — Setup
| Script | Purpose |
|---|---|
| `00_paths.R` | All path variables + ENMTML parameter constants |
| `00_helpers.R` | Loggers, ggplot theme/palette, spatial thinning, IO helpers |

### Phase 1 — Data preparation
| Script | Purpose |
|---|---|
| `01_explore.R` | Explore predictors + presence + conflict points (summary + maps) |
| `02_predictors_present_to_tif.R` | Present predictors → GeoTIFF, reproject to AEA, apply TR mask |
| `03_predictors_future_to_tif.R` | 18 future climate sets → per-scenario GeoTIFF folders |
| `04_points_prep.R` | Clean occurrence points → ENMTML format (`sp / x / y`) |
| `05_accessible_area.R` | Build Türkiye accessible-area mask (sea excluded) for ENMTML |

### Phase 2 — Species distribution model
| Script | Purpose |
|---|---|
| `06_enmtml_run.R` | Single ENMTML call: 8 algorithms × BOOT CV × 18 future projections |
| `06z_ensemble_manual.R` | Manual ensemble reconstruction (ENMTML 1.0.0 ensemble bug workaround) |
| `06z2_ensemble_threshold_fix.R` | Recompute ensemble thresholds (BOOT/path bug workaround) |
| `24_ensemble_no_mah.R` | Drop MAH (random skill) and recompute W_MEAN/MEAN ensembles |

### Phase 3 — SDM post-processing
| Script | Purpose |
|---|---|
| `07_postprocess_eval.R` | Parse `Evaluation_Table.txt`; per-algorithm performance figure |
| `08_postprocess_present.R` | Present W_MEAN ensemble suitability map + binary |
| `09_postprocess_future_each.R` | 18 future suitability maps → standard format |
| `10_gcm_average.R` | GCM averaging → 6 final scenarios + SD uncertainty maps |
| `11_gain_loss_stable.R` | Present vs. future binary change (gain / loss / stable) |

### Phase 4 — Connectivity (resistance + corridors)
| Script | Purpose |
|---|---|
| `12_resistance.R` | Resistance surface from suitability (present + 6 futures) |
| `13_unicor_prep.R` | ASCII conversion + UNICOR `.rsg` configuration |
| `14_unicor_run.slurm` | UNICOR runs (SLURM array, 7 scenarios) — runs the binary directly |
| `15_unicor_post.R` | UNICOR outputs → corridor-strength + cumulative least-cost maps |

### Phase 5 — Prioritisation (patches, graph, overlays)
| Script | Purpose |
|---|---|
| `16_landscapemetrics.R` | Patch-level metrics (area, ENN, proximity) |
| `18_pa_overlay.R` | Core/corridor × protected-area overlap |
| `19_roads_overlay.R` | Road crossings + corridor × road overlay |

### Phase 6 — Human–bear conflict
| Script | Purpose |
|---|---|
| `20_conflict_enmtml.R` | Conflict risk modeling run (`ConflictPoints.txt`) |
| `21_conflict_overlay.R` | Conflict typology × bioregion × connectivity |
| `25_validate_conflict_sdm.R` | Test conflict-vs-bear SDM redundancy (circular-predictor check) |
| `26_compound_risk.R` | Compound-risk hotspots (corridor ∩ roads ∩ high conflict) |

### Phase 7 — Connectivity-manuscript analyses
| Script | Purpose |
|---|---|
| `27_conefor_components.R` | dPC/dIIC intra/flux/connector decomposition across distances |
| `28_dpc_rankshift.R` | Re-ranking of core importance under climate change |
| `29_network_robustness.R` | Sequential node-removal robustness (targeted/random/area; n50) |
| `30_pa_network_gap.R` | Is the PA system a connected network? Unprotected glue cores |
| `31_corridor_centroid_shift.R` | Northward/upslope shift of the corridor backbone |
| `32_corridor_validation.R` | Independent validation of corridors vs. presence/conflict points |
| `33_pinch_road_priority.R` | Ranked corridor × road pinch points (candidate crossings) |
| `34_core_crosswalk.R` | Assign consistent core labels C01…C93 (single source of truth) |

### Phase 8 — Cost-distance connectivity analyses (PRIMARY in final)
| Script | Purpose |
|---|---|
| `35_costdist_conefor.R` | Recompute PC/IIC/dPC/dIIC on **least-cost (effective) distances** |
| `36_c_sensitivity.R` | Sensitivity to resistance shape constant *c* |
| `37_patch_threshold.R` | Sensitivity to source-patch minimum-size threshold (50/83/120 km²) |
| `38_conflict_sdm_figs.R` | Full diagnostic figures/tables for the conflict ENM |
| `S22_paired.R` | Paired dPC vs. dIIC patch importance (Fig. S16 of the published Supporting Information; the `S22` prefix is a legacy filename from an earlier SI numbering) |

### Phase 9 — Sensitivity / robustness strands
`master_submit_all.sh` submits all of these with correct dependencies;
`check_outputs.sh` verifies the resulting files.

| Script | Purpose |
|---|---|
| `C1_maxTSS_sensitivity.R` | Sensitivity of habitat topology to the MAX_TSS binarisation threshold |
| `C2_pseudoabsence_sensitivity.R` | Sensitivity to the pseudo-absence ratio (1:1, 10:1, 15:1) |
| `C3_c_sensitivity_full.R` | Sensitivity to the resistance shape constant *c*, across present + all future climate scenarios |
| `D_conflict_mc_spatial.R` | Human–bear conflict clustering under three Monte Carlo null models (random, torus-shift, rotation) to rule out spatial-autocorrelation inflation |
| `E_source_independence.R` | Spatial-coincidence check confirming presence and conflict records are independently sourced datasets |
| `F_mess_extrapolation.R` | MESS non-analog (extrapolated) climate area per future scenario |

## 5. Notes on reproducibility

* Scripts are deterministic given the inputs, except where random permutations
  are used (e.g. random node-removal envelopes in `29`); those set seeds
  internally.
* Several `06z*` / `24` scripts are documented **workarounds** for known bugs in
  ENMTML 1.0.0 (ensemble assembly, threshold/BOOT paths). They are part of the
  reproducible pipeline and must be run in the order shown.
* The cost-distance graph (scripts `35`–`37`) is the connectivity analysis
  reported in the final manuscript.
* **Stochastic steps and exact numeric reproduction.** `00_paths.R` defines a
  single seed (`TB_SEED <- 42L`), and every script with a stochastic step
  (ENMTML pseudo-absence sampling and BOOT cross-validation; the Monte Carlo
  permutation tests in `29`, `32`, `D_conflict_mc_spatial.R`) calls
  `set.seed(TB_SEED)` immediately before it. This makes reruns of the
  pure-R Monte Carlo scripts exactly reproducible. ENMTML dispatches its
  pseudo-absence/BOOT-replicate fitting to parallel workers internally, so a
  fixed seed narrows but does not strictly guarantee bit-identical habitat-
  suitability output across reruns or machines. Only the *inputs* are archived
  on Zenodo; the outputs are not deposited, so a rerun reproduces the reported
  results to within that stochastic tolerance rather than bit-for-bit.
* **Sensitivity sweeps vs. the main run.** `C1_maxTSS_sensitivity.R` reads the
  archived ensemble rasters directly rather than refitting, so it reproduces the
  main run's baseline exactly (93 cores, 141,049 km2 at MAX_TSS = 0.662) and its
  absolute numbers are directly comparable with the manuscript.
  `C2_pseudoabsence_sensitivity.R` does refit the ensemble internally with a
  fixed seed, whereas the archived main-pipeline run predates the introduction of
  that seed; C2 is therefore internally comparable across its own three ratios
  but does not reproduce the main run's absolute habitat area at its PA = 1:1
  setting, and is reported as a relative sensitivity only.
* **Environment.** The `Dockerfile` is the canonical environment definition
  (§2). R packages install from a dated Posit PPM snapshot (`2023-12-01`;
  `rgdal`/`rgeos` from `2023-10-01`, just before CRAN archival). ENMTML and
  UNICOR default to their upstream `main`; for a rebuild that matches this study
  exactly, pin both to the commits used — `--build-arg ENMTML_REF=<sha>
  --build-arg UNICOR_REF=<sha>` (Docker) or `export ENMTML_REF` / `UNICOR_REF`
  (Apptainer). Two known limitations: the Python packages installed via `pip` are
  **not** version-pinned, and the optional Conefor CLI is fetched from an unstable
  upstream URL, so its presence is not guaranteed (nothing in the pipeline depends
  on it).

## 6. Availability, DOIs and how to cite

This study is distributed as two linked deposits:

| Component | Location | DOI | License |
|---|---|---|---|
| **Analysis code + environment recipe** (this repository) | GitHub | minted on first release (see below) | MIT (see `LICENSE`) |
| **Input data** | Zenodo | [10.5281/zenodo.21358924](https://doi.org/10.5281/zenodo.21358924) | CC-BY-4.0 |

* The code lives on GitHub; a tagged **GitHub Release** is mirrored to **Zenodo**,
  which mints a permanent, citable **software DOI** for that exact snapshot
  (enable the repository in your Zenodo account, then publish a release).
* The **computational environment** is not distributed as a prebuilt image. It is
  defined by the `Dockerfile` in this repository and is rebuilt from that recipe
  (§2), so it is versioned and citable together with the code.
* The **input data** (rasters, occurrence/conflict points, protected-area and
  road layers) are deposited on **Zenodo**
  ([10.5281/zenodo.21358924](https://doi.org/10.5281/zenodo.21358924)) under
  **CC-BY-4.0**, so any reuse must cite this work.
* To reproduce:
  1. Clone this repository.
  2. Download the Zenodo data deposit and unpack it into a folder named `data/`
     beside the scripts (or point `TB_DATA_ROOT` at it).
  3. Build the container from the `Dockerfile` (§2).
  4. Run the scripts inside the container in the numeric order of §4.

### How to cite

Please cite the article together with both deposits.

**Data**

> Sıkdokur, E., Sağlam, İ. K., Şekercioğlu, Ç. H., & Naderi, M. (2026).
> *Severing the Anatolian land bridge: climate change and protection gaps
> threaten brown bear connectivity* [Data set]. Zenodo.
> https://doi.org/10.5281/zenodo.21358924

**Code**

> Sıkdokur, E., Sağlam, İ. K., Şekercioğlu, Ç. H., & Naderi, M. (2026).
> *TR_Bear_Connectivity: analysis code and environment recipe* [Software].
> Zenodo. DOI assigned when the first tagged release is archived.

**Article**

> Sıkdokur, E., Sağlam, İ. K., Şekercioğlu, Ç. H., & Naderi, M.
> *Severing the Anatolian land bridge: climate change and protection gaps
> threaten brown bear connectivity.* (in review)
