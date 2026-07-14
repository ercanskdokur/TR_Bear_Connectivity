# Safeguarding the future ecological network of brown bears across Anatolia: a connectivity-based conservation priority

This repository contains the complete analysis code for the brown bear
species-distribution-model (SDM) + landscape-connectivity study for Türkiye.
The pipeline links an 8-algorithm ENMTML ensemble SDM to a resistance surface,
derives least-cost corridors (UNICOR) and graph-theoretic connectivity indices
(Conefor: PC, IIC, dPC, dIIC), projects all of these onto present and future
(CMIP6) climate, and prioritises core habitats, corridors, pinch points and
protected-area gaps.

All scripts are written in R (analysis) and Bash/SLURM (HPC job wrappers).
Comments and outputs are in English; only proper nouns (e.g. "Türkiye",
Turkish province names) retain their native spelling.

---

## 1. How to read this archive

* Files are numbered in **execution order**. Run them by ascending prefix:
  `00_*` (setup) → `01`…`38`, then the supporting figure scripts
  (`S22_paired.R`, `fig6_costdist.R`).
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

The exact computational environment is pinned by the container definitions, so
the analysis can be reproduced bit-for-bit:

| File | Role |
|---|---|
| `00_Dockerfile_trbear` | Docker image recipe (R + geospatial stack + ENMTML + UNICOR + Conefor) |
| `00_singularity_trbear.def` | Equivalent Apptainer/Singularity definition used on the HPC |
| `00_build_docker_local.ps1` | Helper to build the image locally (Windows/PowerShell) |
| `00_convert_to_sif.slurm` | Convert the Docker image to a `.sif` Apptainer image on the cluster |

Key external tools embedded in the container: **ENMTML** (ensemble SDM),
**UNICOR** (resistant kernel / least-cost corridors), **Conefor**
(graph connectivity indices), and the R `terra`/`sf`/`landscapemetrics` stack.
The Dockerfile clones/installs these at build time; the exact **built image
(`trbear.tar.gz`) is archived on Zenodo** (see §6) so the environment can be
obtained without rebuilding.

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
| `17_conefor.R` | Euclidean graph indices: dPC, dIIC patch prioritisation |
| `18_pa_overlay.R` | Core/corridor × protected-area overlap |
| `19_roads_overlay.R` | Road crossings + corridor × road overlay |

### Phase 6 — Human–bear conflict
| Script | Purpose |
|---|---|
| `20_conflict_enmtml.R` | Conflict risk modeling run (`ConflictPoints.txt`) |
| `21_conflict_overlay.R` | Conflict typology × bioregion × connectivity |
| `25_validate_conflict_sdm.R` | Test conflict-vs-bear SDM redundancy (circular-predictor check) |
| `26_compound_risk.R` | Compound-risk hotspots (corridor ∩ roads ∩ high conflict) |

### Phase 7 — Summary tables
| Script | Purpose |
|---|---|
| `23_tables_master.R` | Consolidated summary tables from the pipeline outputs |

### Phase 8 — Connectivity-manuscript analyses
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

### Phase 9 — Cost-distance connectivity analyses (PRIMARY in final)
| Script | Purpose |
|---|---|
| `35_costdist_conefor.R` | Recompute PC/IIC/dPC/dIIC on **least-cost (effective) distances** |
| `36_c_sensitivity.R` | Sensitivity to resistance shape constant *c* |
| `37_patch_threshold.R` | Sensitivity to source-patch minimum-size threshold (50/83/120 km²) |
| `38_conflict_sdm_figs.R` | Full diagnostic figures/tables for the conflict ENM |
| `S22_paired.R` | Fig. S22 — paired dPC vs. dIIC patch importance |
| `fig6_costdist.R` | Regenerate Fig. 6 with cost-distance PC in panel C |

## 5. Notes on reproducibility

* Scripts are deterministic given the inputs, except where random permutations
  are used (e.g. random node-removal envelopes in `29`); those set seeds
  internally.
* Several `06z*` / `24` scripts are documented **workarounds** for known bugs in
  ENMTML 1.0.0 (ensemble assembly, threshold/BOOT paths). They are part of the
  reproducible pipeline and must be run in the order shown.
* The cost-distance graph (scripts `35`–`38`) is the **primary** connectivity
  analysis in the final; the Euclidean graph (`17`) is retained for
  comparison with the Euclidean-distance formulation.
* **Environment & container.** The **canonical, exactly-versioned environment is
  the archived Docker image (`trbear.tar.gz`) on Zenodo** — pull it to reproduce the analysis
  bit-for-bit. The Dockerfile/def rebuild that environment from source: R packages
  install from a dated Posit PPM snapshot (`2023-12-01`; `rgdal`/`rgeos` from
  `2023-10-01`, just before CRAN archival), while ENMTML and UNICOR default to
  their upstream `main`. For a from-scratch rebuild that matches this study
  exactly, pin those two to the commits used — `--build-arg ENMTML_REF=<sha>
  --build-arg UNICOR_REF=<sha>` (Docker) or `export ENMTML_REF` / `UNICOR_REF`
  (Singularity) — otherwise pull the archived image. The Conefor CLI fetches from a fixed
  URL and the build fails loudly if unavailable; Python packages are unpinned and
  captured exactly only in the archived image.

## 6. Availability, DOIs and how to cite

This study is distributed as two linked deposits:

| Component | Location | DOI | License |
|---|---|---|---|
| **Analysis code** (this repository) | GitHub | `<Zenodo software DOI>` | MIT (see `LICENSE`) |
| **Input data + Docker container image** (`trbear.tar.gz`) | Zenodo | `<Zenodo data DOI>` | CC-BY-4.0 |

* The code lives on GitHub; a tagged **GitHub Release** is mirrored to **Zenodo**,
  which mints a permanent, citable **software DOI** for that exact snapshot
  (enable the repository in your Zenodo account, then publish a release).
* The **input data** (rasters, occurrence/conflict points, protected-area and
  road layers) and the **Docker image `trbear.tar.gz`** are deposited on **Zenodo** under
  **CC-BY-4.0**, so any reuse must cite this work.
* To reproduce: obtain the data deposit, unzip it into a folder named `data/`
  beside the scripts (or point `TB_DATA_ROOT` at it), load `trbear.tar.gz` from
  Zenodo (`docker load`, or build a `.sif` with
  `apptainer build trbear.sif docker-archive://trbear.tar`; or rebuild from
  `00_Dockerfile_trbear`), then run the scripts in numeric order.
* Please cite the published article together with both DOIs. 
