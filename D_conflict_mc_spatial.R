#!/usr/bin/env Rscript
# ============================================================================
# D_conflict_mc_spatial.R
# TR_Bear_Connectivity: Conflict MC Null Model with Spatial Autocorrelation Controls
# Purpose: Test if conflict clustering is real or Monte Carlo artifact (Madde 8)
#
# Perf fix: Step 3's three null loops each called st_distance()
# once per individual point (478 separate GEOS/S2 calls per permutation x
# 999 permutations x 3 nulls), which took well over an hour without finishing
# even the first null. Replaced with a single st_distance() matrix call per
# permutation (all points vs all points at once) followed by a vectorised
# rowwise min -- numerically identical result, ~2 orders of magnitude fewer
# GEOS/S2 invocations.
# ============================================================================

library(tidyverse)
library(raster)
library(sf)
library(spatstat)
select <- dplyr::select  # guard: raster/spatstat mask dplyr::select

source("00_paths.R")

cat("\n=== STRAND D: Conflict MC Spatial Autocorrelation ===\n")

# ---- 1) LOAD DATA ------
cat("Step 1: Loading conflict and occurrence points...\n")

conflict_data <- read_delim(TB_CONFLICT_TXT, delim = "\t", show_col_types = FALSE)
occ_data <- read_delim(TB_PRESENCE_TXT, delim = "\t", show_col_types = FALSE)

conflict_pts <- conflict_data %>% select(x, y) %>% st_as_sf(coords = c("x", "y"), crs = 4326)
occ_pts <- occ_data %>% select(x, y) %>% st_as_sf(coords = c("x", "y"), crs = 4326)

n_conflict <- nrow(conflict_pts)
n_occ <- nrow(occ_pts)

cat(sprintf("✓ Loaded %d conflict points\n", n_conflict))
cat(sprintf("✓ Loaded %d occurrence points\n", n_occ))

# helper: mean of rowwise-nearest-neighbour distance (km) between two sf point sets
mean_nn_dist_km <- function(pts_a, pts_b) {
  D <- st_distance(pts_a, pts_b)                 # n_a x n_b units matrix (m)
  Dn <- matrix(as.numeric(D), nrow = nrow(D)) / 1000
  mean(apply(Dn, 1, min))
}
nn_dist_km <- function(pts_a, pts_b) {
  D <- st_distance(pts_a, pts_b)
  Dn <- matrix(as.numeric(D), nrow = nrow(D)) / 1000
  apply(Dn, 1, min)
}

# ---- 2) OBSERVED CLUSTERING (distance from conflict to nearest occurrence) ------
cat("\nStep 2: Computing observed conflict-occurrence distances...\n")

obs_distances <- nn_dist_km(conflict_pts, occ_pts)

obs_mean_dist <- mean(obs_distances)
obs_median_dist <- median(obs_distances)

cat(sprintf("✓ Mean distance: %.1f km\n", obs_mean_dist))
cat(sprintf("✓ Median distance: %.1f km\n", obs_median_dist))

# Histogram of observed distances
hist_obs <- tibble(distance_km = obs_distances) %>%
  ggplot(aes(x = distance_km)) +
  geom_histogram(bins = 30, fill = "darkred", alpha = 0.6, color = "black") +
  geom_vline(xintercept = obs_mean_dist, linetype = "dashed", color = "red", linewidth = 1.2) +
  labs(title = "Observed: Distance from Conflict to Nearest Occurrence",
       x = "Distance (km)", y = "Frequency",
       subtitle = sprintf("Mean = %.1f km, Median = %.1f km", obs_mean_dist, obs_median_dist)) +
  theme_minimal()

# ---- 3) MONTE CARLO NULLS ------
cat("\nStep 3: Generating null distributions (vectorised)...\n")

n_permutations <- 999
set.seed(42)

# NULL 1: Random resampling (current method; may inflate clustering)
cat("  Null 1: Random resampling (may be inflated)...\n")

extent_box_sfc <- st_as_sfc(st_bbox(occ_pts))
random_null_distances <- vapply(seq_len(n_permutations), function(k) {
  random_pts <- st_sample(extent_box_sfc, size = n_conflict, type = "random")
  mean_nn_dist_km(random_pts, occ_pts)
}, numeric(1))

p_random <- mean(random_null_distances <= obs_mean_dist)
cat(sprintf("    p-value (random null): %.4f\n", p_random))

# NULL 2: Torus shift (preserves spatial autocorrelation structure)
cat("  Null 2: Torus shift (preserves autocorrelation)...\n")

occ_coords <- st_coordinates(occ_pts)
extent <- st_bbox(occ_pts)
width  <- extent["xmax"] - extent["xmin"]
height <- extent["ymax"] - extent["ymin"]

torus_null_distances <- vapply(seq_len(n_permutations), function(k) {
  shift_x <- runif(1, extent["xmin"], extent["xmax"])
  shift_y <- runif(1, extent["ymin"], extent["ymax"])

  shifted_coords <- occ_coords
  shifted_coords[, 1] <- (occ_coords[, 1] - extent["xmin"] + shift_x) %% width + extent["xmin"]
  shifted_coords[, 2] <- (occ_coords[, 2] - extent["ymin"] + shift_y) %% height + extent["ymin"]

  shifted_pts <- st_as_sf(as.data.frame(shifted_coords), coords = c("X", "Y"), crs = 4326)
  mean_nn_dist_km(conflict_pts, shifted_pts)
}, numeric(1))

p_torus <- mean(torus_null_distances <= obs_mean_dist)
cat(sprintf("    p-value (torus shift): %.4f\n", p_torus))

# NULL 3: Random rotation (preserves spatial structure around centroid)
cat("  Null 3: Random rotation (preserves spatial structure)...\n")

centroid <- colMeans(occ_coords)
centered <- t(t(occ_coords) - centroid)

rotation_null_distances <- vapply(seq_len(n_permutations), function(k) {
  angle <- runif(1, 0, 2 * pi)
  rot_matrix <- matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)), nrow = 2)

  rotated <- t(rot_matrix %*% t(centered))
  rotated_coords <- t(t(rotated) + centroid)

  rotated_df <- as.data.frame(rotated_coords)
  names(rotated_df) <- c("X", "Y")
  rotated_pts <- st_as_sf(rotated_df, coords = c("X", "Y"), crs = 4326)
  mean_nn_dist_km(conflict_pts, rotated_pts)
}, numeric(1))

p_rotation <- mean(rotation_null_distances <= obs_mean_dist)
cat(sprintf("    p-value (rotation): %.4f\n", p_rotation))

# ---- 4) COMPARISON ------
cat("\n=== NULL MODEL COMPARISON ===\n")

comparison_df <- tibble(
  Null_Model = c("Random Resampling", "Torus Shift", "Random Rotation"),
  p_value = c(p_random, p_torus, p_rotation),
  Interpretation = c(
    "May be inflated by spatial autocorrelation",
    "Autocorrelation-robust",
    "Autocorrelation-robust"
  )
)

print(comparison_df)

cat("\nInterpretation:\n")
if (p_random < 0.05 & (p_torus > 0.05 | p_rotation > 0.05)) {
  cat("\n⚠ CRITICAL: Random null shows spurious significance.\n")
  cat("  Spatial autocorrelation in conflict pattern inflates clustering signal.\n")
  cat("  RECOMMENDATION: Use torus shift or rotation null as correction.\n")
  cat(sprintf("  Corrected p-value: %.4f (instead of %.4f)\n", max(p_torus, p_rotation), p_random))
} else if (p_random < 0.05 & p_torus < 0.05 & p_rotation < 0.05) {
  cat("\n✓ Conflict clustering is ROBUST:\n")
  cat("  All three nulls show significant p < 0.05.\n")
  cat("  Clustering is real, not an autocorrelation artifact.\n")
} else {
  cat("\n• Clustering is NOT significant (p >= 0.05)\n")
  cat("  Conflicts and occurrences are randomly distributed relative to each other.\n")
}

# ---- 5) VISUALIZATIONS ------
suppressPackageStartupMessages({
  library(ggplot2)
})

# Density plot: all nulls + observed
null_data <- tibble(
  Method = c(
    rep("Observed", length(obs_distances)),
    rep("Random null", length(random_null_distances)),
    rep("Torus shift", length(torus_null_distances)),
    rep("Rotation", length(rotation_null_distances))
  ),
  Mean_Distance = c(obs_distances, random_null_distances, torus_null_distances, rotation_null_distances)
)

p_density <- null_data %>%
  ggplot(aes(x = Mean_Distance, fill = Method, alpha = Method)) +
  geom_density() +
  geom_vline(xintercept = obs_mean_dist, linetype = "dashed", color = "darkred", linewidth = 1.2) +
  scale_alpha_manual(values = c("Observed" = 1.0, "Random null" = 0.4, "Torus shift" = 0.4, "Rotation" = 0.4)) +
  labs(title = "Null Model Comparison: Distance Distributions",
       x = "Mean Distance (km)", y = "Density",
       fill = "Null Model") +
  theme_minimal()

ggsave(file.path(TB_OUT_FIGURES, "D_conflict_MC_null_comparison.png"), p_density, width = 10, height = 6)

# P-value barplot
p_pvalues <- comparison_df %>%
  ggplot(aes(x = Null_Model, y = p_value, fill = ifelse(p_value < 0.05, "Significant", "Not Sig."))) +
  geom_col(alpha = 0.8) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 1) +
  scale_fill_manual(values = c("Significant" = "darkred", "Not Sig." = "gray")) +
  labs(title = "Conflict-Occurrence Clustering: p-values by Null Model",
       x = "Null Model", y = "p-value",
       fill = "Result") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(TB_OUT_FIGURES, "D_conflict_MC_pvalues.png"), p_pvalues, width = 8, height = 6)

# ---- 6) SAVE RESULTS ------
write_csv(comparison_df, file.path(TB_OUT_TABLES, "D_conflict_MC_pvalues.csv"))

results_detail <- tibble(
  null_model = c("Random", "Torus", "Rotation"),
  p_value = c(p_random, p_torus, p_rotation),
  observed_mean_distance_km = obs_mean_dist,
  n_conflicts = n_conflict,
  n_occurrences = n_occ
)

write_csv(results_detail, file.path(TB_OUT_TABLES, "D_conflict_MC_detailed.csv"))

cat("\n✓ Saved: D_conflict_MC_pvalues.csv\n")
cat("✓ Figures: D_conflict_MC_*.png\n")

# ---- 7) RECOMMENDATION ------
cat("\n=== RECOMMENDATION FOR MANUSCRIPT ===\n")

if (p_random < 0.05 & (p_torus > 0.05 | p_rotation > 0.05)) {
  cat("Report corrected p-value using spatial-robust null.\n")
  cat(sprintf("Current text may say: 'p < 0.001' or 'p = 1.000' (from random null)\n"))
  cat(sprintf("Corrected text should say: 'p = %.4f' or 'p > 0.05' (spatial-robust null)\n", max(p_torus, p_rotation)))
} else {
  cat("Random null is valid. Current reporting is appropriate.\n")
}

cat("\n=== STRAND D COMPLETE ===\n")
