#!/bin/bash
# ============================================================================
# check_outputs.sh
# Verify that all critical-analysis strands produced expected output files
# ============================================================================

cd "$(dirname "$0")/.."

echo "================================================================================"
echo "Output Verification Check"
echo "================================================================================"
echo ""

OUT_DIR="outputs"
FIGURES="$OUT_DIR/figures"
TABLES="$OUT_DIR/tables"
RDS="$OUT_DIR/rds_files"
RASTERS="$OUT_DIR/rasters"
CONEFOR="$OUT_DIR/derived/conefor"

# ---- STRAND C ----
echo "STRAND C: Sensitivity Analyses"
echo "  C1 (MAX_TSS):"
[ -f "$TABLES/C1_maxTSS_sensitivity.csv" ] && echo "    OK C1_maxTSS_sensitivity.csv" || echo "    MISSING C1_maxTSS_sensitivity.csv"
[ -f "$FIGURES/C1_maxTSS_n_patches.png" ] && echo "    OK C1_maxTSS_n_patches.png" || echo "    MISSING C1_maxTSS_n_patches.png"

echo "  C2 (PA ratio):"
[ -f "$TABLES/C2_PA_ratio_sensitivity.csv" ] && echo "    OK C2_PA_ratio_sensitivity.csv" || echo "    MISSING C2_PA_ratio_sensitivity.csv"
[ -f "$FIGURES/C2_PA_ratio_TSS.png" ] && echo "    OK C2_PA_ratio_TSS.png" || echo "    MISSING C2_PA_ratio_TSS.png"

echo "  C3 (c-sensitivity):"
[ -f "$TABLES/C3_c_sensitivity_full.csv" ] && echo "    OK C3_c_sensitivity_full.csv" || echo "    MISSING C3_c_sensitivity_full.csv"
[ -f "$TABLES/C3_manuscript_summary.csv" ] && echo "    OK C3_manuscript_summary.csv" || echo "    MISSING C3_manuscript_summary.csv"
[ -f "$FIGURES/C3_c_sensitivity/fig_C3_PC_loss.png" ] && echo "    OK fig_C3_PC_loss.png" || echo "    MISSING fig_C3_PC_loss.png"
[ -f "$FIGURES/C3_c_sensitivity/fig_C3_PC_vs_ECA.png" ] && echo "    OK fig_C3_PC_vs_ECA.png" || echo "    MISSING fig_C3_PC_vs_ECA.png"
[ -f "$FIGURES/C3_c_sensitivity/fig_C3_heatmap.png" ] && echo "    OK fig_C3_heatmap.png" || echo "    MISSING fig_C3_heatmap.png"
echo ""

# ---- STRAND D ----
echo "STRAND D: Conflict MC Spatial Check"
echo "  Expected files:"
[ -f "$TABLES/D_conflict_MC_pvalues.csv" ] && echo "    OK D_conflict_MC_pvalues.csv" || echo "    MISSING D_conflict_MC_pvalues.csv"
[ -f "$TABLES/D_conflict_MC_detailed.csv" ] && echo "    OK D_conflict_MC_detailed.csv" || echo "    MISSING D_conflict_MC_detailed.csv"
[ -f "$FIGURES/D_conflict_MC_null_comparison.png" ] && echo "    OK D_conflict_MC_null_comparison.png" || echo "    MISSING D_conflict_MC_null_comparison.png"
[ -f "$FIGURES/D_conflict_MC_pvalues.png" ] && echo "    OK D_conflict_MC_pvalues.png" || echo "    MISSING D_conflict_MC_pvalues.png"
echo ""

# ---- STRAND E ----
echo "STRAND E: Source Independence"
echo "  Expected files:"
[ -f "$TABLES/E_source_independence_summary.csv" ] && echo "    OK E_source_independence_summary.csv" || echo "    MISSING E_source_independence_summary.csv"
[ -f "$TABLES/E_source_independence_pairs.csv" ] && echo "    OK E_source_independence_pairs.csv" || echo "    MISSING E_source_independence_pairs.csv"
[ -f "$FIGURES/E_source_independence/fig_E_coincidence_map.png" ] && echo "    OK fig_E_coincidence_map.png" || echo "    MISSING fig_E_coincidence_map.png"
echo ""

# ---- STRAND F ----
echo "STRAND F: MESS Extrapolation"
echo "  Expected files:"
[ -f "$TABLES/F_mess_extrapolation.csv" ] && echo "    OK F_mess_extrapolation.csv" || echo "    MISSING F_mess_extrapolation.csv"
[ -f "$TABLES/F_mess_extrapolation_summary.csv" ] && echo "    OK F_mess_extrapolation_summary.csv" || echo "    MISSING F_mess_extrapolation_summary.csv"
[ -f "$FIGURES/F_mess_extrapolation/fig_F_mess_maps.png" ] && echo "    OK fig_F_mess_maps.png" || echo "    MISSING fig_F_mess_maps.png"
[ -f "$FIGURES/F_mess_extrapolation/fig_F_nonanalog_bar.png" ] && echo "    OK fig_F_nonanalog_bar.png" || echo "    MISSING fig_F_nonanalog_bar.png"
echo ""

echo "================================================================================"
echo "Check Job Logs"
echo "================================================================================"
echo ""
echo "Recent log files:"
ls -lht scripts/logs/ | head -10
echo ""

echo "================================================================================"
echo "Next Steps"
echo "================================================================================"
echo ""
echo "When ALL outputs present:"
echo "  1. Download outputs to local machine"
echo "  2. Review key figures and tables"
echo "  3. Prepare manuscript rewrite"
echo ""
echo "Download command:"
echo "  scp -r <user>@<cluster-host>:/path/to/tr_bear_connect/outputs ."
echo ""
echo "================================================================================"
