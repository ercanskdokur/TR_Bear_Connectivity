#!/bin/bash
# ============================================================================
# master_submit_all.sh
# TR_Bear_Connectivity: Parallel Submission of the Robustness/Sensitivity
# Analysis Strands (C, D, E, F)
# ============================================================================

set -e  # Exit on error

cd "$(dirname "$0")"

echo "================================================================================"
echo "TR_Bear_Connectivity: Sensitivity/Robustness Analysis Submission"
echo "================================================================================"
echo ""
echo "Submitting strands:"
echo "  STRAND C:  Sensitivity analyses x3 (MAX_TSS, PA ratio, c-sensitivity)"
echo "  STRAND D:  Conflict MC spatial null"
echo "  STRAND E:  Source independence"
echo "  STRAND F:  MESS extrapolation"
echo ""
echo "================================================================================"

# ---- STRAND C: Sensitivity (parallel) ----
echo ""
echo "[1/4] Submitting STRAND C (sensitivity x3 parallel)..."
echo "      Submitting C1 (MAX_TSS)..."
JOB_C1=$(sbatch C1_maxTSS_sensitivity.slurm | awk '{print $4}')
echo "      Job ID: $JOB_C1"

echo "      Submitting C2 (PA ratio)..."
JOB_C2=$(sbatch C2_pseudoabsence_sensitivity.slurm | awk '{print $4}')
echo "      Job ID: $JOB_C2"

echo "      Submitting C3 (c-sensitivity)..."
JOB_C3=$(sbatch C3_c_sensitivity_full.slurm | awk '{print $4}')
echo "      Job ID: $JOB_C3"

# ---- STRAND D: Conflict MC ----
echo ""
echo "[2/4] Submitting STRAND D (conflict MC spatial null)..."
JOB_D=$(sbatch D_conflict_mc_spatial.slurm | awk '{print $4}')
echo "      Job ID: $JOB_D"

# ---- STRAND E: Source independence ----
echo ""
echo "[3/4] Submitting STRAND E (source independence)..."
JOB_E=$(sbatch E_source_independence.slurm | awk '{print $4}')
echo "      Job ID: $JOB_E"

# ---- STRAND F: MESS extrapolation ----
echo ""
echo "[4/4] Submitting STRAND F (MESS extrapolation)..."
JOB_F=$(sbatch F_mess_extrapolation.slurm | awk '{print $4}')
echo "      Job ID: $JOB_F"

echo ""
echo "================================================================================"
echo "ALL JOBS SUBMITTED"
echo "================================================================================"
echo ""
echo "Job IDs:"
echo "  C1:  $JOB_C1  (MAX_TSS sensitivity)"
echo "  C2:  $JOB_C2  (PA ratio sensitivity)"
echo "  C3:  $JOB_C3  (c-sensitivity)"
echo "  D:   $JOB_D   (Conflict MC)"
echo "  E:   $JOB_E   (Source independence)"
echo "  F:   $JOB_F   (MESS extrapolation)"
echo ""
echo "Monitor progress:"
echo "  squeue --user=\$USER"
echo ""
echo "Check specific job:"
echo "  sstat -j [JOBID]"
echo ""
echo "View logs:"
echo "  tail -f logs/*_{JOBID}.log"
echo ""
echo "When all complete, run:"
echo "  bash check_outputs.sh"
echo ""
echo "================================================================================"
