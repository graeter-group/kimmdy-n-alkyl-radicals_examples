#!/usr/bin/env bash
set -euo pipefail

# ==========================================
# Configuration & Header Variables
# ==========================================
RUN_NAME="HAT-MACE_datasetv2_capdist8"
DRY_RUN=false # Configured for active Slurm job submission

# Internal Relative Directory Paths
TEMPLATES_DIR="templates"
SYSTEMS_DIR="systems"
RUNS_BASE_DIR="runs"

# KIMMDY checkout used by submit_kimmdy.sh (editable install, keep in sync with that script)
KIMMDY_REPO="${HOME}/software/KIMMDY/kimmdy_test"

#SYSTEMS=("propyl" "butyl" "pentyl" "hexyl" "heptyl" "octyl")
SYSTEMS=("heptyl" "octyl")

declare -A SYSTEM_RADICALS=(
    ["propyl"]="3"
    ["butyl"]="4"
    ["pentyl"]="5"
    ["hexyl"]="6"
    ["heptyl"]="7"
    ["octyl"]="8"
)

# Generate Timestamped Run Directory
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
CURRENT_RUN_DIR="${RUNS_BASE_DIR}/${RUN_NAME}_${TIMESTAMP}"

echo "=================================================="
echo " Starting Phase 1 Pipeline: ${TIMESTAMP}_${RUN_NAME}"
echo " DRY_RUN Mode: ${DRY_RUN}"
echo "=================================================="

mkdir -p "${CURRENT_RUN_DIR}"

# Copy global templates directory to current run directory
echo "Copying global templates directory to ${CURRENT_RUN_DIR}..."
cp -r "${TEMPLATES_DIR}" "${CURRENT_RUN_DIR}/"

# Record the KIMMDY code state at submission time (editable install, so the code
# is only fixed when each job starts; submit_kimmdy.sh records it again per task)
echo "Recording KIMMDY version from ${KIMMDY_REPO}..."
if git -C "${KIMMDY_REPO}" rev-parse --git-dir > /dev/null 2>&1; then
    # --submodule=diff: plugins (kimmdy-hat, ...) are submodules; include their code changes
    git -C "${KIMMDY_REPO}" diff HEAD --submodule=diff > "${CURRENT_RUN_DIR}/kimmdy_diff.patch"
    {
        echo "recorded_at: $(date -Iseconds)"
        echo "repo: ${KIMMDY_REPO}"
        echo "branch: $(git -C "${KIMMDY_REPO}" branch --show-current)"
        echo "commit: $(git -C "${KIMMDY_REPO}" rev-parse HEAD)"
        echo "commit_subject: $(git -C "${KIMMDY_REPO}" log -1 --format=%s)"
        echo "diff_sha256: $(sha256sum < "${CURRENT_RUN_DIR}/kimmdy_diff.patch" | cut -d' ' -f1)"
        echo "--- git submodule status (commit per plugin, '+' = differs from pinned)"
        git -C "${KIMMDY_REPO}" submodule status
        echo "--- git status --short"
        git -C "${KIMMDY_REPO}" status --short
    } > "${CURRENT_RUN_DIR}/kimmdy_version.txt"
else
    echo "WARNING: ${KIMMDY_REPO} is not a git repository, KIMMDY version not recorded"
fi

# --------------------------------------------------
# System Setup & Slurm Job Submission
# --------------------------------------------------
echo ""
echo "--- Setting up Systems and Submitting Slurm Jobs ---"

for sys in "${SYSTEMS[@]}"; do
    SYS_RUN_DIR="${CURRENT_RUN_DIR}/${sys}"
    mkdir -p "${SYS_RUN_DIR}"
    # Slurm does not create the --output/--error directory itself
    mkdir -p "${SYS_RUN_DIR}/slurm"

    echo "-> Populating run directory for: ${sys}"
    
    # 1. Copy system topology/assets from systems dir into system run dir
    cp -r "${SYSTEMS_DIR}/${sys}/"* "${SYS_RUN_DIR}/"
    
    # 2. Copy templates directory into newly created system dir
    cp -r "${TEMPLATES_DIR}/"* "${SYS_RUN_DIR}/"

    # 3. Generate kimmdy configuration files
    RAD_VAL="${SYSTEM_RADICALS[$sys]}"
    
    # Replace radicals: 'Y' with radicals: '<RAD_VAL>' to create base kimmdy.yml
    sed "s/radicals: 'Y'/radicals: '${RAD_VAL}'/" "${SYS_RUN_DIR}/template_kimmdy.yml" > "${SYS_RUN_DIR}/kimmdy.yml"
     
    # 4. Submit Slurm jobs with dependency control
    (
        cd "${SYS_RUN_DIR}"
        
        if [ "${DRY_RUN}" = true ]; then
            echo "   [DRY_RUN] Would submit submit_gmx.sh in ${SYS_RUN_DIR}"
            echo "   [DRY_RUN] Would submit array 1-16 of submit_kimmdy.sh with dependency 'afterok' on submit_gmx.sh job ID"
        else
            echo "   Submitting submit_gmx.sh..."
            GMX_JOB_ID=$(sbatch submit_gmx.sh | awk '{print $NF}')
            echo "   -> Submitted submit_gmx.sh with Job ID: ${GMX_JOB_ID}"
            
            echo "   Submitting run_kimmdy.sh array (16 tasks) dependent on Job ID ${GMX_JOB_ID}..."
            sbatch --array=0-15 --dependency=afterok:"${GMX_JOB_ID}" submit_kimmdy.sh
        fi
    )
done

echo ""
echo "=================================================="
echo " Submissions complete."
echo " Run Directory: ${CURRENT_RUN_DIR}"
echo "=================================================="
