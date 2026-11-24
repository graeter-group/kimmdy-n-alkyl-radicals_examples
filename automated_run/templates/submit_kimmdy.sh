#!/bin/bash

#SBATCH --job-name="alkyl_sim"
#SBATCH -t 16:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=160G
#SBATCH --gres=gpu:1
#SBATCH --partition=s.otter
#SBATCH --output=slurm/run_%j.out
#SBATCH --error=slurm/run_%j.error

echo "Start at: $(date)"

echo "Load modules"
module purge
#module load anaconda/3/2023.03
module load gromacs-plumed/2024.6
#eval "$(conda shell.bash hook)"
#conda activate "$HOME/conda-envs/uv"
KIMMDY_REPO="$HOME/software/KIMMDY/kimmdy_test"
source "$KIMMDY_REPO/.venv/bin/activate"

#export OMP_NUM_THREADS="$SLURM_CPUS_PER_TASK"

# Format SLURM_ARRAY_TASK_ID as a 2-digit zero-padded number (e.g., 00, 01, ..., 09)
TASK_DIR=$(printf "sim_%02d" "$SLURM_ARRAY_TASK_ID")

# Store parent directory path
PARENT_DIR=$(pwd)

# Create task subdirectory if it doesn't exist and enter it
mkdir -p "$TASK_DIR"
cd "$TASK_DIR" || exit 1

echo "Current directory: $(pwd)"

# Copy required files (.mdp, .gro, .top, .ndx, and kimmdy.yml) from parent directory
cp "$PARENT_DIR"/*.mdp "$PARENT_DIR"/*.gro "$PARENT_DIR"/*.top "$PARENT_DIR"/*.ndx .
cp "$PARENT_DIR"/kimmdy.yml . 2>/dev/null || true

# Record the KIMMDY code state this task actually imports (editable install).
# Compare diff_sha256 with ../../kimmdy_version.txt to detect changes since submission.
{
    echo "recorded_at: $(date -Iseconds)"
    echo "repo: $KIMMDY_REPO"
    echo "branch: $(git -C "$KIMMDY_REPO" branch --show-current)"
    echo "commit: $(git -C "$KIMMDY_REPO" rev-parse HEAD)"
    echo "diff_sha256: $(git -C "$KIMMDY_REPO" diff HEAD --submodule=diff | sha256sum | cut -d' ' -f1)"
    echo "--- git submodule status"
    git -C "$KIMMDY_REPO" submodule status
    echo "--- git status --short"
    git -C "$KIMMDY_REPO" status --short
} > kimmdy_version.txt 2>&1 || echo "WARNING: could not record KIMMDY version"
cat kimmdy_version.txt

echo "run kimmdy"

pwd
kimmdy --input "kimmdy.yml"

echo "End at: $(date)"
