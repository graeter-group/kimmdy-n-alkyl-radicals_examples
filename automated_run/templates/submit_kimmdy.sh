#!/bin/bash

#SBATCH --job-name="alkyl_sim"
#SBATCH -t 24:00:00
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
source "$HOME/software/KIMMDY/kimmdy/.venv/bin/activate"

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

echo "run kimmdy"

pwd
kimmdy --input "kimmdy.yml"

echo "End at: $(date)"
