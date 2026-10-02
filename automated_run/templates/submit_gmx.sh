#!/bin/bash
#SBATCH --job-name="gmx_sd_sims"
#SBATCH -t 2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=16            # Match this to NUM_SIMS below
#SBATCH --cpus-per-task=1
#SBATCH --mem=160G
#SBATCH --gres=gpu:1
#SBATCH --partition=s.otter
#SBATCH --output=slurm/run_%j.out
#SBATCH --error=slurm/run_%j.error

# Set the total number of parallel simulations
NUM_SIMS=16

echo "Start at: $(date)"
echo "Load modules"

module purge
module load gromacs-plumed/2024.6

# Unset OMP_NUM_THREADS to rely strictly on -ntomp flag
unset OMP_NUM_THREADS
export CUDA_VISIBLE_DEVICES=""

# Ensure slurm log directory exists
mkdir -p slurm

# Define base input files
MDP="sd.mdp"
GRO="mol.gro"
TOP="topol.top"
NDX="index.ndx"

# Array to collect directory paths for -multidir
DIRS=()

echo "=== Preparing ${NUM_SIMS} directories and compiling TPR files ==="

for (( i=0; i<NUM_SIMS; i++ )); do
  # Format directory name with zero-padding (e.g., sim_01, sim_02 ... sim_16)
  DIR=$(printf "sim_%02d" "${i}")
  mkdir -p "${DIR}"
  DIRS+=("${DIR}")

  # Compile equilibrium.tpr directly inside each simulation subdirectory
  gmx grompp \
    -f "${MDP}" \
    -c "${GRO}" \
    -p "${TOP}" \
    -n "${NDX}" \
    -o "${DIR}/equilibrium.tpr" \
    -po "${DIR}/equilibrium.mdp" \
    -maxwarn 1
done

echo "=== Starting ${NUM_SIMS} parallel single-core simulations with mpirun ==="

# Join directory names into a space-separated string
MULTIDIR_ARGS="${DIRS[*]}"

# Run all simulations in parallel across NUM_SIMS ranks
# Note: Pointing -s and -deffnm to 'equilibrium' to match grompp output name
mpirun -np "${NUM_SIMS}" gmx_mpi mdrun \
  -multidir ${MULTIDIR_ARGS} \
  -s equilibrium.tpr \
  -deffnm equilibrium \
  -ntomp 1

echo "=== All ${NUM_SIMS} simulations completed successfully ==="
