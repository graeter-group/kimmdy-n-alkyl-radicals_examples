#!/bin/bash

#SBATCH --job-name="BASENAME"
#SBATCH -t 02:00:00
#SBATCH --nodes=1
#SBATCH -p cascade.p
#SBATCH --ntasks=20
#SBATCH --mem=45G
##SBATCH --gres=gpu:1
#SBATCH --output=slurm/run_%a.out
#SBATCH --error=slurm/run_%a.error
#SBATCH --array=X-Y


echo "Start at: $(date)"

echo "Load modules"
#activate gromacs
ml use /hits/fast/mbm/hartmaec/sw/easybuild/modules/all/
ml load GROMACS/2022.5-plumed2.9_runtime--cuda-11.5	#actually this is gmx 2023.2 with plumed

echo "Activate python"
#activate python
ml load Python/3.10.4-GCCcore-11.3.0

echo "Activate env"
source /hits/fast/mbm/ulanovei/projects/phd/kimmdy_hexyl/.venv_kimmdy_full/bin/activate || exit


echo "run kimmdy"

pwd
kimmdy --input "kimmdy__s${SLURM_ARRAY_TASK_ID}.yml"

echo "End at: $(date)"
