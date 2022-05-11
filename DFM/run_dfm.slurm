#!/bin/bash
#SBATCH --job-name=lp_var_simul
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25
#SBATCH --mem-per-cpu=4G
#SBATCH --time=24:00:00
#SBATCH --array=1

export DGP_TYPE=G
export ESTIMAND_TYPE=ObsShock
export LAG_TYPE=4
export SAVE_PRE=/tigress/dakel/Results

echo "My SLURM_ARRAY_JOB_ID is $SLURM_ARRAY_JOB_ID."
echo "My SLURM_ARRAY_TASK_ID is $SLURM_ARRAY_TASK_ID"
echo "Executing on the machine:" $(hostname)

module purge
module load matlab/R2020b

cd /home/dakel/lp_var_simul/DFM
matlab -nodisplay -nosplash -r "clear all; num_workers=$SLURM_CPUS_PER_TASK; spec_id=$SLURM_ARRAY_TASK_ID; dgp_type='$DGP_TYPE'; estimand_type='$ESTIMAND_TYPE'; lag_type=$LAG_TYPE; save_pre='$SAVE_PRE'; run_dfm;"