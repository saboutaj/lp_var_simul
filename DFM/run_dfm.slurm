#!/bin/bash
#SBATCH --job-name=lp_var_simul
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G
#SBATCH --time=00:30:00
#SBATCH --array=1-2
#SBATCH --mail-type=end
#SBATCH --mail-user=dakel@princeton.edu

export DGP_TYPE=G
export ESTIMAND_TYPE=ObsShock
export LAG_TYPE=4
export SAVE_PRE=/tigress/dakel/Results

echo "My SLURM_ARRAY_JOB_ID is $SLURM_ARRAY_JOB_ID."
echo "My SLURM_ARRAY_TASK_ID is $SLURM_ARRAY_TASK_ID"
echo "Executing on the machine:" $(hostname)

cd /home/dakel/lp_var_simul/DFM
/usr/licensed/bin/matlab -singleCompThread -nodisplay -nosplash -r "num_workers=$SLURM_CPUS_PER_TASK; spec_id=$SLURM_ARRAY_TASK_ID; dgp_type=$DGP_TYPE; estimand_type=$ESTIMAND_TYPE; lag_type=$LAG_TYPE; save_pre='$SAVE_PRE'; run_dfm;"