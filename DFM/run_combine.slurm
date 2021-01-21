#!/bin/bash
#SBATCH --job-name=combine_spec
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1:00:00
#SBATCH --array=1
#SBATCH --mail-type=end
#SBATCH --mail-user=dakel@princeton.edu

export SPEC_ID_ARRAY=1:30
export DGP_TYPE_LIST=MP
export ESTIMAND_TYPE_LIST=IV
export LAG_TYPE_LIST=8
export SAVE_PRE=/tigress/dakel/Results

echo "My SLURM_ARRAY_JOB_ID is $SLURM_ARRAY_JOB_ID."
echo "My SLURM_ARRAY_TASK_ID is $SLURM_ARRAY_TASK_ID"
echo "Executing on the machine:" $(hostname)

cd /home/dakel/lp_var_simul/DFM/Subroutines
/usr/licensed/bin/matlab -singleCompThread -nodisplay -nosplash -r "clear all; spec_id_array=$SPEC_ID_ARRAY; dgp_type_list={'$DGP_TYPE_LIST'}; estimand_type_list={'$ESTIMAND_TYPE_LIST'}; lag_type_list={$LAG_TYPE_LIST}; save_pre='$SAVE_PRE'; run_combine;"