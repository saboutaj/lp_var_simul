#!/bin/bash
#SBATCH --job-name=report_fig
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=128G
#SBATCH --time=24:00:00
#SBATCH --array=1

export LAGS_SELECT=1
export EXPER_SELECT=1
export ROOTFOLDER=/tigress/dakel/Results

echo "My SLURM_ARRAY_JOB_ID is $SLURM_ARRAY_JOB_ID."
echo "My SLURM_ARRAY_TASK_ID is $SLURM_ARRAY_TASK_ID"
echo "Executing on the machine:" $(hostname)

module purge
module load matlab/R2020b

cd /home/dakel/lp_var_simul/DFM/Reporting
matlab -singleCompThread -nodisplay -nosplash -r "clear all; lags_select=$LAGS_SELECT; exper_select=$EXPER_SELECT; rootfolder='$ROOTFOLDER'; run_plot_dgp;"
matlab -singleCompThread -nodisplay -nosplash -r "clear all; lags_select=$LAGS_SELECT; exper_select=$EXPER_SELECT; rootfolder='$ROOTFOLDER'; run_plot_diagn;"
matlab -singleCompThread -nodisplay -nosplash -r "clear all; lags_select=$LAGS_SELECT; exper_select=$EXPER_SELECT; rootfolder='$ROOTFOLDER'; run_plot_loss;"
matlab -singleCompThread -nodisplay -nosplash -r "clear all; lags_select=$LAGS_SELECT; exper_select=$EXPER_SELECT; rootfolder='$ROOTFOLDER'; run_plot_tradeoff;"