#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH --time=71:59:00
 
cd /home/dakel/4/lp_var_simul/DFM
/usr/licensed/bin/matlab -singleCompThread -nodisplay -nosplash -r run_dfm_G_Recursive