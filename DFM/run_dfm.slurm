#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=00:50:00
 
cd /home/dakel/4/lp_var_simul/DFM
/usr/licensed/bin/matlab -singleCompThread -nodisplay -nosplash -r run_dfm