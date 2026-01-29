#!/bin/bash -l

#SBATCH --qos=regular
#SBATCH --job-name=diag_feedback_mmf_cpl
#SBATCH --nodes=1
#SBATCH --output=diag_feedback_mmf_cpl.o%j
#SBATCH --time=24:00:00
#SBATCH --account=m3312
#SBATCH --constraint=cpu
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=hannah6@llnl.gov

# conda create -n diag_feedback -c conda-forge -c cdat/label/v8.2.1 python=3.7.12 cdtime cdms2 genutil cdutil psutil cartopy nco scikit-learn statsmodels

module load python/3.9-anaconda-2021.11
conda activate diag_feedback

cd /global/homes/w/whannah/Research/E3SM/pub_figs/cpl/diag_feedback_E3SM

python -u run_diag_feedback_mmf_cpl_main.py > diag_feedback_mmf_cpl.log 2>&1
