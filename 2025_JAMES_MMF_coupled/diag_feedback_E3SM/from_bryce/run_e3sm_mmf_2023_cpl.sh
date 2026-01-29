#!/bin/bash -l

#SBATCH --qos=regular
#SBATCH --job-name=mmf_cpl2023
#SBATCH --nodes=1
#SBATCH --output=mmf_cpl2023.o%j
#SBATCH --time=15:00:00
#SBATCH --account=m3312
#SBATCH --constraint=cpu
#SBATCH --mail-type=end,fail
#SBATCH --mail-user=bryce.harrop@pnnl.gov

module load python/3.9-anaconda-2021.11
conda activate diag_feedback

cd /global/homes/b/beharrop/diag_feedback_e3sm/diag_feedback_E3SM/

python main_mmf_2023cpl.py >LOG.MMF 2>&1
