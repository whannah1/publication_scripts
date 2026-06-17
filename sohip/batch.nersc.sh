#!/bin/bash
#SBATCH --account=m4842
#SBATCH --constraint=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --qos=regular
#SBATCH --time=02:00:00
###SBATCH --job-name=sohip_batch_curtain_wgt_method_inv_dist
#SBATCH --job-name=sohip_batch_curtain_wgt_method_area
#SBATCH --output=/global/homes/w/whannah/publication_scripts/sohip/slurm_logs/batch.slurm-%A.out
#SBATCH --mail-user=hannah6@llnl.gov
#SBATCH --mail-type=END,FAIL

# command to submit:  sbatch batch.nersc.sh

source activate ux_env

date

SRC_FILE=/global/homes/w/whannah/publication_scripts/sohip/code/calc.curtain.v2.py
# LOG_FILE=/global/homes/w/whannah/publication_scripts/sohip/code/calc.curtain.v2.wgt_method_inv_dist.out

# LOG_FILE=/global/homes/w/whannah/publication_scripts/sohip/code/calc.curtain.v2.wgt_method_inv_dist.out

time python -u ${SRC_FILE}

date
