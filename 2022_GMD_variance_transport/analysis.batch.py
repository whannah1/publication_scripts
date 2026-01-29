#!/bin/bash
#SBATCH -A CLI115
###SBATCH -J analysis_batch
#SBATCH -N 1
#SBATCH -p gpu
#SBATCH -t 6:00:00
#SBATCH -o /ccs/home/hannah6/Research/E3SM/pub_figs/vtval/slurm_logs/analysis.batch.slurm-%A.out
#SBATCH --mail-user=hannah6@llnl.gov
#SBATCH --mail-type=END,FAIL

# commands to submit:  
# sbatch -J analysis_batch     analysis.batch.py
# sbatch -J analysis_batch_chx analysis.batch.py
# sbatch -J analysis_batch_cli analysis.batch.py
# sbatch -J analysis_batch_var analysis.batch.py
# sbatch -J analysis_batch_wfs analysis.batch.py

source activate pyn_env

TOPDIR=${HOME}/Research/E3SM/pub_figs/vtval

date

# time python -u code/calculate.pattern_frequency.v1.py
# time python -u code/calculate.pattern_frequency.v2.py
# time python -u code/calculate.pattern_frequency.v3.py

# time python -u ${TOPDIR}/code/calculate.pattern_frequency.v1.py > ${TOPDIR}/calc.chx.freq.v1.gpm.out 
# time python -u ${TOPDIR}/code/calculate.pattern_frequency.v2.py > ${TOPDIR}/calc.chx.freq.v2.gpm.out

# sbatch -J analysis_batch_cli_prc analysis.batch.py
# time python -u  ${TOPDIR}/code/F01-clim-map.py -v PRECT > ${TOPDIR}/clim_map_prc.out

# sbatch -J analysis_batch_cli_lwp analysis.batch.py
# time python -u  ${TOPDIR}/code/F01-clim-map.py -v TGCLDLWP > ${TOPDIR}/clim_map_lwp.out

# time python -u  ${TOPDIR}/code/FXX-variance-map.py    > ${TOPDIR}/variance_map.out
# time python -u  ${TOPDIR}/code/FXX-wk-wave-spectra.py > ${TOPDIR}/wk.out


# BVARIN=MMF_VT_T     ; sbatch -J analysis_batch_var1 --export=BVAR=$BVARIN analysis.batch.py
# BVARIN=OMEGA        ; sbatch -J analysis_batch_var2 --export=BVAR=$BVARIN analysis.batch.py
# BVARIN=MMF_VT_TEND_T; sbatch -J analysis_batch_var3 --export=BVAR=$BVARIN analysis.batch.py
# BVARIN=MMF_VT_TLS   ; sbatch -J analysis_batch_var4 --export=BVAR=$BVARIN analysis.batch.py

# BVARIN=T ; sbatch -J analysis_batch_var_T --export=BVAR=$BVARIN analysis.batch.py
# BVARIN=Q ; sbatch -J analysis_batch_var_Q --export=BVAR=$BVARIN analysis.batch.py
# BVARIN=U ; sbatch -J analysis_batch_var_U --export=BVAR=$BVARIN analysis.batch.py
# time python -u code/FXX-bin-VT.py -v $BVAR > ${TOPDIR}/bin_${BVAR}.out

# sbatch -J analysis_batch_bin analysis.batch.py
time python -u code/F08-bin-VT-tracer.py > ${TOPDIR}/bin_var.out

# sbatch -J analysis_batch_zmv analysis.batch.py
# time python -u code/FXX-zonal-mean-VT-tracer.py > ${TOPDIR}/zonal-mean-VT.out

date

