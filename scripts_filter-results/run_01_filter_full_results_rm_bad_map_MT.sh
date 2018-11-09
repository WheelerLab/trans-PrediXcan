#!/bin/bash
#PBS -N filt.mappa
#PBS -S /bin/bash
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=32gb
#PBS -o logs/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e logs/${PBS_JOBNAME}.e${PBS_JOBID}.err
cd $PBS_O_WORKDIR

module load gcc/5.4.0
module load python/2.7.13

time python 01_filter_full_results_rm_bad_map.py -f multi-trans-px_FHS_diff_chrs_overall_results_2017-12-11.txt.gz -p /group/im-lab/nas40t2/hwheeler/trans-px/multi-transpx_results/results_allgenes/results_PCs/
