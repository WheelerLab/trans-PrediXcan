#!/bin/bash
#PBS -N filt.mappa
#PBS -S /bin/bash
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=16gb
#PBS -o logs/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e logs/${PBS_JOBNAME}.e${PBS_JOBID}.err
cd $PBS_O_WORKDIR

module load gcc/5.4.0
module load python/2.7.13

time python 01_filter_full_results_rm_bad_map.py -f FHSobs_v_GTExWBpred.meqtl.trans.diffchr.allres.txt.gz -p /group/im-lab/nas40t2/hwheeler/trans-px/
