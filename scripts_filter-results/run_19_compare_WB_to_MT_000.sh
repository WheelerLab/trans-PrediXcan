#!/bin/bash
#PBS -N wb-mt.000
#PBS -S /bin/bash
#PBS -l walltime=240:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=16gb
#PBS -o logs/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e logs/${PBS_JOBNAME}.e${PBS_JOBID}.err
cd $PBS_O_WORKDIR

module load gcc/6.2.0
module load R/3.3.2

time R --no-save < 19_compare_WB_to_MT.R --args 000 FHS
