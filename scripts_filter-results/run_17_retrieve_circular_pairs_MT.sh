#!/bin/bash
#PBS -N circleMT
#PBS -S /bin/bash
#PBS -l walltime=240:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=64gb
#PBS -o logs/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e logs/${PBS_JOBNAME}.e${PBS_JOBID}.err
cd $PBS_O_WORKDIR

module load gcc/5.4.0
module load python/2.7.13

time python 17_retrieve_circular_pairs_MT.py
