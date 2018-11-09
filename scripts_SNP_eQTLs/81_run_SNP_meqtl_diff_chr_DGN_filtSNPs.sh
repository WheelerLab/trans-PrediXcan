#!/bin/bash
#PBS -N filt-meqtl
#PBS -S /bin/bash
#PBS -l walltime=96:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=20gb
#PBS -o logs/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e logs/${PBS_JOBNAME}.e${PBS_JOBID}.err
cd $PBS_O_WORKDIR

module load gcc/6.2.0
module load R/3.3.2

time R --no-save < snp_meqtl_diffchr.R --args /group/im-lab/nas40t2/hwheeler/trans-px/DGN_dosages/meqtl_input/ chr1-22_pruned_FHStrans_meqtl_SNPfile.txt.gz chr1-22_pruned_FHStrans_meqtl_LOCfile.txt /group/im-lab/nas40t2/hwheeler/trans-px/meqtl_format/ obsEXP_DGN_from_GTExWBdb.txt obsLOC_DGN_from_GTExWBdb.txt SNP_DGN_prunedFHStrans
