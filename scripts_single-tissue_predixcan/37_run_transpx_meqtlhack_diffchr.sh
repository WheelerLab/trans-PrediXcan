#!/bin/bash
#PBS -N transpx
#PBS -S /bin/bash
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=12gb
#PBS -o logs/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e logs/${PBS_JOBNAME}.e${PBS_JOBID}.err
cd $PBS_O_WORKDIR

module load gcc/6.2.0
module load R/3.3.2

#time R --no-save < transpx_meqtlhack_diffchr.R --args /group/im-lab/nas40t2/hwheeler/trans-px/meqtl_format/ predEXP_FHS_from_GTExWBdb.txt predLOC_gtex.txt obsEXP_FHS_gtex.txt obsLOC_gtex.txt FHSobs_v_GTExWBpred
#time R --no-save < transpx_meqtlhack_diffchr.R --args /group/im-lab/nas40t2/hwheeler/trans-px/meqtl_format/ predEXP_DGN_from_GTExWBdb.txt predLOC_DGN_from_GTExWBdb.txt obsEXP_DGN_from_GTExWBdb.txt obsLOC_DGN_from_GTExWBdb.txt DGNobs_v_GTExWBpred


#time R --no-save < transpx_meqtlhack_diffchr.R --args /group/im-lab/nas40t2/hwheeler/trans-px/meqtl_format/ predEXP_FHS_GENIII_from_DGNdb.txt predLOC.txt obsEXP_FHS_GENIII.txt obsLOC.txt FHS-GENIIIobs_v_DGNpred
#time R --no-save < transpx_meqtlhack_diffchr.R --args /group/im-lab/nas40t2/hwheeler/trans-px/meqtl_format/ predEXP_FHS_OFF_from_DGNdb.txt predLOC.txt obsEXP_FHS_OFF.txt obsLOC.txt FHS-OFFobs_v_DGNpred
#time R --no-save < transpx_meqtlhack_diffchr.R --args /group/im-lab/nas40t2/hwheeler/trans-px/meqtl_format/ predEXP_FHS_GENIII_from_GTExWB.txt predLOC_gtex.txt obsEXP_FHS_GENIII_gtex.txt obsLOC_gtex.txt FHS-GENIIIobs_v_GTExWBpred
#time R --no-save < transpx_meqtlhack_diffchr.R --args /group/im-lab/nas40t2/hwheeler/trans-px/meqtl_format/ predEXP_FHS_OFF_from_GTExWB.txt predLOC_gtex.txt obsEXP_FHS_OFF_gtex.txt obsLOC_gtex.txt FHS-OFFobs_v_GTExWBpred

time R --no-save < transpx_meqtlhack_diffchr_covariates.R --args /group/im-lab/nas40t2/hwheeler/trans-px/meqtl_format/ predEXP_GTExWB_from_DGNdb.txt predLOC_GTExWB_from_DGNdb.txt obsEXP_GTExWB_from_DGNdb.txt obsLOC_GTExWB_from_DGNdb.txt GTExWBobs20PCs_v_DGNpred covfile_GTExWB_20PCs.txt
