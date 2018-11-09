#!/bin/bash
#PBS -N transpx-cor
#PBS -S /bin/bash
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=12gb
#PBS -o logs/${PBS_JOBNAME}.o${PBS_JOBID}.log
#PBS -e logs/${PBS_JOBNAME}.e${PBS_JOBID}.err
cd $PBS_O_WORKDIR

module load gcc/6.2.0
module load R/3.3.2

#fix FDRs 2017-10-05 in trans-diffChrs_FHSobs_v_GTExWBpred_pval0.05_2017-07-03.txt.gz and trans-diffChrs_DGNobs_v_GTExWBpred_pval0.05_2017-07-03.txt.gz
#need to rerun meqtl to fix FDRs in OFF and GENIII analyses 

#time R --no-save < transpx_calc_cormat.R --args /group/im-lab/nas40t2/hwheeler/trans-px/meqtl_format/ predEXP_FHS_from_GTExWBdb.txt obsEXP_FHS_gtex.txt FHS_gene_annot_gtex.txt FHSobs_v_GTExWBpred.meqtl.trans.diffchr.allres.txt.gz FHSobs_v_GTExWBpred
#time R --no-save < transpx_calc_cormat.R --args /group/im-lab/nas40t2/hwheeler/trans-px/meqtl_format/ predEXP_DGN_from_GTExWBdb.txt obsEXP_DGN_from_GTExWBdb.txt gene_annot_DGN_from_GTExWBdb.txt DGNobs_v_GTExWBpred.meqtl.trans.diffchr.allres.txt.gz DGNobs_v_GTExWBpred

#time R --no-save < transpx_calc_cormat.R --args /group/im-lab/nas40t2/hwheeler/trans-px/meqtl_format/ predEXP_FHS_GENIII_from_DGNdb.txt obsEXP_FHS_GENIII.txt FHS_gene_annot.txt FHS-GENIIIobs_v_DGNpred.meqtl.trans.diffchr.allres.txt.gz FHS-GENIIIobs_v_DGNpred
#time R --no-save < transpx_calc_cormat.R --args /group/im-lab/nas40t2/hwheeler/trans-px/meqtl_format/ predEXP_FHS_OFF_from_DGNdb.txt obsEXP_FHS_OFF.txt FHS_gene_annot.txt FHS-OFFobs_v_DGNpred.meqtl.trans.diffchr.allres.txt.gz FHS-OFFobs_v_DGNpred 
time R --no-save < transpx_calc_cormat.R --args /group/im-lab/nas40t2/hwheeler/trans-px/meqtl_format/ predEXP_FHS_GENIII_from_GTExWB.txt obsEXP_FHS_GENIII_gtex.txt FHS_gene_annot_gtex.txt FHS-GENIIIobs_v_GTExWBpred.meqtl.trans.diffchr.allres.txt.gz FHS-GENIIIobs_v_GTExWBpred
time R --no-save < transpx_calc_cormat.R --args /group/im-lab/nas40t2/hwheeler/trans-px/meqtl_format/ predEXP_FHS_OFF_from_GTExWB.txt obsEXP_FHS_OFF_gtex.txt FHS_gene_annot_gtex.txt FHS-OFFobs_v_GTExWBpred.meqtl.trans.diffchr.allres.txt.gz FHS-OFFobs_v_GTExWBpred
