#!/usr/bin/python
'''
Retrive all data in Table S3 format for Zenodo in trans-px paper revision
'''
import gzip
mydir =  "/gpfs/data/im-lab/nas40t2/hwheeler/trans-px/"

#mtfile = mydir + "multi-transpx_results/results_allgenes/results_PCs/multi-trans-px_FHS_diff_chrs_overall_results_2017-12-11.txt.gz"
mtfile = mydir + "NCBI_Gene_Summaries/rm_suspect_NCBIgenes_Results_rm_bad_map_multi-trans-px_FHS_diff_chrs_overall_results_2017-12-11.txt.gz"

outfile = gzip.open(mydir + "NCBI_Gene_Summaries/DatasetS2_trans-MultiXcan_multi-tissue_FHS_results.txt.gz","wb")
outfile.write("predgene\tpredname\tpredChr\tpredS1\tpredS2\tobsgene\tobsname\tobsChr\tobsS1\tobsS2\tFHS_Rsq\tFHS_Fstat\tFHS_pval\n")

for line in gzip.open(mtfile):
  (obs, pred, Rsq, Fstat, pval, ntis, ntislm, predChr, predS1, predS2, predname, obsChr, obsS1, obsS2, obsname) = line.strip().split()
  pred = pred[:15] #rm .version number
  obs = obs[:15] #rm .version number
  outfile.write(pred + '\t' + predname + '\t' + predChr + '\t' + predS1 + '\t' + predS2 + '\t' + obs + '\t' + obsname + '\t' + obsChr + '\t' + obsS1 + '\t' + obsS2 + '\t' + Rsq + '\t' + Fstat + '\t' + pval + '\n')

outfile.close()
