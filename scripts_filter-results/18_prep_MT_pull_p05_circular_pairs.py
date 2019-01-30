#!/usr/bin/python
'''
Retrive circular pairs and p-values for trans-px paper revision
'''
import gzip
mydir =  "/gpfs/data/im-lab/nas40t2/hwheeler/trans-px/"

#mtfile = mydir + "multi-transpx_results/results_allgenes/results_PCs/multi-trans-px_FHS_diff_chrs_overall_results_2017-12-11.txt.gz"
mtfile = mydir + "NCBI_Gene_Summaries/circular_pairs/circular_pairs_multi-trans-px_FHS_diff_chrs_overall_results_2017-12-11.txt.gz"

outfile = gzip.open(mydir + "NCBI_Gene_Summaries/circular_pairs/circular_pairs_p05_multi-trans-px_FHS_diff_chrs_overall_results_2017-12-11.txt.gz","wb")
outfile.write("pair1\tpval1\tpair2\tpval2\n")

for line in gzip.open(mtfile):
  (pair1, pval1, pair2, pval2) = line.strip().split()
  if pval1 == "pval1":
    next
  elif float(pval1) < 0.05 or float(pval2) < 0.05:
    outfile.write(pair1 + '\t' + pval1 + '\t' + pair2 + '\t' + pval2 + '\n')
    
outfile.close()
