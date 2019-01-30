#!/usr/bin/python
'''
Retrive circular pairs and p-values for trans-px paper revision
'''
import gzip
mydir =  "/gpfs/data/im-lab/nas40t2/hwheeler/trans-px/"

wbfile = mydir + "NCBI_Gene_Summaries/rm_suspect_NCBIgenes_Results_rm_bad_map_FHSobs_v_GTExWBpred.meqtl.trans.diffchr.allres.txt.gz"
annotfile = mydir + "meqtl_format/gencode.v18.genes.patched_contigs.summary.protein"

genepos = dict()
for line in open(annotfile):
  (chr,strand,s1,s2,ensembl_id,gene_id,pc,k) = line.strip().split()
  genepos[ensembl_id[:15]] = gene_id + '\t' + chr[3:] + '\t' + s1 + '\t' + s2

outfile = gzip.open(mydir + "NCBI_Gene_Summaries/DatasetS1_trans-PrediXcan_whole_blood_FHS_results.txt.gz","wb")
outfile.write("predgene\tpredname\tpredChr\tpredS1\tpredS2\tobsgene\tobsname\tobsChr\tobsS1\tobsS2\tFHS_t-stat\tFHS_beta\tFHS_pval\n")

for line in gzip.open(wbfile):
  if(line.startswith('snps')):
    continue
  (pred, obs, tstat, pval, fdr, beta) = line.strip().split()
  pred = pred[:15] #rm .version number
  obs = obs[:15] #rm .version number
  if pred in genepos:
    predinfo = genepos[pred]
  else: 
    predinfo = "NA\tNA\tNA\tNA"
  if obs in genepos:
    obsinfo = genepos[obs]
  else:
    obsinfo = "NA\tNA\tNA\tNA"
  outfile.write(pred + '\t' + predinfo + '\t' + obs + '\t' + obsinfo + '\t' + tstat + '\t' + beta + '\t' + pval + '\n')

outfile.close()
