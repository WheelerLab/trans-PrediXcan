#!/usr/bin/python
'''
Find top FHS hits in DGN, add chr positions and gene names
'''

import gzip
import sys
import os

mydir = "/gpfs/data/im-lab/nas40t2/hwheeler/trans-px/"

annotfile = mydir + "meqtl_format/gencode.v18.genes.patched_contigs.summary.protein"
dgnfile = mydir + "DGNobs_v_GTExWBpred.meqtl.trans.diffchr.allres.txt.gz"
fhsfile = mydir + "NCBI_Gene_Summaries/rm_suspect_NCBIgenes_Results_rm_bad_map_FHSobs_v_GTExWBpred.meqtl.trans.diffchr.FDR0.05.txt"

genepos = dict()
for line in open(annotfile):
  (chr,strand,s1,s2,ensembl_id,gene_id,pc,k) = line.strip().split()
  genepos[ensembl_id[:15]] = gene_id + '\t' + chr[3:] + '\t' + s1 + '\t' + s2

dgnres = dict()
for line in gzip.open(dgnfile):
  (snps, gene, statistic, pvalue, FDR, beta) = line.strip().split()
  gene1 = snps[:15]
  gene2 = gene[:15]
  genes = str(gene1) + ':' + str(gene2)
  dgnres[genes] = statistic + '\t' + beta + '\t' + pvalue + '\t' + FDR
  
outfile = open(mydir + "NCBI_Gene_Summaries/paper_figures_2018-10-26/TableS1_WHLBLD_results_2018-10-29.txt","w")
outfile.write("tissue\tpredgene\tpredname\tpredChr\tpredS1\tpredS2\tobsgene\tobsname\tobsChr\tobsS1\tobsS2\tFHS_stat\tFHS_beta\tFHS_pval\tFHS_FDR\tDGN_stat\tDGN_beta\tDGN_pval\tDGN_FDR\n")

for line in open(fhsfile):
  if(line.startswith('snps')):
    continue
  (FHS_pred,FHS_obs,FHS_stat,FHS_pval,FHS_FDR,FHS_beta,newFDR,BONF) = line.strip().split()
  if FHS_pred in genepos:
    predinfo = genepos[FHS_pred]
  else: 
    predinfo = "NA\tNA\tNA\tNA"
  if FHS_obs in genepos:
    obsinfo = genepos[FHS_obs]
  else:
    obsinfo = "NA\tNA\tNA\tNA"
  genepair = FHS_pred + ':' + FHS_obs
  if genepair in dgnres:
    dgninfo = dgnres[genepair]
  else:
    dgninfo = "NA\tNA\tNA\tNA"
  outline = "WHLBLD\t" + FHS_pred + "\t" + predinfo + "\t" + FHS_obs + "\t" + obsinfo + "\t" + FHS_stat + "\t" + FHS_beta + "\t" + FHS_pval + "\t" + newFDR + "\t" + dgninfo
  outfile.write(outline + "\n")
