#!/usr/bin/python
'''
Find top FHS MultiXcan hits in DGN
'''

import gzip
import sys
import os

mydir = "/gpfs/data/im-lab/nas40t2/hwheeler/trans-px/multi-transpx_results/results_allgenes/results_PCs/"
rmdir = "/gpfs/data/im-lab/nas40t2/hwheeler/trans-px/NCBI_Gene_Summaries/"
annotfile = "/gpfs/data/im-lab/nas40t2/hwheeler/trans-px/meqtl_format/gencode.v18.genes.patched_contigs.summary.protein"
dgnfile = mydir + "multi-trans-px_DGN_diff_chrs_overall_results_2017-12-11.txt.gz"
#dgnfile = rmdir + "testdgn.gz"
fhsfile = rmdir + "rm_suspect_NCBIgenes_Results_rm_bad_map_multi-trans-px_FHS_overall_results_FDR0.05_2017-12-11.txt"

genepos = dict()
for line in open(annotfile):
  (chr,strand,s1,s2,ensembl_id,gene_id,pc,k) = line.strip().split()
  genepos[ensembl_id[:15]] = gene_id + '\t' + chr[3:] + '\t' + s1 + '\t' + s2

fhsset = set()
fhsinfile = open(fhsfile)
for line in fhsinfile:
  if(line.startswith('obsgene')):
    continue
  arr = line.strip().split()
  predobs = arr[1][:15] + ":" + arr[0][:15]
  fhsset.add(predobs)
fhsinfile.close()

dgnres = dict()
for line in gzip.open(dgnfile):
  if(line.startswith('obsgene')):
    continue
  arr = line.strip().split()
  predobs = arr[1][:15] + ":" + arr[0][:15]
  if predobs in fhsset:
    dgnres[predobs] = "\t".join(arr[2:5]) #keep Rsq Fstat Pval

outfile = open(rmdir + "paper_figures_2018-10-26/TableS2_MultiXcan_results_2018-10-29.txt","w")
outfile.write("predgene\tpredname\tpredChr\tpredS1\tpredS2\tobsgene\tobsname\tobsChr\tobsS1\tobsS2\tFHS_Rsq\tFHS_Fstat\tFHS_Pval\tFHS_FDR\tFHS_BONF\tDGN_Rsq\tDGN_Fstat\tDGN_Pval\n")

for line in open(fhsfile):
  if(line.startswith('obsgene')):
    continue
  arr = line.strip().split()
  (obsgene, predgene, Rsq, Fstat, Pval, ntis, ntislm, predChr, predS1, predS2, predname, obsChr, obsS1, obsS2, obsname, FDR, BONF) = arr
  obsgene = obsgene[:15]
  predgene = predgene[:15]
  if predgene in genepos:
    predinfo = genepos[predgene]
  else: 
    predinfo = "NA\tNA\tNA\tNA"
  if obsgene in genepos:
    obsinfo = genepos[obsgene]
  else:
    obsinfo = "NA\tNA\tNA\tNA"
  po = predgene[:15] + ":" + obsgene[:15]
  fhsline = predgene[:15] + "\t" + predinfo + "\t" + obsgene + "\t" + obsinfo + "\t" + Rsq + "\t" + Fstat + "\t" + Pval + "\t" + FDR + "\t" + BONF
  if po in dgnres:
    dgnline = dgnres[po]
  else:
    dgnline = "NA\tNA\tNA"
  outfile.write(fhsline + "\t" + dgnline + "\n")
    
outfile.close()
