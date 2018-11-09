#!/usr/bin/python
'''
filter Matrix eQTL SNP file using SNPlist_FHS.meqtl.trans.diffchr.FDR0.05.txt 
'''

import gzip
import sys
import argparse
import os

mydir = "/group/im-lab/nas40t2/hwheeler/trans-px/DGN_dosages/"
#prunefile ="/group/im-lab/nas40t2/hwheeler/trans-px/FHS_dosages/chr1-22.prune.50_5_0.5.in"
prunefile = "/group/im-lab/nas40t2/hwheeler/trans-px/SNPlist_FHS.meqtl.trans.diffchr.FDR0.05.txt"

locfile = mydir + "meqtl_input/chr1-22_hapmap_meqtl_LOCfile.txt"
snpfile = mydir + "meqtl_input/chr1-22_hapmap_meqtl_SNPfile.txt.gz"

prunedsnps = set()
for line in open(prunefile):
    arr = line.strip().split()
    prunedsnps.add(arr[0])

outsnp = gzip.open(mydir + "meqtl_input/chr1-22_pruned_FHStrans_meqtl_SNPfile.txt.gz","wb")
outloc = open(mydir + "meqtl_input/chr1-22_pruned_FHStrans_meqtl_LOCfile.txt","w")

for line in open(locfile):
    arr = line.strip().split()
    rs = arr[0]
    if rs == 'snp':
        outloc.write(' '.join(arr) + '\n')
    elif rs in prunedsnps:
        outloc.write(' '.join(arr) + '\n')

for line in gzip.open(snpfile):
    arr = line.strip().split()
    rs = arr[0]
    if rs == 'id':
        outsnp.write(' '.join(arr) + '\n')
    elif rs in prunedsnps:
        outsnp.write(' '.join(arr) + '\n')

            
outsnp.close()
outloc.close()
