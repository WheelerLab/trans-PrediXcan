#!/usr/bin/python
'''
filter MulTiXcan and PrediXcan WB rm bad map results Mappability>0.8 and no cross-mappability between pairs
to also rm likely false positives based on NCBI gene summaries
grep retro, pseudogene, or paralog
needed for Table 1 tested count and Fig 2 QQ plot
'''

import gzip
import sys
import argparse
import os

def check_arg(args=None):
    parser = argparse.ArgumentParser(description='Script to filter poor mappability genes')
    parser.add_argument('-f', '--infile',
                        help='input file',
                        required='True'
                        )
    parser.add_argument('-p', '--infilepath',
                        help='path to input file',
                        required='True'
                        )
    return parser.parse_args(args)

#retrieve command line arguments

args = check_arg(sys.argv[1:])
resfilestring = args.infile
resfilepath = args.infilepath

mydir = "/gpfs/data/im-lab/nas40t2/hwheeler/trans-px/"
rmdir = mydir + "rm_bad_mapping_genes/"
sumdir = mydir + "NCBI_Gene_Summaries/"

suspectfile = sumdir + "hgFixed.refSeqSummary_gencode.v19.annotations_suspect_list.txt.gz"

resfile = resfilepath + resfilestring


suspect = dict()
for line in gzip.open(suspectfile):
    arr = line.strip().split('\t')
    ensgene = arr[2][:15]
    suspect[ensgene] = line

outres = gzip.open(sumdir + "rm_suspect_NCBIgenes_" + resfilestring,"wb")

for line in gzip.open(resfile):
    arr = line.strip().split()
    gene1 = arr[0][:15]
    gene2 = arr[1][:15]
    if gene1 == 'snps':
        outres.write('\t'.join(arr) + '\n')
    elif gene1 not in suspect and gene2 not in suspect:
        outres.write('\t'.join(arr) + '\n')

outres.close()
