#!/usr/bin/python
'''
filter MulTiXcan and PrediXcan WB results Mappability>0.8 and no cross-mappability between pairs
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

mydir = "/group/im-lab/nas40t2/hwheeler/trans-px/"
mapdir = mydir + "hg19_GENCODE19/"
rmdir = mydir + "rm_bad_mapping_genes/"

genemapfile = mapdir + "hg19_gene_mappability.txt.gz"
xmapfile = mapdir + "hg19_cross_mappability_strength.txt.gz"

#resfilestring = "FHSobs_v_GTExWBpred.meqtl.trans.diffchr.allres.txt.gz"
resfile = resfilepath + resfilestring


genemap = dict()
for line in gzip.open(genemapfile):
    arr = line.strip().split()
    gene = arr[0][:15]
    if arr[1] != 'NA': #get rid of NAs
        mapscore = float(arr[1]) #convert values to float to use threshold
        if mapscore > 0.8:
            genemap[gene] = mapscore
    
xmap = dict()
for line in gzip.open(xmapfile):
    arr = line.strip().split()
    gene1 = arr[0][:15]
    gene2 = arr[1][:15]
    genes = str(gene1) + ':' + str(gene2)
    xmap[genes] = arr[2]

outres = gzip.open(rmdir + "Results_rm_bad_map_" + resfilestring,"wb")

for line in gzip.open(resfile):
    arr = line.strip().split()
    gene1 = arr[0][:15]
    gene2 = arr[1][:15]
    genepair = str(gene1) + ':' + str(gene2)
    if gene1 == 'snps':
        outres.write('\t'.join(arr) + '\n')
    elif gene1 in genemap and gene2 in genemap and genepair not in xmap:
        outres.write('\t'.join(arr) + '\n')

outres.close()
