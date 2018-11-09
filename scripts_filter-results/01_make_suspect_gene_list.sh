#!/bin/bash

#remove genes with NCBI gene summaries that include the following terms:
# retro
# pseudogene
# paralog


zcat hgFixed.refSeqSummary_gencode.v19.annotations.txt.gz |grep -v pseudogene | grep -v retro |grep -v paralog > hgFixed.refSeqSummary_gencode.v19.annotations_rm_suspect.txt
gzip hgFixed.refSeqSummary_gencode.v19.annotations_rm_suspect.txt

zcat hgFixed.refSeqSummary_gencode.v19.annotations.txt.gz |grep pseudogene > hgFixed.refSeqSummary_gencode.v19.annotations_suspect_list.txt
zcat hgFixed.refSeqSummary_gencode.v19.annotations.txt.gz |grep -v pseudogene > complement
grep retro complement >> hgFixed.refSeqSummary_gencode.v19.annotations_suspect_list.txt
grep -v retro complement >> complement2
grep paralog complement2 >> hgFixed.refSeqSummary_gencode.v19.annotations_suspect_list.txt
gzip hgFixed.refSeqSummary_gencode.v19.annotations_suspect_list.txt
rm complement*
