#!/usr/bin/python
'''
Retrive circular pairs and p-values for trans-px paper revision
'''
import gzip
mydir =  "/gpfs/data/im-lab/nas40t2/hwheeler/trans-px/"

#mtfile = mydir + "multi-transpx_results/results_allgenes/results_PCs/multi-trans-px_FHS_diff_chrs_overall_results_2017-12-11.txt.gz"
mtfile = mydir + "NCBI_Gene_Summaries/rm_suspect_NCBIgenes_Results_rm_bad_map_multi-trans-px_FHS_diff_chrs_overall_results_2017-12-11.txt.gz"

#dictionary with (pred, obs) as key, pval as value.
resdict = dict()

#set of all pred genes (trans-acting)
predset = set()

#set of all obs genes (targets)
obsset = set()

for line in gzip.open(mtfile):
  (obs, pred, Rsq, Fstat, pval, ntis, ntislm, predChr, predS1, predS2, predname, obsChr, obsS1, obsS2, obsname) = line.strip().split()
  pred = pred[:15] #rm .version number
  obs = obs[:15] #rm .version number
  resdict[(pred, obs)] = pval
  predset.add(pred)
  obsset.add(obs)

#build output set without duplicates
outset = set()

for p in predset:
  for o in obsset:
    if (p,o) in resdict and (o,p) in resdict:
      ab = p + '->' + o + '\t' + resdict[(p,o)]
      ba = o + '->' + p + '\t' + resdict[(o,p)]
      if (ba, ab) not in outset:
        outset.add((ab, ba))

outfile = gzip.open(mydir + "NCBI_Gene_Summaries/circular_pairs/circular_pairs_multi-trans-px_FHS_diff_chrs_overall_results_2017-12-11.txt.gz","wb")
outfile.write("pair1\tpval1\tpair2\tpval2\n")


for pair in outset:
  outfile.write(pair[0] + '\t' + pair[1] + '\n')

outfile.close()
