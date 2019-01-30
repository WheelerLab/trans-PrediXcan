#!/usr/bin/python
'''
Retrive circular pairs and p-values for trans-px paper revision
'''
import gzip
mydir =  "/gpfs/data/im-lab/nas40t2/hwheeler/trans-px/NCBI_Gene_Summaries/"

wbfile = mydir + "rm_suspect_NCBIgenes_Results_rm_bad_map_FHSobs_v_GTExWBpred.meqtl.trans.diffchr.allres.txt.gz"
#wbfile = mydir + "test.gz"

#dictionary with (pred, obs) as key, pval as value.
resdict = dict()

#set of all pred genes (trans-acting)
predset = set()

#set of all obs genes (targets)
obsset = set()

for line in gzip.open(wbfile):
  (pred, obs, tstat, pval, fdr, beta) = line.strip().split()
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

outfile = gzip.open(mydir + "circular_pairs/circular_pairs_FHSobs_v_GTExWBpred.meqtl.trans.diffchr.allres.txt.gz","wb")
outfile.write("pair1\tpval1\tpair2\tpval2\n")


for pair in outset:
  print(pair[0] + '\t' + pair[1])
  outfile.write(pair[0] + '\t' + pair[1] + '\n')

outfile.close()
