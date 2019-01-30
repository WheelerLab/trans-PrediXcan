#prep pred and obs gene lists to compare WB to MT
#using rm_suspect_NCBIgenes_Results_rm_bad_map_* genes

library(dplyr)
library(data.table)
library(tibble)
"%&%" = function(a,b) paste(a,b,sep="")

my.dir <- "/Volumes/im-lab/nas40t2/hwheeler/trans-px/"
predwbfile <- my.dir %&% "predFHSexp/GTEx-WB_predicted_expression.txt.gz"
predwb <- fread("gzcat " %&% predwbfile)

#fuma lists
fhstestedobs <- fread(my.dir %&% "NCBI_Gene_Summaries/FUMA_newFDR_2018-10-29/FHS_tested_obsgene_list.txt",header=F)
fhstestedpred <- fread(my.dir %&% "NCBI_Gene_Summaries/FUMA_newFDR_2018-10-29/FHS_tested_predgene_list.txt",header=F)
#add chr annotation
annot <- fread(my.dir %&% "meqtl_format/gencode.v18.genes.patched_contigs.summary.protein") %>% mutate(short=substr(V5,1,15))

#get shortened predname
wbpredgenes <- data.frame(fullpred=colnames(predwb[,-1:-2])) %>% mutate(shortpred=substr(fullpred,1,15))

#filtered pred list (in WB and rm suspect/bad map)
predall <- inner_join(wbpredgenes, fhstestedpred, by=c("shortpred"="V1")) 
predall <- inner_join(predall,annot,by=c("shortpred"="short")) 
obsall <- inner_join(fhstestedpred,annot,by=c("V1"="short")) #SPOTTED mistake on 2019-01-17 in case I need to rerun-switch to fhstestedobs

predlist <- dplyr::select(predall, fullpred, V1)
obslist <- dplyr::select(obsall,V1,V1.y)
colnames(predlist) <- c("pred","chr")
colnames(obslist) <- c("obs","chr")

fwrite(obslist,file=my.dir %&% "NCBI_Gene_Summaries/compare_WB_MT_models/genelists/FHS_tested_obsgene_chr_list.txt",sep=' ')
fwrite(predlist,file=my.dir %&% "NCBI_Gene_Summaries/compare_WB_MT_models/genelists/FHS_tested_predgene_chr_list.txt",sep=' ')


##from commandline:
# cd /gpfs/data/im-lab/nas40t2/hwheeler/trans-px/NCBI_Gene_Summaries/compare_WB_MT_models/genelists/

##makes 101 genefilelists from genefilelist000 - genefilelist100
# split -d -l 58 -a 3 FHS_tested_predgene_chr_list.txt genefilelist