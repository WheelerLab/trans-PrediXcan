#prep pred and obs gene lists to compare to eQTLGen
#using rm_suspect_NCBIgenes_Results_rm_bad_map_* genes

library(dplyr)
library(data.table)
library(tibble)
"%&%" = function(a,b) paste(a,b,sep="")

my.dir <- "/Volumes/im-lab/nas40t2/hwheeler/trans-px/"
#get tested (NCBI/bad map filtered) WB genes
wbfile <- my.dir %&% "NCBI_Gene_Summaries/rm_suspect_NCBIgenes_Results_rm_bad_map_FHSobs_v_GTExWBpred.meqtl.trans.diffchr.allres.txt.gz"
predwb <- fread("gzcat " %&% wbfile)

wbpredtested <- data.frame(pred=unique(predwb$snps))
wbobstested <- data.frame(obs=unique(predwb$gene))
#add chr annotation
annot <- fread(my.dir %&% "meqtl_format/gencode.v18.genes.patched_contigs.summary.protein") %>% mutate(short=substr(V5,1,15))

wbpred <- left_join(wbpredtested, annot, by = c("pred"="short")) %>% 
  dplyr::select(pred,chr=V1)
fwrite(wbpred, my.dir %&% "NCBI_Gene_Summaries/tested_genelists_eQTLGen_comp/FHS_GTExWB_pred_genelist.txt", sep=" ")

#fuma lists
fhstestedobs <- fread(my.dir %&% "NCBI_Gene_Summaries/FUMA_newFDR_2018-10-29/FHS_tested_obsgene_list.txt",header=F)
fhstestedpred <- fread(my.dir %&% "NCBI_Gene_Summaries/FUMA_newFDR_2018-10-29/FHS_tested_predgene_list.txt",header=F)
#add chr annotation
annot <- fread(my.dir %&% "meqtl_format/gencode.v18.genes.patched_contigs.summary.protein") %>% mutate(short=substr(V5,1,15))

#get shortened predname
wbpredgenes <- data.frame(fullpred=colnames(predwb[,-1:-2])) %>% mutate(shortpred=substr(fullpred,1,15))

#filtered pred list (in WB and rm suspect/bad map)
mtpred <- inner_join(fhstestedpred,annot,by=c("V1"="short")) %>% dplyr::select(pred=V1,chr=V1.y)
obsall <- inner_join(fhstestedobs,annot,by=c("V1"="short")) %>% dplyr::select(obs=V1,chr=V1.y)

fwrite(mtpred,file=my.dir %&% "NCBI_Gene_Summaries/tested_genelists_eQTLGen_comp/FHS_Multi_pred_genelist.txt",sep=' ')
fwrite(obsall,file=my.dir %&% "NCBI_Gene_Summaries/tested_genelists_eQTLGen_comp/FHS_obs_genelist.txt",sep=' ')
