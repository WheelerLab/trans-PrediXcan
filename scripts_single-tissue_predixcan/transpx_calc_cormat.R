library(dplyr)
library(data.table)
library(dtplyr)
"%&%" = function(a,b) paste(a,b,sep="")
args <- commandArgs(trailingOnly=T)
date <- Sys.Date()

#calc pred v obs corr matrix, combine top hits (P<0.05) with MatrixEQTL hack output
#doesn't output the corr matrix--faster to regenerate when needed than write/read

filedir <- args[1] #file directory
predfile <- filedir %&% args[2]
obsfile <- filedir %&% args[3]
annotfile <- filedir %&% args[4]
meqtlresfile <- args[5] #assumes file is in working dir
outprefix <- args[6] #ouputs to working dir, e.g. GEUobs_v_DGNpred

#files for testing
#filedir <- "/Volumes/im-lab/nas40t2/hwheeler/trans-px/meqtl_format/"
#obsfile <- filedir %&% "obsEXP_GEU_from_DGNdb.txt"
#predfile <- filedir %&% "predEXP_GEU_from_DGNdb.txt"
#annotfile <- filedir %&% "gene_annot_GEU_from_DGNdb.txt"
#meqtlresfile <- "GEUobs_v_DGNpred.meqtl.trans.txt"
#outprefix <- "testcormat"

annot <- fread(annotfile) %>% dplyr::select(gene_id,strand,chrnum,s1,s2,ensembl_id)

obs <- fread(obsfile,header=T)
pred <- fread(predfile,header=T)
obsmat <- t(as.matrix(obs[,-1]))
predmat <- t(as.matrix(pred[,-1]))
colnames(obsmat) <- obs$gene
colnames(predmat) <- pred$gene

n <- dim(obsmat)[1]

obsmat <- scale(obsmat,center=TRUE,scale=TRUE)/sqrt(n-1)
predmat <- scale(predmat,center=TRUE,scale=TRUE)/sqrt(n-1)

cor2<- t(predmat) %*% obsmat

#combine with MatrixEQTL output
long <- as.data.frame.table(cor2)
colnames(long) <- c('pred_ensembl', 'obs_ensembl', 'R')
#mres <- fread(meqtlresfile)
mres <- fread('zcat ' %&% meqtlresfile)
combine <- right_join(long,mres,by=c("pred_ensembl"="snps","obs_ensembl"="gene"))

#get position info
a<-left_join(combine,annot,by=c('pred_ensembl'='ensembl_id')) %>% 
  rename(pred_gene_id=gene_id,pred_strand=strand,pred_chr=chrnum,pred_s1=s1,pred_s2=s2)
combine<-left_join(a,annot,by=c('obs_ensembl'='ensembl_id')) %>% 
  rename(obs_gene_id=gene_id,obs_strand=strand,obs_chr=chrnum,obs_s1=s1,obs_s2=s2)

#output top results
output <- dplyr::filter(combine,pvalue<0.05)
write.table(output, "trans_" %&% outprefix %&% "_pval0.05_" %&% date %&% ".txt",quote=F,row.names = F,sep = '\t')

outputdiff <- dplyr::filter(output,pred_chr != obs_chr)
write.table(outputdiff, "trans-diffChrs_" %&% outprefix %&% "_pval0.05_" %&% date %&% ".txt",quote=F,row.names = F,sep = '\t')
