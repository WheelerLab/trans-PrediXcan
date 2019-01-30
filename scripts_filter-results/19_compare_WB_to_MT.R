start.time <- Sys.time()
library(dplyr)
library(data.table)
library(tibble)
"%&%" = function(a,b) paste(a,b,sep="")
date <- Sys.Date()
args <- commandArgs(trailingOnly=T)
#args <- c("000","FHS") #for testing

listnum <- args[1]
pop <- args[2]

my.dir <- "/gpfs/data/im-lab/nas40t2/hwheeler/trans-px/"
#my.dir <- "/Volumes/im-lab/nas40t2/hwheeler/trans-px/"
predwbfile <- my.dir %&% "predFHSexp/GTEx-WB_predicted_expression.txt.gz"
predwb <- fread("zcat " %&% predwbfile)

obsfile <- my.dir %&% "obs" %&% pop %&% "exp/" %&% pop %&% "_obs_expression.txt.gz"

obs <- fread("zcat " %&% obsfile,header=TRUE)

#add those in predicted_expression files without obs expression
samples <- fread(my.dir %&% pop %&% "_dosages/samples.txt",header=FALSE) %>% mutate(IID=V2) %>% dplyr::select(IID)
obs <- left_join(samples,obs,by="IID")

obsgenetable <- fread(my.dir %&% "NCBI_Gene_Summaries/compare_WB_MT_models/genelists/FHS_tested_obsgene_chr_list.txt")
obsgenelist <- obsgenetable$obs
#obsgenelist <- obsgenelist[1:3] #for testing

pred.dir <- my.dir %&% "pred" %&% pop %&% "exp/" %&% pop %&% "_PredictedExpression_GTExdb/"
predgenetable <- fread(my.dir %&% "NCBI_Gene_Summaries/compare_WB_MT_models/genelists/genefilelist" %&% listnum)
#predgenefilelist <- predgenefilelist[1:3] #for testing


fwrite(list('obsgene','predgene','WB_adjRsq','WB_Pval','MT_Pval','MT_adjRsq','an_Fstat','an_Pval'),file=my.dir %&% 
         "NCBI_Gene_Summaries/compare_WB_MT_models/WB_v_MT_anova_" %&% 
         pop %&% "_list" %&% listnum %&% "_" %&% date %&% 
         ".txt",sep=" ")

for(j in 1:dim(obsgenetable)[[1]]){
  obsgene <- obsgenetable[j,1][[1]]
  if(obsgene %in% colnames(obs)){
    obsgenedf <- dplyr::select(obs,starts_with(obsgene))
    cat(j,"/",dim(obsgenetable)[[1]],"\n")
  }else{
    next
  }
  for(i in 1:dim(predgenetable)[[1]]){
    #skip if chr's match
    if(predgenetable[i,2][[1]] != obsgenetable[j,2][[1]]){
      predgene <- predgenetable[i,1][[1]]
      genefile <- predgene %&% ".txt"
      #pull GTEx-WB predicted exp
      predwbdf <- data.frame(dplyr::select(predwb,starts_with(predgene)))
      if(dim(predwbdf)[2]==1){
        fitwb <- lm(obsgenedf[,1]~predwbdf[,1])
        wbrsq <- summary(fitwb)$adj.r.squared
        wbp <- anova(fitwb)$P[1]
      }
      else{
        next
      }
      if(file.exists(pred.dir %&% "by_gene/PCs_by_gene/" %&% genefile)){
        genedf <- fread(pred.dir %&% "by_gene/PCs_by_gene/" %&% genefile, header=TRUE)
        genemat <- t(as.matrix(genedf[,-1]))
        colnames(genemat) <- genedf$tissue
        genemat <- genemat[,colSums(genemat)!=0,drop=FALSE] #rm any tissues where sum = 0
        ntis <- dim(genemat)[2] #num tissues in model
        if(ntis > 0){ #only run if sum of >=1 gene is > 0
          fitmt <- lm(obsgenedf[,1]~genemat)
          mtrsq <- summary(fitmt)$adj.r.squared
          mtp <- anova(fitmt)$P[1]
          anfit <- anova(fitwb,fitmt)
          anF <- anfit$F[2]
          anP <- anfit$`Pr(>F)`[2]
          gr <- list(obsgene,predgene,wbrsq,wbp,mtrsq,mtp,anF,anP)
          fwrite(gr,file=my.dir %&% "NCBI_Gene_Summaries/compare_WB_MT_models/WB_v_MT_anova_" %&% 
                   pop %&% "_list" %&% listnum %&% "_" %&% date %&% 
                   ".txt",sep=" ",append=TRUE,na="NA")
        }
      }
    }
  }
}
end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)


