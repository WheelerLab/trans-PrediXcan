start.time <- Sys.time()
library(dplyr)
library(data.table)
library(tibble)
"%&%" = function(a,b) paste(a,b,sep="")
date <- Sys.Date()
args <- commandArgs(trailingOnly=T)
#args <- c("277","GEU") #for testing

listnum <- args[1]
pop <- args[2]

#make test fam file
my.dir <- "/group/im-lab/nas40t2/hwheeler/trans-px/"
obsfile <- my.dir %&% "obs" %&% pop %&% "exp/" %&% pop %&% "_obs_expression.txt.gz"

obs <- fread("zcat " %&% obsfile,header=TRUE)

#add those in predicted_expression files without obs expression
samples <- fread(my.dir %&% pop %&% "_dosages/samples.txt",header=FALSE) %>% mutate(IID=V2) %>% dplyr::select(IID)
obs <- left_join(samples,obs,by="IID")

obsgenelist <- colnames(obs)[-1]
#obsgenelist <- obsgenelist[1:3] #for testing

pred.dir <- my.dir %&% "pred" %&% pop %&% "exp/" %&% pop %&% "_PredictedExpression_GTExdb/"
predgenefilelist <- scan(my.dir %&% "predFHSexp/FHS_PredictedExpression_GTExdb/by_gene/genefilelist" %&% listnum,"c")
#predgenefilelist <- predgenefilelist[1:3] #for testing

fwrite(list('obsgene','predgene','tissue','beta','se','tval','pval'),file=my.dir %&% 
         "multi-transpx_results/results_allgenes/results_PCs/multi-trans-px_" %&% 
         pop %&% "_tissue_results_list" %&% listnum %&% "_" %&% date %&% 
         ".txt",sep=" ")

fwrite(list('obsgene','predgene','Rsq','Fstat','Pval','ntis','ntislm'),file=my.dir %&% 
         "multi-transpx_results/results_allgenes/results_PCs/multi-trans-px_" %&% 
         pop %&% "_overall_results_list" %&% listnum %&% "_" %&% date %&% 
         ".txt",sep=" ")

for(j in 1:length(obsgenelist)){
  obsgene <- obsgenelist[j]
  obsgenedf <- dplyr::select(obs,starts_with(obsgene))
  cat(j,"/",length(obsgenelist),"\n")
  for(i in 1:length(predgenefilelist)){
    genefile <- predgenefilelist[i]
    predgene <- strsplit(genefile,".txt")[[1]][1]
    if(file.exists(pred.dir %&% "by_gene/PCs_by_gene/" %&% genefile)){
      genedf <- fread(pred.dir %&% "by_gene/PCs_by_gene/" %&% genefile, header=TRUE)
      genemat <- t(as.matrix(genedf[,-1]))
      colnames(genemat) <- genedf$tissue
      genemat <- genemat[,colSums(genemat)!=0,drop=FALSE] #rm any tissues where sum = 0
      ntis <- dim(genemat)[2] #num tissues in model
      if(ntis > 0){ #only run if sum of >=1 gene is > 0
        fit <- summary(lm(obsgenedf[,1]~genemat))
        rsq <- signif(fit$r.squared,5)
        an <- anova(lm(obsgenedf[,1]~genemat))
        pval <- signif(an$P[1],5)
        fstat <- signif(an$F[1],5)
        a <- signif(coef(fit),5)
        rows <- c("Intercept",colnames(genemat))[!fit$aliased] #rm any variables not defined because of singularities
        rownames(a) <- c(rows)
        ntislm <- length(rows) - 1 #num tissues in final model 
        er <- cbind(obsgene,predgene,rownames_to_column(data.frame(a),var="tissue"))
        gr <- list(obsgene,predgene,rsq,fstat,pval,ntis,ntislm)
        fwrite(gr,file=my.dir %&% "multi-transpx_results/results_allgenes/results_PCs/multi-trans-px_" %&% 
                 pop %&% "_overall_results_list" %&% listnum %&% "_" %&% date %&% 
                 ".txt",sep=" ",append=TRUE)
        fwrite(er,file=my.dir %&% "multi-transpx_results/results_allgenes/results_PCs/multi-trans-px_" %&% 
                 pop %&% "_tissue_results_list" %&% listnum %&% "_" %&% date %&% 
                 ".txt",sep=" ",append = TRUE)
      }
    }
  }
}

end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)


##setup: make predgene lists from output of 35_generate_predEXP_genefiles.R
##do this once

##from commandline:
# cd /group/im-lab/nas40t2/hwheeler/trans-px/predFHSexp/FHS_PredictedExpression_GTExdb/by_gene
# ls *txt > genefilelist
## limit to genes in meqtl_format/obsEXP_FHS_gtex.txt
# split -d -l 58 -a 3 genefilelist genefilelist

##makes 301 genefilelists from genefilelist000 - genefilelist300
