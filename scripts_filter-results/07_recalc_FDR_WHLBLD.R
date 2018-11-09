library(dplyr)
library(data.table)
library(ggplot2)
"%&%" = function(a,b) paste(a,b,sep="")
date <- Sys.Date()
my.dir <- '/gpfs/data/im-lab/nas40t2/hwheeler/trans-px/NCBI_Gene_Summaries/'

pop <- "FHS"

#plot all results
f <- fread('zcat ' %&% my.dir %&% 'rm_suspect_NCBIgenes_Results_rm_bad_map_FHSobs_v_GTExWBpred.meqtl.trans.diffchr.allres.txt.gz')
pdf(file = my.dir %&% pop %&% '_rm_bad_map_WHLBLD_hists.pdf')
hist(f$pvalue, main=pop %&% " WHLBLD trans-PrediXcan ALL", xlab="P-value",n=100)

#sort by P-value and get top 1e6 hits
f <- arrange(f,pvalue)
topf <- f[1:1e6,]

#calc B-H FDR based on number of tests:
nt <- dim(f)[1]
fres <- dplyr::mutate(topf, newFDR=p.adjust(pvalue,method="BH",n=nt),BONF=p.adjust(pvalue,method="bonferroni",n=nt))
#print the number of trans-px gene pairs with FDR < 0.05
fres_fdr0.05 <- dplyr::filter(fres,newFDR<0.05) %>% arrange(newFDR)
print(dim(fres_fdr0.05)[1])
#print the number of trans-px gene pairs with BONF P < 0.05
print(dim(dplyr::filter(fres_fdr0.05,BONF<0.05))[1])

dev.off()
fwrite(fres_fdr0.05,file=my.dir %&% 'rm_suspect_NCBIgenes_Results_rm_bad_map_FHSobs_v_GTExWBpred.meqtl.trans.diffchr.FDR0.05.txt',sep=" ")

