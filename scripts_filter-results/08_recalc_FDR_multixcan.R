library(dplyr)
library(data.table)
library(ggplot2)
"%&%" = function(a,b) paste(a,b,sep="")
date <- Sys.Date()
my.dir <- '/gpfs/data/im-lab/nas40t2/hwheeler/trans-px/NCBI_Gene_Summaries/'

pop <- "FHS"

#plot all results
f <- fread('zcat ' %&% my.dir %&% 'rm_suspect_NCBIgenes_Results_rm_bad_map_multi-trans-px_FHS_diff_chrs_overall_results_2017-12-11.txt.gz',header=F)
colnames(f) <- c("obsgene", "predgene", "Rsq", "Fstat", "Pval", "ntis", "ntislm", "predChr", "predS1", "predS2", "predname", "obsChr", "obsS1", "obsS2", "obsname")
pdf(file = my.dir %&% pop %&% '_rm_bad_map_multi-trans-px_hists.pdf')
hist(f$Pval, main=pop %&% " multi-tissue trans-PrediXcan ALL", xlab="P-value",n=100)
hist(f$Rsq, main=pop %&% " multi-tissue trans-PrediXcan ALL", xlab="Rsq", n=100, ylim=c(0,10000))
hist(f$ntislm, main=pop %&% " multi-tissue trans-PrediXcan ALL", xlab="Number of terms in gene model",n=50)

#sort by P-value and print top 1e6 hits
f <- arrange(f,Pval)
topf <- f[1:1e6,]
fwrite(topf, my.dir %&% 'rm_suspect_NCBIgenes_Results_rm_bad_map_multi-trans-px_' %&% pop %&% '_overall_results_top1e6_2017-12-11.txt', sep=" " )

#calc B-H FDR based on number of tests:
nt <- dim(f)[1]
#zcat Results_rm_bad_map_multi-trans-px_FHS_diff_chrs_overall_results_2017-12-11.txt.gz |wc
fres <- dplyr::mutate(topf, FDR=p.adjust(Pval,method="BH",n=nt),BONF=p.adjust(Pval,method="bonferroni",n=nt))
#print the number of trans-px gene pairs with FDR < 0.05
fres_fdr0.05 <- dplyr::filter(fres,FDR<0.05) %>% arrange(FDR)
print(dim(fres_fdr0.05)[1])
#print the number of trans-px gene pairs with BONF P < 0.05
print(dim(dplyr::filter(fres_fdr0.05,BONF<0.05))[1])

hist(fres_fdr0.05$Rsq, main=pop %&% " multi-tissue trans-PrediXcan FDR < 0.05\n n tests = " %&% nt, xlab="Rsq", n=100)
hist(fres_fdr0.05$ntislm, main=pop %&% " multi-tissue trans-PrediXcan FDR < 0.05\n n = " %&% dim(fres_fdr0.05)[1], xlab="Number of terms in gene model",n=50)
dev.off()
fwrite(fres_fdr0.05,file=my.dir %&% 'rm_suspect_NCBIgenes_Results_rm_bad_map_multi-trans-px_' %&% pop %&% '_overall_results_FDR0.05_2017-12-11.txt',sep=" ")

