#read in diff chr results
library(dplyr)
library(data.table)
library(ggplot2)
"%&%" = function(a,b) paste(a,b,sep="")
date <- Sys.Date()
my.dir <- '/Volumes/im-lab/nas40t2/hwheeler/trans-px/NCBI_Gene_Summaries/'
####load data
#zcat rm_suspect_NCBIgenes_Results_rm_bad_map_FHSobs_v_GTExWBpred.meqtl.trans.diffchr.allres.txt.gz |wc
#24180873
#zcat rm_suspect_NCBIgenes_Results_rm_bad_map_multi-trans-px_FHS_diff_chrs_overall_results_2017-12-11.txt.gz |wc
#204944305

wbres <- fread("gzcat " %&% my.dir %&% "rm_suspect_NCBIgenes_Results_rm_bad_map_FHSobs_v_GTExWBpred.meqtl.trans.diffchr.allres.txt.gz")
mtres <- fread("gzcat " %&% my.dir %&% "rm_suspect_NCBIgenes_Results_rm_bad_map_multi-trans-px_FHS_overall_results_top1e6_2017-12-11.txt.gz")
header <- scan(my.dir %&% 'header','c')
colnames(mtres) <- header
mtres <- mutate(mtres,predgene=substr(predgene,1,15))
all <- full_join(wbres,mtres,by=c("snps"="predgene","gene"="obsgene"))
all <- mutate(all,gtexWB=ifelse(pvalue<1e-30,30,-log10(pvalue)),gtexMULTI=ifelse(Pval<1e-30,30,-log10(Pval)))


####qq plot compare GTEx-WB to all multi-tissue results
nnwb <- table(is.na(all$gtexWB))[[1]]
xxwb =  -log10((1:nnwb)/(nnwb+1))
wb <- dplyr::select(all,gtexWB) %>% dplyr::filter(!is.na(gtexWB)) %>% 
  mutate(model="Whole Blood", obsP=sort(gtexWB,decreasing=TRUE), expP=xxwb) %>%
  dplyr::select(-gtexWB)
nnmt <- 204944305 #see FHS_rm_bad_map_multi-trans-px_hists.pdf
xxmt =  -log10((1:nnmt)/(nnmt+1))
mt <- dplyr::select(all,gtexMULTI) %>% dplyr::filter(!is.na(gtexMULTI)) %>% 
  mutate(model="Multi-Tissue", obsP=sort(gtexMULTI,decreasing=TRUE,na.last=TRUE), expP=xxmt[1:dim(mtres)[1]]) %>% 
  dplyr::select(-gtexMULTI)
wbmt <- rbind(wb[1:1e6,],mt[1:1e6,])

qq1 <- ggplot(wbmt,aes(x=expP,y=obsP,col=model)) + geom_point(shape=1) + coord_cartesian(xlim=c(-0.05,8.05),ylim=c(-0.05,30.05)) + theme_bw(12) +
  scale_colour_manual(values=c(rgb(163,0,66,maxColorValue = 255),"dark gray")) + theme(legend.position = c(0.01,0.99),legend.justification = c(0,1)) +
  geom_abline(slope=1,intercept = 0) + xlab(expression(Expected~~-log[10](italic(p)))) + ylab(expression(Observed~~-log[10](italic(p))))

png(filename = my.dir %&% 'paper_figures_2018-10-26/Fig3_QQ-mt-PCs_v_wb-all.png',width=480,height=480)
qq1
dev.off()

tiff(filename=my.dir %&% 'paper_figures_2018-10-26/Fig3_QQ-mt-PCs_v_wb-all.tiff', width = 3.75, height = 3.75, units = 'in', res = 300, compression = 'lzw')
qq1
dev.off()

pdf(file = my.dir %&% 'paper_figures_2018-10-26/Fig3_QQ-mt-PCs_v_wb-all.pdf',width = 3.75, height = 3.75)
qq1
dev.off()

