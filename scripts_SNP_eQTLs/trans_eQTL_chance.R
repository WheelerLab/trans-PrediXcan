#tests how often the gene pairs w/ trans-eQTLs occur by chance
#Angela Andaleon (aandaleon@luc.edu)
library(data.table)
library(dplyr)
set.seed(1337)
setwd("/home/angela/trans_px/chance_proportions/")

db_SNPs_trans_eQTLGen <- fread("db_SNPs_trans_eQTLGen.csv")
GTExWB_pred <- fread("FHS_GTExWB_pred_genelist.txt", stringsAsFactors = F)  
Multi_pred <- fread("FHS_Multi_pred_genelist.txt", stringsAsFactors = F)
obs <- fread("FHS_obs_genelist.txt", stringsAsFactors = F)

#take random samples of predicted and observed gene lists and see if they have trans-eQTLs
gene_pairs_w_trans_eQTL <- c() #store # of overlap from output
for(i in 1:1000){
  pred_sample <- sample(GTExWB_pred$pred, 55)
  obs_sample <- sample(obs$obs, 55)
  sample <- as.data.frame(cbind(pred_sample, obs_sample))
  sample <- sample %>% mutate_all(as.character) #why are you all factors
  colnames(sample) <- c("GTEx_gene", "eQTLGen_gene")
  overlap_samples <- dplyr::inner_join(sample, db_SNPs_trans_eQTLGen, by = c("GTEx_gene", "eQTLGen_gene")) #whats the overlap between randomized sample and trans-eQTLGen?
  gene_pairs_w_trans_eQTL <- c(gene_pairs_w_trans_eQTL, nrow(overlap_samples)) #add number of overlaps to gene_pairs_w_trans_eQTL
}

pdf("GTExWB_chance_trans_eQTLGen_hist.pdf")
hist(gene_pairs_w_trans_eQTL, main = "# of randomized gene pairs with trans-eQTLs")
abline(v = mean(gene_pairs_w_trans_eQTL), col = "blue")
abline(v = 15, col = "red")
dev.off()

pdf("GTExWB_chance_trans_eQTLGen_boxplot.pdf")
boxplot(gene_pairs_w_trans_eQTL, main = "# of randomized gene pairs with trans-eQTLs")
abline(h = mean(gene_pairs_w_trans_eQTL), col = "blue")
abline(h = 15, col = "red")
dev.off()

#take random samples of predicted and observed gene lists and see if they have trans-eQTLs
gene_pairs_w_trans_eQTL <- c() #store # of overlap from output
for(i in 1:1000){
  pred_sample <- sample(Multi_pred$pred, 2356)
  obs_sample <- sample(obs$obs, 2356)
  sample <- as.data.frame(cbind(pred_sample, obs_sample))
  sample <- sample %>% mutate_all(as.character) #why are you all factors
  colnames(sample) <- c("GTEx_gene", "eQTLGen_gene")
  overlap_samples <- dplyr::inner_join(sample, db_SNPs_trans_eQTLGen, by = c("GTEx_gene", "eQTLGen_gene")) #whats the overlap between randomized sample and trans-eQTLGen?
  gene_pairs_w_trans_eQTL <- c(gene_pairs_w_trans_eQTL, nrow(overlap_samples)) #add number of overlaps to gene_pairs_w_trans_eQTL
}

pdf("Multi_chance_trans_eQTLGen_hist.pdf")
hist(gene_pairs_w_trans_eQTL, main = "# of randomized gene pairs with trans-eQTLs")
abline(v = mean(gene_pairs_w_trans_eQTL), col = "blue")
abline(v = 814, col = "red")
dev.off()

pdf("Multi_chance_trans_eQTLGen_boxplot.pdf")
boxplot(gene_pairs_w_trans_eQTL, main = "# of randomized gene pairs with trans-eQTLs")
abline(h = mean(gene_pairs_w_trans_eQTL), col = "blue")
abline(h = 814, col = "red")
dev.off()
