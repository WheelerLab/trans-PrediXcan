#tests how often the gene pairs w/ trans-eQTLs occur by chance
#Angela Andaleon (aandaleon@luc.edu)
library(data.table)
library(dplyr)
set.seed(1000)
setwd("/home/angela/trans_px/chance_proportions/")
"%&%" = function(a, b) paste(a, b, sep = "")

db_SNPs_trans_eQTLGen <- fread("db_SNPs_trans_eQTLGen.csv")
cis_eQTLGen_FDR_05 <- fread("cis-eQTLGen_FDR_0.05_SNPs_only.txt", header = F)$V1
GTExWB_pred <- fread("FHS_GTExWB_pred_genelist.txt", stringsAsFactors = F)  
Multi_pred <- fread("FHS_Multi_pred_genelist.txt", stringsAsFactors = F)
obs <- fread("FHS_obs_genelist.txt", stringsAsFactors = F)

#take random samples of predicted and observed gene lists and see if they have trans-eQTLs and if those trans-eQTLs are also cis-eQTLs
num_samples_list <- c(55, 2356)
for(num_samples in num_samples_list){
  gene_pairs_w_trans_eQTL <- c() #store # of overlap from output
  trans_eQTL_also_cis_eQTL <- c()
  
  if(num_samples == 55){
    db_SNPs_trans_eQTLGen <- subset(db_SNPs_trans_eQTLGen, GTEx_tissue == "Whole_Blood") #restrict to just WB models
  }#else{
  #  db_SNPs_trans_eQTLGen$GTEx_tissue <- NA
  #  db_SNPs_trans_eQTLGen <- db_SNPs_trans_eQTLGen[!duplicated(db_SNPs_trans_eQTLGen), ] #collapse w/o tiss.
  #}
  
  for(i in 1:1000){
    pred_sample <- GTExWB_pred[sample(nrow(GTExWB_pred), num_samples),] #https://stackoverflow.com/questions/8273313/sample-random-rows-in-dataframe
    obs_row_nums <- sample(nrow(obs), num_samples) #without this extra step it throws: "Error: invalid first argument"
    obs_sample <- obs[obs_row_nums,]
    sample <- as.data.frame(cbind(pred_sample, obs_sample))
    colnames(sample) <- c("GTEx_gene", "GTEx_gene_chr", "eQTLGen_gene", "eQTLGen_gene_chr")
    sample <- sample %>% mutate_all(as.character) #why are you all factors
    sample <- subset(sample, GTEx_gene_chr != eQTLGen_gene_chr) #keep sample if genes are not on the same chr.
    repeat{ #keep adding randomized samples until its the proper number of pairs
      if(nrow(sample) < num_samples){ #if not a full sample set (if first sampling had cis-genes)
        pred_sample <- GTExWB_pred[sample(nrow(GTExWB_pred), num_samples - nrow(sample)), ]
        obs_row_nums <- sample(nrow(obs), num_samples - nrow(sample))
        obs_sample <- obs[obs_row_nums, ]
        sample2 <- as.data.frame(cbind(pred_sample, obs_sample))
        colnames(sample2) <- c("GTEx_gene", "GTEx_gene_chr", "eQTLGen_gene", "eQTLGen_gene_chr")
        sample2 <- sample2 %>% mutate_all(as.character)
        sample2 <- subset(sample2, GTEx_gene_chr != eQTLGen_gene_chr)
        sample <- rbind(sample, sample2)
      }else if(nrow(sample) == num_samples){ #once you have the correct number of gene pairs that are all on diff. chrs.
        break
      }
    }
    sample$GTEx_gene_chr <- NULL #we don't need these anymore
    sample$eQTLGen_gene_chr <- NULL
    overlap_samples <- dplyr::inner_join(sample, db_SNPs_trans_eQTLGen, by = c("GTEx_gene", "eQTLGen_gene")) #whats the overlap between randomized sample and trans-eQTLGen?
    if(nrow(overlap_samples) > 0){ #if there are overlaps, see if they are also cis-eQTLs
      gene_pairs_w_trans_eQTL_is_cis_eQTL <- subset(overlap_samples, rsid %in% cis_eQTLGen_FDR_05)
      trans_eQTL_also_cis_eQTL <- c(trans_eQTL_also_cis_eQTL, nrow(gene_pairs_w_trans_eQTL_is_cis_eQTL))
    }
    gene_pairs_w_trans_eQTL <- c(gene_pairs_w_trans_eQTL, nrow(overlap_samples)) #add number of overlaps to gene_pairs_w_trans_eQTL
  }
  
  if(num_samples == 55){ #whole blood
    pdf("GTExWB_chance_trans_eQTLGen_hist.pdf")
    num_gene_pairs_over_threshold <- length(which(gene_pairs_w_trans_eQTL >= 15))
    hist(gene_pairs_w_trans_eQTL, main = "# of randomized trans-gene pairs (GTExWB) with trans-eQTLs", xlab = num_gene_pairs_over_threshold %&% " pairs have 15 or more trans-eQTLs.\nThe average number of trans-eQTLs is " %&% round(mean(gene_pairs_w_trans_eQTL), 3) %&% ".\nOf all the observed trans-eQTLs, " %&% sum(trans_eQTL_also_cis_eQTL) %&% " out of " %&% sum(gene_pairs_w_trans_eQTL) %&% " are also cis-eQTLs.")
    abline(v = mean(gene_pairs_w_trans_eQTL), col = "blue")
    abline(v = 15, col = "red")
    dev.off()
  
    pdf("GTExWB_chance_trans_eQTLGen_boxplot.pdf")
    boxplot(gene_pairs_w_trans_eQTL, main = "# of randomized trans-gene pairs (GTExWB) with trans-eQTLs", xlab = num_gene_pairs_over_threshold %&% " pairs have 15 or more trans-eQTLs.\nThe average number of trans-eQTLs is " %&% round(mean(gene_pairs_w_trans_eQTL), 3) %&% ".\nOf all the observed trans-eQTLs, " %&% sum(trans_eQTL_also_cis_eQTL) %&% " out of " %&% sum(gene_pairs_w_trans_eQTL) %&% " are also cis-eQTLs.")
    abline(h = mean(gene_pairs_w_trans_eQTL), col = "blue")
    abline(h = 15, col = "red")
    dev.off()
  }
  
  if(num_samples == 2356){ #Multi
    pdf("Multi_chance_trans_eQTLGen_hist.pdf")
    num_gene_pairs_over_threshold <- length(which(gene_pairs_w_trans_eQTL >= 814))
    hist(gene_pairs_w_trans_eQTL, main = "# of randomized trans-gene pairs (Multi) with trans-eQTLs", xlab = num_gene_pairs_over_threshold %&% " pairs have 814 or more trans-eQTLs.\nThe average number of trans-eQTLs is " %&% round(mean(gene_pairs_w_trans_eQTL), 3) %&% ".\nOf all the observed trans-eQTLs, " %&% sum(trans_eQTL_also_cis_eQTL) %&% " out of " %&% sum(gene_pairs_w_trans_eQTL) %&% " are also cis-eQTLs.")
    abline(v = mean(gene_pairs_w_trans_eQTL), col = "blue")
    abline(v = 814, col = "red")
    dev.off()
    
    pdf("Multi_chance_trans_eQTLGen_boxplot.pdf")
    boxplot(gene_pairs_w_trans_eQTL, main = "# of randomized trans-gene pairs (Multi) with trans-eQTLs", xlab = num_gene_pairs_over_threshold %&% " pairs have 814 or more trans-eQTLs.\nThe average number of trans-eQTLs is " %&% round(mean(gene_pairs_w_trans_eQTL), 3) %&% ".\nOf all the observed trans-eQTLs, " %&% sum(trans_eQTL_also_cis_eQTL) %&% " out of " %&% sum(gene_pairs_w_trans_eQTL) %&% " are also cis-eQTLs.")
    abline(h = mean(gene_pairs_w_trans_eQTL), col = "blue")
    abline(h = 814, col = "red")
    dev.off()
  }
}
