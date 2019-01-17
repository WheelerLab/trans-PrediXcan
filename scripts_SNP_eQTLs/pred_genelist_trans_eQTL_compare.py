#do the SNPs in these lists have a trans-eQTL?
#Angela Andaleon (aandaleon@luc.edu)

#import numpy as np
import pandas as pd
import re
pd.options.mode.chained_assignment = None

#read in trans-eQTLs and .db SNP overlaps
db_SNPs_trans_eQTLGen = pd.read_csv("db_SNPs_trans_eQTLGen.csv")
db_SNPs_trans_eQTLGen = db_SNPs_trans_eQTLGen[["rsid_chr", "rsid", "GTEx_gene", "eQTLGen_gene_chr", "eQTLGen_gene", "eQTLGen_P", "eQTLGen_FDR"]]

#restrict to just the observed genes
obs_genelist = pd.read_csv("FHS_obs_genelist.txt", delim_whitespace = True) #obs = eQTLGen
db_SNPs_trans_eQTLGen_obs = db_SNPs_trans_eQTLGen.loc[db_SNPs_trans_eQTLGen['eQTLGen_gene'].isin(obs_genelist.obs.tolist())]

#match to input list
GTExWB_genelist = pd.read_csv("FHS_GTExWB_pred_genelist.txt", delim_whitespace = True)
GTExWB_genelist["chr"] = GTExWB_genelist["chr"].apply(lambda x: re.sub("chr*", "", x)) #why do these genelist chr have "chr"s in them? Oh well, remove
db_SNPs_trans_eQTLGen_GTExWB = db_SNPs_trans_eQTLGen_obs.loc[db_SNPs_trans_eQTLGen_obs['GTEx_gene'].isin(GTExWB_genelist.pred.tolist())]
db_SNPs_trans_eQTLGen_GTExWB = db_SNPs_trans_eQTLGen_GTExWB.drop(labels = "rsid", axis = 1).drop_duplicates() #drop SNP column and collapse duplicates
db_SNPs_trans_eQTLGen_GTExWB.to_csv("db_SNPs_trans_eQTLGen_GTExWB.csv", index = False)

Multi_genelist = pd.read_csv("FHS_Multi_pred_genelist.txt", delim_whitespace = True)
Multi_genelist["chr"] = Multi_genelist["chr"].apply(lambda x: re.sub("chr*", "", x)) #why do these genelist chr have "chr"s in them? Oh well, remove
db_SNPs_trans_eQTLGen_Multi = db_SNPs_trans_eQTLGen_obs.loc[db_SNPs_trans_eQTLGen_obs['GTEx_gene'].isin(Multi_genelist.pred.tolist())]
db_SNPs_trans_eQTLGen_Multi = db_SNPs_trans_eQTLGen_Multi.drop(labels = "rsid", axis = 1).drop_duplicates() #drop SNP column and collapse duplicates
db_SNPs_trans_eQTLGen_Multi.to_csv("db_SNPs_trans_eQTLGen_Multi.csv", index = False)
