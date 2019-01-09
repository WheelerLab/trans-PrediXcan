# -*- coding: utf-8 -*-
"""
Created on Tue Jan 08 14:47:20 2019
#pull all SNPs from predicted expression models and check if they're trans-eQTLs in eQTLGen
@author: Angela Andaleon (aandaleon@luc.edu)
"""

import numpy as np
import sqlite3
import pandas as pd
import re
pd.options.mode.chained_assignment = None

tissues = ["Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Adrenal_Gland", "Artery_Aorta", "Artery_Coronary", "Artery_Tibial", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Caudate_basal_ganglia", "Brain_Cerebellar_Hemisphere", "Brain_Cortex", "Brain_Frontal_Cortex_BA9", "Brain_Hippocampus", "Brain_Hypothalamus", "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia", "Breast_Mammary_Tissue", "Cells_EBV-transformed_lymphocytes", "Cells_Transformed_fibroblasts", "Colon_Sigmoid", "Colon_Transverse", "Esophagus_Gastroesophageal_Junction", "Esophagus_Mucosa", "Esophagus_Muscularis", "Heart_Atrial_Appendage", "Heart_Left_Ventricle", "Liver", "Lung", "Muscle_Skeletal", "Nerve_Tibial", "Ovary", "Pancreas", "Pituitary", "Prostate", "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg", "Small_Intestine_Terminal_Ileum", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole_Blood"]
#tissues = ["Whole_Blood"]

#pull all SNPs from predicted expression models
db_SNPs = [] #store SNPs
for tissue in tissues:
  conn = sqlite3.connect("/home/wheelerlab3/Data/PrediXcan_db/GTEx-V6p-HapMap-2016-09-08/TW_" + tissue + "_0.5.db") #connect to .db file
  c = conn.cursor()
  c.execute('select * from weights;') #run search throguh python
  for row in c:
     (rsid, gene, weight, ref_allele, eff_allele) = row[0:5]
     db_SNPs.append([rsid, gene, tissue])
  conn.close()
db_SNPs = pd.DataFrame(db_SNPs, columns = ['rsid', 'gene', 'tissue'])
db_SNPs.columns = ["rsid", "GTEx_gene", "GTEx_tissue"]
db_SNPs['GTEx_gene'] = db_SNPs['GTEx_gene'].apply(lambda x: re.sub("\\..*", "", x)) #remove suffix to genes https://stackoverflow.com/questions/32902837/replace-substring-in-pandas-data-frame-column

#read in trans-eQTLs
trans_eQTLGen_FDR_05 = pd.read_csv("eQTLGen_FDR_0.05.txt", delim_whitespace = True, header = None)
trans_eQTLGen_FDR_05.columns = ["eQTLGen_P", "rsid", "rsid_chr", "rsid_pos", "Zscore", "A1", "A0", "eQTLGen_gene", "eQTLGen_gene_symbol", "eQTLGen_gene_chr", "eQTLGen_gene_pos", "eQTLGen_nr_cohorts", "eQTLGen_nr_samples", "eQTLGen_FDR"]

#Are SNPs in pred exp models trans-eQTLs?
db_SNPs_trans_eQTLGen = pd.merge(db_SNPs, trans_eQTLGen_FDR_05, how = "outer", on = "rsid")
db_SNPs_trans_eQTLGen_noNA = db_SNPs_trans_eQTLGen.dropna() #drop SNPs that didn't overlap
db_SNPs_trans_eQTLGen_noNA = db_SNPs_trans_eQTLGen_noNA[["rsid", "rsid_chr", "rsid_pos", "GTEx_gene", "GTEx_tissue", "eQTLGen_gene_chr", "eQTLGen_gene_pos", "eQTLGen_gene", "eQTLGen_gene_symbol", "eQTLGen_P", "eQTLGen_FDR"]]
db_SNPs_trans_eQTLGen_noNA.to_csv("db_SNPs_trans_eQTLGen.csv", index = False)
print("Completed writing .db SNP and trans-eQTL overlap to db_SNPs_trans_eQTLGen.csv.")
db_SNPs_trans_eQTLGen_obsgene_predgene = db_SNPs_trans_eQTLGen_noNA[["GTEx_gene", "eQTLGen_gene", "rsid"]] #these genes have a cis-eQTL in a .db gene that acts as a trans-eQTL to another gene
db_SNPs_trans_eQTLGen_obsgene_predgene = db_SNPs_trans_eQTLGen_obsgene_predgene.rename(index = str, columns = {"GTEx_gene": "predgene", "eQTLGen_gene": "obsgene", "rsid": "rsid"})
db_SNPs_trans_eQTLGen_obsgene_predgene["trans_eQTL_in_eQTLGen"] = "Yes"

#list of all cis-eQTLs in eQTLGen w/ FDR < 0.05
cis_eQTL = set(np.loadtxt("cis-eQTLGen_FDR_0.05_SNPs_only.txt", dtype = str))
db_SNPs_trans_eQTLGen_obsgene_predgene["trans_eQTL_also_cis_eQTL"] = "No"
db_SNPs_trans_eQTLGen_obsgene_predgene.loc[db_SNPs_trans_eQTLGen_obsgene_predgene.rsid.isin(cis_eQTL), ['trans_eQTL_also_cis_eQTL']] = "Yes" #Add another column in TableS1 and TableS2 - is the trans-eQTL you found (keep this column) also a cis-eQTL?

#subset to WHLBLD results
WHLBLD = pd.read_csv("TableS1_WHLBLD_results_2018-10-29.csv")
WHLBLD.columns = ['predgene', 'predname', 'predChr', 'predS1', 'predS2', 'obsgene', 'obsname', 'obsChr', 'obsS1', 'obsS2', 'FHS_stat', 'FHS_beta', 'FHS_pval', 'FHS_FDR', 'DGN_stat', 'DGN_beta', 'DGN_pval']
WHLBLD_trans_eQTL = pd.merge(WHLBLD, db_SNPs_trans_eQTLGen_obsgene_predgene, how = "left", on = ["predgene", "obsgene"]).drop_duplicates() #force trans-eQTLGen results to order of WHLBLD
WHLBLD_trans_eQTL["trans_eQTL_in_eQTLGen"].fillna("No", inplace = True) #if there's no matching in trans-eQTLGen
WHLBLD_trans_eQTL.to_csv("WHLBLD_trans_eQTL.csv", index = False, na_rep = "NA")
print("Completed writing WHLBLD trans-PX and trans-eQTLGen overlapping results to WHLBLD_trans_eQTL.csv")

#subset to MultiXcan results
MultiXcan = pd.read_csv("TableS2_MultiXcan_results_2018-10-29.csv")
MultiXcan.columns = ['predgene', 'predname', 'predChr', 'predS1', 'predS2', 'obsgene', 'obsname', 'obsChr', 'obsS1', 'obsS2', 'FHS_stat', 'FHS_beta', 'FHS_pval', 'FHS_FDR', 'DGN_stat', 'DGN_beta', 'DGN_pval']
MultiXcan_trans_eQTL = pd.merge(MultiXcan, db_SNPs_trans_eQTLGen_obsgene_predgene, how = "left", on = ["predgene", "obsgene"]).drop_duplicates() #force trans-eQTLGen results to order of WHLBLD
MultiXcan_trans_eQTL["trans_eQTL_in_eQTLGen"].fillna("No", inplace = True) #if there's no matching in trans-eQTLGen
MultiXcan_trans_eQTL.to_csv("MultiXcan_trans_eQTL.csv", index = False, na_rep = "NA")
print("Completed writing MultiXcan trans-PX and trans-eQTLGen overlapping results to MultiXcan_trans_eQTL.csv")

#Recreate Table S3 - Trans-acting/target gene pairs from FHS MulTiXcan (FDR < 0.05) that were previously identified as trans-eQTLs (FDR < 0.05) in whole blood in Westra et al. (https://molgenis58.target.rug.nl/bloodeqtlbrowser/)
#\textit{trans}-eQTLs within 1Mb of our PrediXcan and MulTiXcan trans-acting genes with the same target genes are shown in Table S3. 
#db_SNPs = pd.read_csv("db_SNPs.txt")
#okay so WHLBLD has predgene and obsgene pairs
#and trans_eQTLGen_FDR_05 has SNP and obsgene pairs
#match on: Westra.GENEname == obsname; Westra.SNPchr == predChr
    #and then restrict Westra.SNPpos to - 1Mb of predS1 and + 1 Mb of predS2
trans_eQTLGen_FDR_05.columns = ['eQTLGen_P', 'rsid', 'predChr', 'rsid_pos', 'Zscore', 'A1', 'A0', 'obsgene', 'eQTLGen_gene_symbol', 'eQTLGen_gene_chr', 'eQTLGen_gene_pos', 'eQTLGen_nr_cohorts', 'eQTLGen_nr_samples', 'eQTLGen_FDR']
WHLBLD_match_eQTLGen = pd.merge(WHLBLD, trans_eQTLGen_FDR_05, on = ["predChr", "obsgene"])
WHLBLD_match_eQTLGen.loc[WHLBLD_match_eQTLGen.rsid_pos <= (WHLBLD_match_eQTLGen.predS1 - 1000000), ['rsid_pos']] = None #keep if w/in +/- 1 Mb https://stackoverflow.com/questions/19226488/change-one-value-based-on-another-value-in-pandas
WHLBLD_match_eQTLGen.loc[WHLBLD_match_eQTLGen.rsid_pos >= (WHLBLD_match_eQTLGen.predS2 + 1000000), ['rsid_pos']] = None
WHLBLD_match_eQTLGen = WHLBLD_match_eQTLGen[np.isfinite(WHLBLD_match_eQTLGen['rsid_pos'])] #keep SNPs that pass
WHLBLD_match_eQTLGen = WHLBLD_match_eQTLGen[["eQTLGen_P", "rsid", "predChr", "rsid_pos", "eQTLGen_gene_chr", "eQTLGen_gene_pos", "eQTLGen_gene_symbol", "eQTLGen_FDR", "predgene", "predname", "predS1", "obsgene", "obsname", "obsS1", "obsS2", "FHS_pval", "FHS_FDR", "DGN_pval"]] #keep important cols.
WHLBLD_match_eQTLGen.to_csv("WHLBLD_match_eQTLGen.csv", index = False)

MultiXcan_match_eQTLGen = pd.merge(MultiXcan, trans_eQTLGen_FDR_05, on = ["predChr", "obsgene"])
MultiXcan_match_eQTLGen.loc[MultiXcan_match_eQTLGen.rsid_pos <= (MultiXcan_match_eQTLGen.predS1 - 1000000), ['rsid_pos']] = None #keep if w/in +/- 1 Mb https://stackoverflow.com/questions/19226488/change-one-value-based-on-another-value-in-pandas
MultiXcan_match_eQTLGen.loc[MultiXcan_match_eQTLGen.rsid_pos >= (MultiXcan_match_eQTLGen.predS2 + 1000000), ['rsid_pos']] = None
MultiXcan_match_eQTLGen = MultiXcan_match_eQTLGen[np.isfinite(MultiXcan_match_eQTLGen['rsid_pos'])] #keep SNPs that pass
MultiXcan_match_eQTLGen = MultiXcan_match_eQTLGen[["eQTLGen_P", "rsid", "predChr", "rsid_pos", "eQTLGen_gene_chr", "eQTLGen_gene_pos", "eQTLGen_gene_symbol", "eQTLGen_FDR", "predgene", "predname", "predS1", "obsgene", "obsname", "obsS1", "obsS2", "FHS_pval", "FHS_FDR", "DGN_pval"]]
MultiXcan_match_eQTLGen.to_csv("MultiXcan_match_eQTLGen.csv", index = False)
print("Completed comparing trans-acting/target gene pairs b/w trans-PX and trans-eQTLGen.")

