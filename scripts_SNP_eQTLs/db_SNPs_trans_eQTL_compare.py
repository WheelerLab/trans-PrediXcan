'''
Created on Tue Jan 08 14:47:20 2019
#pull all SNPs from predicted expression models and check if they're cis- or trans-eQTLs in eQTLGen
@author: Angela Andaleon (aandaleon@luc.edu)
'''

import sqlite3
import pandas as pd
import re
pd.options.mode.chained_assignment = None

tissues = ["Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Adrenal_Gland", "Artery_Aorta", "Artery_Coronary", "Artery_Tibial", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Caudate_basal_ganglia", "Brain_Cerebellar_Hemisphere", "Brain_Cortex", "Brain_Frontal_Cortex_BA9", "Brain_Hippocampus", "Brain_Hypothalamus", "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia", "Breast_Mammary_Tissue", "Cells_EBV-transformed_lymphocytes", "Cells_Transformed_fibroblasts", "Colon_Sigmoid", "Colon_Transverse", "Esophagus_Gastroesophageal_Junction", "Esophagus_Mucosa", "Esophagus_Muscularis", "Heart_Atrial_Appendage", "Heart_Left_Ventricle", "Liver", "Lung", "Muscle_Skeletal", "Nerve_Tibial", "Ovary", "Pancreas", "Pituitary", "Prostate", "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg", "Small_Intestine_Terminal_Ileum", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole_Blood"]

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
#db_SNPs_trans_eQTLGen_noNA = pd.read_csv("db_SNPs_trans_eQTLGen.csv")
print("Completed writing .db SNP and trans-eQTL overlap to db_SNPs_trans_eQTLGen.csv.")
db_SNPs_trans_eQTLGen_noNA = db_SNPs_trans_eQTLGen_noNA.rename(columns = {"GTEx_gene":"predgene", "eQTLGen_gene":"obsgene", "eQTLGen_P":"eQTLGen_trans_P", "eQTLGen_FDR":"eQTLGen_trans_FDR"})

#list of all cis-eQTLs in eQTLGen w/ FDR < 0.05
cis_eQTL = pd.read_csv("cis-eQTLGen_P_FDR.txt", delim_whitespace = True, header = None)
cis_eQTL.columns = ["eQTLGen_cis_P", "rsid", "predgene", "eQTLGen_cis_FDR"] #predgene = eQTLGen_cis_gene

#subset to WHLBLD results
db_SNPs_trans_eQTLGen_WHLBLD = db_SNPs_trans_eQTLGen_noNA.query('GTEx_tissue == "Whole_Blood"')
WHLBLD = pd.read_csv("TableS1_WHLBLD_results_2018-10-29.csv")
WHLBLD.columns = ['predgene', 'predname', 'predChr', 'predS1', 'predS2', 'obsgene', 'obsname', 'obsChr', 'obsS1', 'obsS2', 'FHS_stat', 'FHS_beta', 'FHS_pval', 'FHS_FDR', 'DGN_stat', 'DGN_beta', 'DGN_pval']
WHLBLD_trans_eQTL = pd.merge(WHLBLD, db_SNPs_trans_eQTLGen_WHLBLD, on = ["predgene", "obsgene"], how = "left")
WHLBLD_trans_eQTL = pd.merge(WHLBLD_trans_eQTL, cis_eQTL, on = ["rsid", "predgene"], how = "left")
WHLBLD_trans_eQTL = WHLBLD_trans_eQTL[['predgene', 'predname', 'predChr', 'predS1', 'predS2', 'obsgene', 'obsname', 'obsChr', 'obsS1', 'obsS2', 'FHS_stat', 'FHS_beta', 'FHS_pval', 'FHS_FDR', 'DGN_stat', 'DGN_beta', 'DGN_pval', 'eQTLGen_trans_FDR', 'eQTLGen_cis_FDR']]
WHLBLD_trans_eQTL.loc[(WHLBLD_trans_eQTL.eQTLGen_trans_FDR >= 0), 'eQTLGen_trans_FDR'] = "Yes"
WHLBLD_trans_eQTL.loc[(WHLBLD_trans_eQTL.eQTLGen_trans_FDR.isnull()), 'eQTLGen_trans_FDR'] = "No"
WHLBLD_trans_eQTL.loc[(WHLBLD_trans_eQTL.eQTLGen_trans_FDR == "Yes") & (WHLBLD_trans_eQTL.eQTLGen_cis_FDR >= 0), 'eQTLGen_cis_FDR'] = "Yes"
WHLBLD_trans_eQTL.loc[(WHLBLD_trans_eQTL.eQTLGen_trans_FDR == "Yes") & (WHLBLD_trans_eQTL.eQTLGen_cis_FDR.isnull()), 'eQTLGen_cis_FDR'] = "No"
WHLBLD_trans_eQTL = WHLBLD_trans_eQTL.sort_values('eQTLGen_cis_FDR').drop_duplicates(subset = ['predname', 'obsname'], keep = "first") # https://stackoverflow.com/questions/32093829/pythonpandas-removing-duplicates-based-on-two-columns-keeping-row-with-max-va
WHLBLD_trans_eQTL = WHLBLD_trans_eQTL.rename(columns = {"eQTLGen_trans_FDR":"trans_eQTL_in_trans_eQTLGen", "eQTLGen_cis_FDR":"trans_eQTL_is_cis_eQTL"})
WHLBLD_trans_eQTL.to_csv("S1_WHLBLD_trans_eQTL.csv", index = False, na_rep = "NA")
print("Completed writing WHLBLD trans-PX and trans-eQTLGen overlapping results to WHLBLD_trans_eQTL.csv")

MultiXcan = pd.read_csv("TableS2_MultiXcan_results_2018-10-29.csv")
MultiXcan.columns = ['predgene', 'predname', 'predChr', 'predS1', 'predS2', 'obsgene', 'obsname', 'obsChr', 'obsS1', 'obsS2', 'FHS_stat', 'FHS_beta', 'FHS_pval', 'FHS_FDR', 'DGN_stat', 'DGN_beta', 'DGN_pval']
MultiXcan_trans_eQTL = pd.merge(MultiXcan, db_SNPs_trans_eQTLGen_noNA, on = ["predgene", "obsgene"], how = "left")
MultiXcan_trans_eQTL = pd.merge(MultiXcan_trans_eQTL, cis_eQTL, on = ["rsid", "predgene"], how = "left")
MultiXcan_trans_eQTL = MultiXcan_trans_eQTL[['predgene', 'predname', 'predChr', 'predS1', 'predS2', 'obsgene', 'obsname', 'obsChr', 'obsS1', 'obsS2', 'FHS_stat', 'FHS_beta', 'FHS_pval', 'FHS_FDR', 'DGN_stat', 'DGN_beta', 'DGN_pval', 'eQTLGen_trans_FDR', 'eQTLGen_cis_FDR']]
MultiXcan_trans_eQTL.loc[(MultiXcan_trans_eQTL.eQTLGen_trans_FDR >= 0), 'eQTLGen_trans_FDR'] = "Yes"
MultiXcan_trans_eQTL.loc[(MultiXcan_trans_eQTL.eQTLGen_trans_FDR.isnull()), 'eQTLGen_trans_FDR'] = "No"
MultiXcan_trans_eQTL.loc[(MultiXcan_trans_eQTL.eQTLGen_trans_FDR == "Yes") & (MultiXcan_trans_eQTL.eQTLGen_cis_FDR >= 0), 'eQTLGen_cis_FDR'] = "Yes"
MultiXcan_trans_eQTL.loc[(MultiXcan_trans_eQTL.eQTLGen_trans_FDR == "Yes") & (MultiXcan_trans_eQTL.eQTLGen_cis_FDR.isnull()), 'eQTLGen_cis_FDR'] = "No"
MultiXcan_trans_eQTL = MultiXcan_trans_eQTL.sort_values('eQTLGen_cis_FDR').drop_duplicates(subset = ['predname', 'obsname'], keep = "first") # https://stackoverflow.com/questions/32093829/pythonpandas-removing-duplicates-based-on-two-columns-keeping-row-with-max-va
MultiXcan_trans_eQTL = MultiXcan_trans_eQTL.rename(columns = {"eQTLGen_trans_FDR":"trans_eQTL_in_trans_eQTLGen", "eQTLGen_cis_FDR":"trans_eQTL_is_cis_eQTL"})
MultiXcan_trans_eQTL.to_csv("S3_MultiXcan_trans_eQTL.csv", index = False, na_rep = "NA")

print("Completed writing MultiXcan trans-PX and trans-eQTLGen overlapping results to MultiXcan_trans_eQTL.csv")
trans_eQTLGen_FDR_05.columns = ['eQTLGen_trans_P', 'rsid', 'predChr', 'rsid_pos', 'Zscore', 'A1', 'A0', 'obsgene', 'eQTLGen_gene_symbol', 'eQTLGen_gene_chr', 'eQTLGen_gene_pos', 'eQTLGen_nr_cohorts', 'eQTLGen_nr_samples', 'eQTLGen_trans_FDR']

WHLBLD_match_eQTLGen = pd.merge(WHLBLD, db_SNPs_trans_eQTLGen_WHLBLD, on = ["predgene", "obsgene"], how = "left")
WHLBLD_match_eQTLGen = pd.merge(WHLBLD_match_eQTLGen, cis_eQTL, on = ["rsid", "predgene"], how = "left")
WHLBLD_match_eQTLGen = WHLBLD_match_eQTLGen[['predgene', 'predname', 'predChr', 'predS1', 'predS2', 'obsgene', 'obsname', 'obsChr', 'obsS1', 'obsS2', 'FHS_stat', 'FHS_beta', 'FHS_pval', 'FHS_FDR', 'DGN_stat', 'DGN_beta', 'DGN_pval', 'rsid', 'eQTLGen_trans_P', 'eQTLGen_trans_FDR', 'eQTLGen_cis_P', 'eQTLGen_cis_FDR']].drop_duplicates()
WHLBLD_match_eQTLGen = WHLBLD_match_eQTLGen.dropna(subset = ['rsid'])
WHLBLD_match_eQTLGen.to_csv("S2_WHLBLD_match_eQTLGen.csv", index = False, na_rep = "NA")

MultiXcan_match_eQTLGen = pd.merge(MultiXcan, db_SNPs_trans_eQTLGen_noNA, on = ["predgene", "obsgene"], how = "left")
MultiXcan_match_eQTLGen = pd.merge(MultiXcan_match_eQTLGen, cis_eQTL, on = ["rsid", "predgene"], how = "left")
MultiXcan_match_eQTLGen = MultiXcan_match_eQTLGen[['predgene', 'predname', 'predChr', 'predS1', 'predS2', 'obsgene', 'obsname', 'obsChr', 'obsS1', 'obsS2', 'FHS_stat', 'FHS_beta', 'FHS_pval', 'FHS_FDR', 'DGN_stat', 'DGN_beta', 'DGN_pval', 'rsid', 'eQTLGen_trans_P', 'eQTLGen_trans_FDR', 'eQTLGen_cis_P', 'eQTLGen_cis_FDR', 'GTEx_tissue']].drop_duplicates()
MultiXcan_match_eQTLGen = MultiXcan_match_eQTLGen.dropna(subset = ['rsid'])#Old: 1 Mb
MultiXcan_match_eQTLGen.to_csv("S4_MultiXcan_match_eQTLGen.csv", index = False, na_rep = "NA")
print("Completed comparing trans-acting/target gene pairs b/w trans-PX and trans-eQTLGen.")

