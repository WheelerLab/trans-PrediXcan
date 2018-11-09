# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
# 
# Be sure to use an up to date version of R and Matrix eQTL.

library(MatrixEQTL)
"%&%" = function(a,b) paste(a,b,sep="")
args <- commandArgs(trailingOnly=T)
snpfiledir <- args[1] #snpfile directory
snpfile <- args[2]
snplocfile <- args[3]
expfiledir <- args[4] #obs exp file dir
obsfile <- args[5]
obslocfile <- args[6]
outprefix <- args[7] #output file prefix, e.g. SNP_GEU

## Location of the package with the data files.
#base.dir = find.package('MatrixEQTL');
# base.dir = '.';

## Settings

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Genotype file name
SNP_file_name = (snpfiledir %&% snpfile)
snps_location_file_name = (snpfiledir %&% snplocfile)

# Gene expression file name
expression_file_name = (expfiledir %&% obsfile)
gene_location_file_name = (expfiledir %&% obslocfile)

# Covariates file name
# Set to character() for no covariates
covariates_file_name = character() 

# Output file name
output_file_name_cis = tempfile();
output_file_name_tra = tempfile();

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 1;
pvOutputThreshold_tra = 1; #save all diff chr trans

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

# Distance for local gene-SNP pairs
cisDist = 3e8; #makes cis the same chromosome

## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = " ";      # the SPACE character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

me = Matrix_eQTL_main(
  snps = snps, 
  gene = gene, 
  #cvrt = cvrt,
  output_file_name     = output_file_name_tra,
  pvOutputThreshold     = pvOutputThreshold_tra,
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE, 
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos, 
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = 100,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

unlink(output_file_name_tra);
unlink(output_file_name_cis);

## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');

## Make the histogram plot of local and distant p-values
#png(filename = outprefix %&% '.meqtl.png', width = 650, height = 650)
#plot(me)
#dev.off()
write.table(me$cis$eqtls, outprefix %&% '.meqtl.cis.samechr.txt', quote = FALSE, row.names = FALSE, sep = "\t")
write.table(me$trans$eqtls, outprefix %&% '.meqtl.trans.diffchr.allres.txt', quote = FALSE, row.names = FALSE, sep = "\t")
