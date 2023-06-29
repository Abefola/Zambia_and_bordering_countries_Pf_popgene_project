
####################################################
# R script to calculate MAF, SNP and Sample missingness rate from VCF file 
#   Written by Abebe Fola    05/26/2022  
####################################################

##  Please install the R packages first
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("vcfR")
#BiocManager::install("SeqArray")
#BiocManager::install("adegenet")
#BiocManager::install("gdsfmt")
#BiocManager::install("SNPRelate")

# Load the R packages
library(gdsfmt)
library(SNPRelate)
library(vcfR)
library(adegenet)
library(ade4)
library(SeqArray)
library(SeqVarTools)
library(dplyr)
library(ggplot2)

#) Set your working directory

setwd("C:/Users/afola/Desktop/Zambia_2018_WGS_project_ms/Zambia_Popgene_ms/Suplement_figures/FigS2")


##!1) Please change the name of vcf and pop data files to the one you want to use

pop<-read.csv('population_metadata.csv', header=T,  sep=",") 

head(pop)#use province as pop

##!2) read in the VCF file with vcfR and check the file

vcf_238samples <- read.vcfR("Pf_capture238_final.SNPeff.vcf", verbose = FALSE)

##!3) vcf Summarization and check metadata

head(vcf_238samples) 
queryMETA(vcf_238samples) 
queryMETA(vcf_238samples, element = 'DP')

head(is.polymorphic(vcf_238samples, na.omit = TRUE))
head(is.biallelic(vcf_238samples))
queryMETA(vcf_238samples, element = 'FORMAT=<ID=DP')
strwrap(vcf_238samples@meta[1:7])

# If you needed you can Subset samples
vcf_238samples[,1:10]

# The fix region
# The fix region contains information for each variant which is sometimes summarized over all samples. The first eight columns of the fixed region and are titled CHROM, POS, ID, REF, ALT, QUAL, FILTER and INFO. This is per variant information which is ‘fixed’, or the same, over all samples. The first two columns indicate the location of the variant by chromosome and position within that chromosome. Here, the ID field has not been used, so it consists of missing data (NA). The REF and ALT columns indicate the reference and alternate allelic states. When multiple alternate allelic states are present they are delimited with commas. The QUAL column attempts to summarize the quality of each variant over all samples. The FILTER field is not used here but could contain information on whether a variant has passed some form of quality assessment.
head(getFIX(vcf_238samples))

# The gt region
#The gt (genotype) region contains information about each variant for each sample. The values for each variant and each sample are colon delimited. Multiple types of data for each genotype may be stored in this manner. The format of the data is specified by the FORMAT column (column nine). Here we see that we have information for GT, AD, DP, GQ and PL. The definition of these acronyms can be referenced by querying the the meta region, as demonstrated previously. Every variant does not necessarily have the same information (e.g., SNPs and indels may be handled differently), so the rows are best treated independently. Different variant callers may include different information in this region.
vcf_238samples@gt[1:6, 1:4]



#######################
# Calculate missingness 
#############

vcf238.fn <- "Pf_capture238_final.SNPeff.vcf"

# Reformat (change vcf to gds file )

snpgdsVCF2GDS(vcf238.fn, "samples238.gds", method="biallelic.only")

# Summary
snpgdsSummary("samples238.gds")

genofile_238samples <- snpgdsOpen("samples238.gds")

## get sample id and SNP id

sample.id <- read.gdsn(index.gdsn(genofile_238samples, "sample.id"))

snp.id <- read.gdsn(index.gdsn(genofile_238samples, "snp.id"))

# calculate sample missingness rate - checks samples containing not calls across the genome

samplemissr <- snpgdsSampMissRate(genofile_238samples, sample.id=sample.id, snp.id=snp.id, with.id=FALSE)

write.csv(samplemissr,"sample_missingness_238samples.csv")

figpdf = paste('FigS2A.pdf', sep="")
pdf(file = figpdf)
hist(samplemissr, breaks= 30, col="grey",   xlab="Missingness", main= "Sample Missingness")
abline(v=0.2, col="red", lwd=2, lty=2)

dev.off()

# calculate SNP missingness rate - checks individual missingness rate across samples
SNPmissrate<-snpgdsSNPRateFreq(genofile_238samples, sample.id=sample.id, snp.id=snp.id, with.id=FALSE,
                               with.sample.id=FALSE, with.snp.id=FALSE)

write.csv(SNPmissrate,"SNP_missingness_238samples.csv")

figpdf = paste('FigS2B.pdf', sep="")
pdf(file = figpdf)
hist(SNPmissrate$MissingRate, breaks= 30, col="grey",   xlab="Missingness", main= "SNP Missingness")
abline(v=0.2, col="red", lwd=2, lty=2)

dev.off()

#!6 Determine Minor allele frequency (MAF)

MAF<-SNPmissrate$MinorFreq

figpdf = paste('FigS2C.pdf', sep="")
pdf(file = figpdf)
hist(MAF, breaks= 30, col="grey",   ylab= "Frequency", xlab="MAF", main= "Minor Alelle frequency")
abline(v=0.02, col="red", lwd=2, lty=2)

dev.off()


#!NOTE: MAF filtering - we used 0.02 maf threshold for downstream data analysis
 
#Command line/ unix used to filter MAF using the following command 
#module load bioinfo vcftools/0.1.16

#vcftools --vcf <input vcf file> --recode --recode-INFO-all --maf 0.02 --out <output prefix>
  
sessionInfo()

