####################################################
# FigS11
# R script to do PCA for Zambia Pf population from VCF file using SNPRelate R package. 
# Written by Abebe Fola        
# Date 05/26/2022 
####################################################

##  Please install the R packages first
#if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("gdsfmt")
BiocManager::install("SNPRelate")

#Load the library

library(gdsfmt)
library(SNPRelate)
library(FactoClass)
library(scatterplot3d)
library(ComplexHeatmap)
library(RColorBrewer)
library(gplots)
library(vcfR)
library(dplyr)
library(plotly)
library(stats)
library(MASS)


##################################
# set working directory

setwd("C:/Users/afola/Desktop/Abefola_github/Zambia_2018_MIS_Pf_WGS_project/FigS11")

# Load your vcf file
# Please change the name of vcf file to the one you want to use
input='Pf_capture238_final_MAFfilter_002.recode.vcf'


# Please change the name of the prefix of output files if needed
name = "Zambia_maf0.02"


# Please change the file name if needed
# File format: tab delimited
# 1st column is sample name
# 2nd column is country name 
poplist<-read.table("meta_info.sample_id_province.tsv", sep="\t")
colnames(poplist) <- c("sample", "province")

# Please change the colors if needed, here are my color lists per province
col.list <- c("red","orange", "tomato4", "mediumvioletred","darkgreen", "darkblue", "dodgerblue")

vcf.fn <-input

Zambiamfa.02 = paste('Zambiamfa.02', name, '.gds', sep = "")
Zambiamfa.02 

## Read the VCF file and save it as GDS format file
snpgdsVCF2GDS(vcf.fn, Zambiamfa.02,  method="biallelic.only")

## Open the SNP GDS file
Zambiagenofile <- snpgdsOpen(Zambiamfa.02 )

#### PCA analysis
## To calculate the eigenvectors and eigenvalues for principal component analysis.
Zambia_pcamaf.02<-snpgdsPCA(Zambiagenofile, autosome.only=FALSE)

summary(Zambia_pcamaf.02)

## Get data "sample.id" from a GDS node
sample.id <- read.gdsn(index.gdsn(Zambiagenofile, "sample.id"))

EV1 = Zambia_pcamaf.02$eigenvect[,1]    # the first eigenvector
EV2 = Zambia_pcamaf.02$eigenvect[,2]    # the second eigenvector
EV3 = Zambia_pcamaf.02$eigenvect[,3]
## Get a list of countries following the order of the sample IDs
pop = factor(poplist$province)[match(Zambia_pcamaf.02$sample.id, poplist$sample)]


# Percentage of variance explained calculated as Eval(x)/sum(Eval)*100
#https://psico.fcep.urv.cat/utilitats/factor/documentation/Percentage_of_explained_common_variance.pdf
# https://support.minitab.com/en-us/minitab/18/help-and-how-to/modeling-statistics/multivariate/how-to/principal-components/interpret-the-results/key-results/

list(Zambia_pcamaf.02$eigenval[1:10])
# 3.375674 3.059228 2.779306 2.273331 2.240875 2.014508 1.910268 1.830409 1.783038 1.753333

y= as.numeric(Zambia_pcamaf.02$eigenval[1:10]/23*100)

y= as.numeric(format(round(y, 2), nsmall = 2))

figpdf = paste('varianceexp.pdf', sep="")
pdf(file = figpdf)
PVE=barplot(Zambia_pcamaf.02$eigenval[1:10]/23*100, ylim= c(0,20), type = "o", col = "grey", xlab = "PCA", ylab = "Percentage of variance explained",
            main = "Percentage of variance explained Zambia maf0.02")
#Add text at top of bars
text(x = PVE, y = y, label = y, pos = 3,cex = 0.8,  col = "black")

dev.off()


