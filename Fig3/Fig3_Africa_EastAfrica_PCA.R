
####################################################
####################################################
# R script to generate PCA plot from VCF file to parasite clustering 
# Written by Abebe Fola 
# Date 05/27/2022
####################################################
####################################################


####################################################
# FigA - African Map highlighting parasite populations included in this study
# We down loaded editable Africa map from the the following site (https://www.mapchart.net/africa.html) and colored to match PCA colors. 
# Written by Abebe Fola -05/27/2022
####################################################


##################################
##Fig3B- Pf3k-zambia-tanz data set (West and East Africa)
##################################


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


# set working directory 

setwd("C:/Users/afola/Desktop/Zambia_2018_WGS_project_ms/Zambia_Popgene_ms/Main_figures_final/Fig3")
# Please change the name of vcf file to the one you want to use
input='Pf_capture1001_final_MAFfilter_002.recode.vcf'


# Please change the name of the prefix of output files if needed
name = "pf3k_Zambia_maf0.02"


# Please change the file name if needed
# File format: tab delimited
# 1st column is sample name
# 2nd column is country name 
poplist <-read.table("meta_info.sample_id.country.txt", sep="\t")
colnames(poplist) <- c("sample", "country")

# Please change the colors if needed
col.list <- c("blue", "forestgreen", "darkslategray3", "red",
              "orchid1", "tan3","wheat4","purple","black") 

vcf <-input

pf3k_zambia = paste('pf3k_zambi_', name, '.gds', sep = "")
pf3k_zambia

## Read the VCF file and save it as GDS format file
snpgdsVCF2GDS(vcf, pf3k_zambia,  method="biallelic.only")

## Open the SNP GDS file
genofile <- snpgdsOpen(pf3k_zambia)

## To calculate the eigenvectors and eigenvalues for principal component analysis.
pf3k_zambia_pca<-snpgdsPCA(genofile, autosome.only=FALSE)

summary(pf3k_zambia_pca)

## Get data "sample.id" from a GDS node
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

EV1 = pf3k_zambia_pca$eigenvect[,1]    # the first eigenvector
EV2 = pf3k_zambia_pca$eigenvect[,2]    # the second eigenvector
EV3 = pf3k_zambia_pca$eigenvect[,3]
## Get a list of countries following the order of the sample IDs
pop = factor(poplist$country)[match(pf3k_zambia_pca$sample.id, poplist$sample)]

tail (pop)

# eigenvalues
EV=plot(pf3k_zambia_pca$eigenval[1:10], type = "o", col = "red", xlab = "PCA", ylab = "Eigenvalues",
        main = "Eigenvalues for 238 Zambia PCA analysis")

# Percentage of variance explained calculated as Eval(x)/sum(Eval)*100
#https://psico.fcep.urv.cat/utilitats/factor/documentation/Percentage_of_explained_common_variance.pdf
# https://support.minitab.com/en-us/minitab/18/help-and-how-to/modeling-statistics/multivariate/how-to/principal-components/interpret-the-results/key-results/


list(pf3k_zambia_pca$eigenval[1:10])
# 13.481939	8.833697	5.166541	3.984513	3.626854	3.424197	3.324931	3.275663	3.219526	3.152994


# variance proportion (%)
PV=plot(ccm_pca_eastA$varprop[1:10]*100, type = "o", col = "black", xlab = "PCA", ylab = "Variance proportion (%)",
        main = "Variance proportion for East Africa PCA analysis maf0.02")
        
PVE=barplot(pf3k_zambia_pca$eigenval[1:10]/51*100, type = "o", col = "lightblue", xlab = "PC", ylab = "Percentage of variance explained",
         main = "Percentage of variance explained for Pf3k-Zambia")

tab <- data.frame(pf3k_zambia_pca$sample.id,pop, EV1, EV2, stringsAsFactors = FALSE)
write.table(tab, file=paste(name,'-PCA.tsv', sep=""), 
            quote = F, row.names = F, sep="\t")
head(tab)
figpdf = paste('Fig4B.pdf', sep="")
pdf(file = figpdf)
plot(EV1, EV2,xlab="PC1 (26%)", ylab="PC2 (17)",xlim= c(-0.06, 0.06),ylim= c(-0.05, 0.05), main="Pf3k-Zambia PCA (outlier removed)", col=col.list[as.integer(pop)], pch=19,cex=0.3)
legend("bottomleft", legend=levels(pop), bg="transparent",pch=19, cex=0.7,col=col.list,text.col=col.list)
abline(v=0, h= 0, col="grey", cex= 0.3, lty = "dashed")
#text(EV1, EV2, pos = 1, labels = sample.id, cex =0.7 )
dev.off()
sessionInfo()


##########################
#Fig3C East Africa-Zambia data set
##########################

#) Please change the name of vcf file to the one you want to use
input='Pf_capture.EastAfrica.missness.MAFfilter_0.02.recode.vcf'

#) Please change the name of the prefix of output files if needed
name = "EastAfrica_maf0.02"


#) Please change the file name if needed
# File format: tab delimited
# 1st column is sample name
# 2nd column is country name 
poplist_east <-read.table("EA_meta_info.sample_id.country.txt", sep="\t")
colnames(poplist_east) <- c("sample", "country")

#) Please change the colors if needed
coleast <- c("blue", "red",
              "orchid1", "tan3","wheat4","purple","black") 


vcf.fn <-input

eastA= paste('eastA_', name, '.gds', sep = "")
eastA

## Read the VCF file and save it as GDS format file
snpgdsVCF2GDS(vcf.fn, eastA,  method="biallelic.only")

## Open the SNP GDS file
genofile <- snpgdsOpen(eastA)

#### PCA analysis
## To calculate the eigenvectors and eigenvalues for principal component analysis.
eastA_pca<-snpgdsPCA(genofile, autosome.only=FALSE)

summary(eastA_pca)

## Get data "sample.id" from a GDS node
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

EV1 = eastA_pca$eigenvect[,1]    # the first eigenvector
EV2 = eastA_pca$eigenvect[,2]    # the second eigenvector
EV3 = eastA_pca$eigenvect[,3]
popeA = factor(poplist_east$country)[match(eastA_pca$sample.id, poplist_east$sample)]


# eigenvalues
EV=plot(eastA_pca$eigenval[1:10], type = "o", col = "red", xlab = "PC", ylab = "Eigenvalues",
        main = "Eigenvalues for East Africa PCA analysis")

# Percentage of variance explained calculated as Eval(x)/sum(Eval)*100
#https://psico.fcep.urv.cat/utilitats/factor/documentation/Percentage_of_explained_common_variance.pdf
# https://support.minitab.com/en-us/minitab/18/help-and-how-to/modeling-statistics/multivariate/how-to/principal-components/interpret-the-results/key-results/


list(eastA_pca$eigenval[1:10])
#  5.401045 5.206530 3.531962 3.516118 3.341880 3.295215 3.190034 3.157343 3.143553 3.093130

PVE=barplot(eastA_pca$eigenval[1:10]/37*100, type = "o", col = "lightblue", xlab = "PC", ylab = "Percentage of variance explained",
            main = "Percentage of variance explained for EastAfrica")

tab <- data.frame(eastA_pca$sample.id,popeA, EV1, EV2, stringsAsFactors = FALSE)
write.table(tab, file=paste(name,'-PCA.tsv', sep=""), 
            quote = F, row.names = F, sep="\t")
head(tab)
figpdf = paste('Fig4C.pdf', sep="")
pdf(file = figpdf)
plot(EV1, EV2,xlab="PC1 (15%)", ylab="PC2 (13%)",  main= "Central/East Africa PCA (outlier removed) ",xlim=c(-0.08, 0.2),col=coleast[as.integer(popeA)], pch=19,cex=0.4)
legend("bottomleft", legend=levels(popeA), bg="transparent",pch=19, cex=0.7,col=col.list,text.col=col.list)
abline(v=0, h= 0, col="grey", cex= 0.3, lty = "dashed")
#text(EV1, EV2, pos = 1, labels = sample.id, cex =0.7 )
dev.off()
sessionInfo()


