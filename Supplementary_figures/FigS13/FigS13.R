#############################################
#############################################
#############################################
# Identify by descent (IBD) analysis using isoRelate R package- Only focusing monogenomic samples (n=76)
# See the Fig4 code for all samples 
# Write by Abebe Fola
# Date 06-06-2023
#############################################
#############################################
#############################################

#############################################
#############################################
# PART I- Created ped and map files using moimix R package  for isoRelate IBD analysis-moimix version 0.0.2.9001- R version 3.6.3
# If you already have ped and map file please ignore part I and start from Part II
#############################################
#############################################


# load libraries 
# install using devtools packages
# first install bioc dependencies
install.packages("BiocManager")
BiocManager::install("bahlolab/moimix", build_vignettes = TRUE)
BiocManager::install("gdsfmt")
BiocManager::install("SNPRelate")
BiocManager::install("SeqArray")
BiocManager::install("SeqVarTools")

library(moimix)
library(flexmix)
library(lattice)
library(tidyr)
library(ggplot2)
library(knitr)
library(kableExtra)
library(dplyr)
library(readr)
library(SeqVarTools)
library(gdsfmt)
library(SNPRelate)
library(SeqArray)
library(vcfR)
library(ape)

# Version used 

# R version 3.6.3 openintro_2.0.0   CMplot_3.6.2     
# MASS_7.3-51.5       qqman_0.1.4         tibble_3.1.0        ape_5.4-1           ggnetwork_0.5.8     ggplot2_3.3.3       igraph_1.2.6       
# dplyr_1.0.5         network_1.16.1      isoRelate_0.1.0     adegenet_2.1.3      ade4_1.7-16        

# Set working directory 

setwd("C:/Users/afola/Desktop/Zambia_MIS_addressing_NatMedcom_comments")

# load vcf file only monogenomic samples 
# bcftools view -S only_monogenomic.tsv Pf_capture241_final_MAFfilter_002.recode.vcf > Pf_capture241_final_MAFfilter_002onlymonogenomic.vcf

Zambia_monogenomicsamples<- ("Pf_capture238_final.SNPeffonlymonogenomic.vcf") # NB three samples where removed from downstream analysis
vcf_header <- seqVCF_Header(Zambia_monogenomicsamples)

# recode header format for AD from R to .
vcf_header$format$Number[vcf_header$ID == "AD"] <- "."

# info columns to retain
info.import <- c("AC", "AF", "AN", "BaseQRankSum", "DP", "DS",
                 "ExcessHet", "FS", "InbreedingCoeff", "MQ",
                 "MQRankSum", "QD", "ReadPosRankSum", "SOR", "ANN")

# format columns to retain
format.import <- c("AD", "DP", "GQ", "GT", "PL", "RGQ", "SB")

# convert VCF to GDS

Zambia_monogenomicsamples_gds<-seqVCF2GDS(Zambia_monogenomicsamples,
                              "Zambia_monogenomicsamples.gds",
                              header=vcf_header, info.import=info.import,
                              fmt.import=format.import)

#save sample identifiers

seqSummary(Zambia_monogenomicsamples_gds)

monoclonal_isolates <- seqOpen(Zambia_monogenomicsamples_gds) # Use isolates file for downstream analysis 

seqSummary(monoclonal_isolates)

# Get sample ids

sample.id <- seqGetData(monoclonal_isolates, "sample.id")

coords <- getCoordinates(monoclonal_isolates)

head(coords)

# Creating map and ped files for isorelate analysis 
#Extract PED files from gds. This function writes a plink .ped and .map file for a given gdsfile


zambia_ped_map_file_monoclonalonoclonal_isolates<-extractPED(monoclonal_isolates, moi.estimates = NULL, use.hets = F,
                                        outfile = NULL) # create map and ped file for isoRelate 
# Save Ped_map as.rds file
saveRDS(zambia_ped_map_file_monoclonalonoclonal_isolates, 
        file = "zambia_ped_map_file_monoclonal1_isolates.rds")
# Restore the object
zambia_ped_map_file_monoclonal_isolates<-readRDS(file = "zambia_ped_map_file_monoclonal1_isolates.rds")

#Uses big data format and load the following library() 
library(doMC)
registerDoMC(cores = 4)

# You can write ped and map file separately for further formating 

#write.table(zambia_ped_map_file_monoclonal_isolates$ped,file ="zambia_ped_file_monogenomic.txt",  sep = "\t",row.names = F, col.names = F ) # modify this file accordingly to added some variables and MIO at column 5

#write.table(zambia_ped_map_file_monoclonal_isolates$map,file ="zambia_map_file_monogenomic.txt", sep = "\t",row.names = F, col.names = F ) # SNP numbers should be equal to ped file and check the out put for centimorgan vs morgan.defoult is cM

#############################################
#############################################
#IBD analysis 
#############################################
#############################################

# install using devtools packages

#devtools::install_github("hadley/devtools"), devtools::install_github("bahlolab/isoRelate") , devtools::install_github("tidyverse/tibble")

library(isoRelate)
library(network)
library(dplyr)
library(igraph)
library(ggnetwork)
library(ggplot2)
library(ape)
library(stats)
library(tibble)

# Restore the object
zambia_ped_map_file_monoclonal_isolates<-readRDS(file = "zambia_ped_map_file_monoclonal_isolates.rds")

# reformat and filter genotypes
samplesmonogenomic_final_genotypes<-getGenotypes(ped.map = zambia_ped_map_file_monoclonal_isolates,
                                 reference.ped.map = NULL,
                                 maf = 0.02, # you can modify the threshold based on you data type
                                 isolate.max.missing = 0.35,# you can modify the threshold based on you data type
                                 snp.max.missing = 0.2,# you can modify the threshold based on you data type
                                 chromosomes = NULL,
                                input.map.distance = "cM",
                                reference.map.distance = "cM")#  centimorgan 


# Estimate Parameters


##Next we estimate the model parameters. These parameters can be useful as an initial measure of relatedness.
# Estimates the number of meiosis and the probabilities of sharing 0, 1 and 2 alleles IBD between all pairwise combinations of isolates

samplesmonogenomic_final_parameters <- getIBDparameters(ped.genotypes = samplesmonogenomic_final_genotypes,
                                               number.cores = 4) 

 head(samplesmonogenomic_final_parameters)

write.csv(samplesmonogenomic_final_parameters, "samplesmonogenomic_final_parameters.csv") # modify this file accordingly to reatian highly_related_pairs (pairs sharing>0.75)

#Assess probability of allele sharing

figpdf = paste('Probability of sharing 1 allele IBD.pdf', sep="")
pdf(file = figpdf)
hist(samplesmonogenomic_final_parameters$ibd1,breaks = 50, col="grey", xlab= " one allele sharing", main= "Probability of sharing 1 allele IBD between 
                   all pairwise combinations of isolates")
dev.off()

figpdf = paste('Probability of sharing zero allele IBD.pdf', sep="")
pdf(file = figpdf)
hist(samplesmonogenomic_final_parameters$ibd0,breaks = 50, col="grey", xlab= " two alleles sharing", main= "Probability of sharing 1 allele IBD between 
                   all pairwise combinations of isolates")
dev.off()


####################################

### Detect IBD Segments

####################################

#Following parameter estimation we can detect IBD segments
#detects genomic regions shared IBD between all pairwise combinations of isolates.
#We have set thresholds on the minimum number of SNPs and length of IBD segments reported in order to reduce false positive IBD calls that are most likely due to population linkage disequilibrium from extremely distant relatedness.

# infer IBD
monogenicsamples_final_ibd_segments<-getIBDsegments(ped.genotypes= samplesmonogenomic_final_genotypes,                    
                                             parameters = samplesmonogenomic_final_parameters, 
                                             number.cores = 4, #The number of cores used for parallel execution
                                             minimum.snps = 20, #the minimum number of SNPs in an IBD segment for it to be reported
                                             minimum.length.bp = 50000,# The minimum length of a reported IBD segment
                                             error = 0.001) #The genotyping error rate




####################################

### Exploring IBD Results

####################################

#isoRelate provides a number of graphical functions to help with exploring the vast number of IBD segments detected.
#We begin by plotting the IBD segments as colored blocks across the genome to get an idea of their size and distribution.

# plot IBD segments

pdf("Distribution of IBD segments across Zambia_final.pdf",width=15,height=20,paper='special')
plotIBDsegments(ped.genotypes = samplesmonogenomic_final_genotypes, 
                ibd.segments = monogenicsamples_final_ibd_segments, 
                interval = NULL,
                annotation.genes = NULL,
                annotation.genes.color = NULL,
                highlight.genes = NULL, 
                highlight.genes.labels = FALSE,
                highlight.genes.color = "red",
                highlight.genes.alpha = 0.1,
                segment.height = 0.5,
                number.per.page = NULL, # you can add page numbers to view it better 
                fid.label = F, 
                iid.label = F, 
                ylabel.size = 2, 
                add.rug =T,
                plot.title = "Distribution of IBD segments Only monoclonal samples", 
                add.legend = F, 
                segment.color = NULL) 

dev.off()

pdf("Distribution of IBD segments across and SNP density Zambia_final.pdf",width=15,height=20,paper='special')
plotIBDsegments(ped.genotypes = samplesmonogenomic_final_genotypes, 
                ibd.segments = monogenicsamples_final_ibd_segments, 
                interval = NULL,
                annotation.genes = NULL,
                annotation.genes.color = NULL,
                highlight.genes = NULL, 
                highlight.genes.labels = FALSE,
                highlight.genes.color = "red",
                highlight.genes.alpha = 0.1,
                segment.height = 0.5,
                number.per.page = NULL, # you can add page numbers to view it better 
                fid.label = F, 
                iid.label = F, 
                ylabel.size = 2, 
                add.rug =T,
                plot.title = "Distribution of IBD segments across Zambia", 
                add.legend = F, 
                segment.color = NULL) 

# generate a binary IBD matrix

monogenomicsamples_final_IBD_matrix<- getIBDmatrix(ped.genotypes = samplesmonogenomic_final_genotypes, 
                                             ibd.segments = monogenicsamples_final_ibd_segments)



# calculate the proportion of pairs IBD at each SNP for all samples together 

monogenomicsamples_final_proportion <- getIBDproportion(ped.genotypes = samplesmonogenomic_final_genotypes,
                                               ibd.matrix = monogenomicsamples_final_IBD_matrix, 
                                               groups = NULL)


head(monogenomicsamples_final_proportion)


# plot the proportion of pairs IBD

pdf("Proportion of pairs who are IBD at each SNP across across all isolates in Zambia_final.pdf", width=15,height=20,paper='special')
plotIBDproportions(ibd.proportions =monogenomicsamples_final_proportion, 
                   interval = NULL, 
                   annotation.genes = NULL,
                   annotation.genes.color = NULL,
                   highlight.genes = NULL,
                   highlight.genes.labels = TRUE,
                   highlight.genes.color = NULL,
                   highlight.genes.alpha = 0.1,
                   add.rug = F, 
                   plot.title = "Proportion of pairs who are IBD at each SNP across all isolates", 
                   add.legend = FALSE,
                   line.color = NULL, 
                   facet.label = F, 
                   facet.scales = "fixed", 
                   subpop.facet = FALSE)

dev.off()




###############
###############
# Fig13D - IBD distribution only monogenomic
###############
###############

IBDzambia_mono <-read.csv("samplesmonogenomic_IBD_parameters.csv")

mainplot2 <- IBDzambia_mono %>%
  ggplot() +
  geom_histogram(aes(x=ibd1, y = (..count../sum(..count..))*100),
                 color = "#000000", fill = "lightblue") +
  xlab("Pairwise relatedness (IBD)") + ylab("frequency (%)") +
  theme_classic()

insetplot <- IBDzambia_mono %>%
  ggplot() +
  geom_histogram(aes(x=ibd1, y = (..count../sum(..count..))*100),
                 color = "#000000", fill = "lightblue") +
  xlab("IBD>0.255") + ylab("frequency (%)") +
  theme_classic() +
  coord_cartesian(xlim = c(0.25,1), ylim = c(0,5)) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"))
# cowplot
cowplot::ggdraw() +
  cowplot::draw_plot(mainplot2, x = 0, y = 0, width = 1, height = 1, scale = 1) +
  cowplot::draw_plot(insetplot, x = 0.5, y= 0.3, width = 0.4, height = 0.4)

getwd()
# Save plot 
ggsave("FigonlymonoIBDdistribustion.svg", dpi=600, width=7.5, height=7)
ggsave("FigonlymonoIBDdistribustion.pdf", dpi=600, width=7.5, height=7)


###############
###############
# Fig1C - IBD distribution all samples 
###############
###############


IBDzambia<-read.csv("All_samples_samples_IBD_parameters.csv")

mainplot2 <- IBDzambia %>%
  ggplot() +
  geom_histogram(aes(x=ibd1, y = (..count../sum(..count..))*100),
                 color = "#000000", fill = "lightblue") +
  xlab("Pairwise relatedness (IBD)") + ylab("frequency (%)") +
  theme_classic()

insetplot <- IBDzambia %>%
  ggplot() +
  geom_histogram(aes(x=ibd1, y = (..count../sum(..count..))*100),
                 color = "#000000", fill = "lightblue") +
  xlab("IBD>0.25") + ylab("frequency (%)") +
  theme_classic() +
  coord_cartesian(xlim = c(0.25,1), ylim = c(0,5)) +
  theme_bw() +
  theme(panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"))
# cowplot
cowplot::ggdraw() +
  cowplot::draw_plot(mainplot2, x = 0, y = 0, width = 1, height = 1, scale = 1) +
  cowplot::draw_plot(insetplot, x = 0.5, y= 0.3, width = 0.4, height = 0.4)

getwd()
# Save plot 
ggsave("FigIBDdistribustion.svg", dpi=600, width=7.5, height=7)
ggsave("FigoIBDdistribustion.pdf", dpi=600, width=7.5, height=7)


