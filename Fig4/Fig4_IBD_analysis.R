#############################################
#############################################
#############################################
# Identify by descent (IBD) analysis using isoRelate R package- please check original pipeline https://github.com/bahlolab/isoRelate/blob/master/vignettes/introduction.Rmd
# Writen by Abebe Fola 06/03/2022
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

setwd("C:/Users/afola/Desktop/Zambia_projects/Zambia_2018_WGS_project_ms/Zambia_2018_MIS_Data_Analysis/IBD_analysis")

# load vcf file

Pf238_samples<- ("Pf_capture283_final.SNPeff.vcf") # NB three samples where removed from downstream analysis
vcf_header <- seqVCF_Header(Pf238_samples)

# recode header format for AD from R to .
vcf_header$format$Number[vcf_header$ID == "AD"] <- "."

# info columns to retain
info.import <- c("AC", "AF", "AN", "BaseQRankSum", "DP", "DS",
                 "ExcessHet", "FS", "InbreedingCoeff", "MQ",
                 "MQRankSum", "QD", "ReadPosRankSum", "SOR", "ANN")

# format columns to retain
format.import <- c("AD", "DP", "GQ", "GT", "PL", "RGQ", "SB")

# convert VCF to GDS

Pf238_samples_gds<-seqVCF2GDS(Pf238_samples,
                              "Pf238_samples.gds",
                              header=vcf_header, info.import=info.import,
                              fmt.import=format.import)

#save sample identifiers

seqSummary(Pf238_samples_gds)

isolates238samples_coregenome <- seqOpen(Pf238_samples_gds) # Use isolates file for downstream analysis 

seqSummary(isolates238samples_coregenome)

# Get sample ids

sample.id <- seqGetData(isolates238samples_coregenome, "sample.id")

coords <- getCoordinates(isolates238samples_coregenome)

head(coords)

# Creating map and ped files for isorelate analysis 
#Extract PED files from gds. This function writes a plink .ped and .map file for a given gdsfile


zambia_ped_map_file_238samples<-extractPED(isolates238samples_coregenome, moi.estimates = NULL, use.hets = T,
                                           outfile = NULL) # create map and ped file for isoRelate 
# Save Ped_map as.rds file
saveRDS(ped_map238samples_Nuegenome, file = "zambia_ped_map_file_238samples.rds")
# Restore the object
readRDS(file = "zambia_ped_map_file_238samples.rds")

#Uses big data format and load the following library() 
library(doMC)
registerDoMC(cores = 4)

# You can write ped and map file separately for further formating 

write.table(ped_map238samples_Nuegenome$ped,file ="zambia_ped_file_238samples.txt",  sep = "/t",row.names = F, col.names = F ) # modify this file accordingly to added some variables and MIO at column 5

write.table(ped_map238samples_Nuegenome$map,file ="zambia_map_file_238samples.txt", sep = "/t",row.names = F, col.names = F ) # SNP numbers should be equal to ped file and check the out put for centimorgan vs morgan.defoult is cM

#############################################
#############################################
#PART II - IBD analysis 
#############################################
#############################################

# install using devtools packages

devtools::install_github("hadley/devtools")

# load libraries 

devtools::install_github("bahlolab/isoRelate") 
devtools::install_github("tidyverse/tibble")

library(isoRelate)
library(network)
library(dplyr)
library(igraph)
library(ggnetwork)
library(ggplot2)
library(ape)
library(stats)
library(tibble)

###########################
# Check isoRelate package using software inbuilt test files - it takes few minutes to complete test file analysis
###########################

# Test files

str(png_pedmap)

my_genotypes <- getGenotypes(ped.map = png_pedmap,
                             reference.ped.map = NULL,
                             maf = 0.01,
                             isolate.max.missing = 0.1,
                             snp.max.missing = 0.1,
                             chromosomes = NULL,
                             input.map.distance = "cM",
                             reference.map.distance = "cM")



my_parameters <- getIBDparameters(ped.genotypes = my_genotypes, 
                                  number.cores = 4)
head(my_parameters)

my_ibd <- getIBDsegments(ped.genotypes = my_genotypes,
                         parameters = my_parameters, 
                         number.cores = 4, 
                         minimum.snps = 20, 
                         minimum.length.bp = 50000,
                         error = 0.001)
head(my_ibd)

getIBDsummary(ped.genotypes = my_genotypes, 
              ibd.segments = my_ibd) 


str(png_pedmap)


my_matrix <- getIBDmatrix(ped.genotypes = my_genotypes, 
                          ibd.segments = my_ibd)
# calculate the proportion of pairs IBD at each SNP
my_proportion <- getIBDproportion(ped.genotypes = my_genotypes, 
                                  ibd.matrix = my_matrix, 
                                  groups = NULL)
# plot the proportion of pairs IBD
plotIBDproportions(ibd.proportions = my_proportion, 
                   interval = NULL, 
                   annotation.genes = NULL,
                   annotation.genes.color = NULL,
                   highlight.genes = NULL,
                   highlight.genes.labels = TRUE,
                   highlight.genes.color = NULL,
                   highlight.genes.alpha = 0.1,
                   add.rug = FALSE, 
                   plot.title = "Proportion of pairs IBD in PNG", 
                   add.legend = FALSE,
                   line.color = NULL, 
                   facet.label = TRUE, 
                   facet.scales = "fixed", 
                   subpop.facet = FALSE)

summary(my_proportion$prop_ibd)


# generate the isolates who are IBD over the Plasmodium falciparum CRT gene
my_i_clusters <- getIBDiclusters(ped.genotypes = my_genotypes, 
                                 ibd.segments = my_ibd, 
                                 interval = c("Pf3D7_07_v3", 403222, 406317), 
                                 prop=0, 
                                 hi.clust = FALSE)
str(my_i_clusters)

my_groups <- my_genotypes[[1]][,1:3]
my_groups[1:10,"pid"] <- "a"
my_groups[11:25,"pid"] <- "b"
my_groups[26:38,"pid"] <- "c"

# plot the network of clusters
plotIBDclusters(ped.genotypes = my_genotypes, 
                clusters = my_i_clusters, 
                groups = my_groups, 
                vertex.color = NULL, 
                vertex.frame.color = "white",
                vertex.size = 4, 
                vertex.name = FALSE, 
                edge.color = "gray60", 
                edge.width = 0.8, 
                mark.border = "white",
                mark.col = "gray94", 
                add.legend = TRUE, 
                legend.x = -1.5, 
                legend.y = -0.25, 
                layout = NULL, 
                return.layout = FALSE)


####################################
# Once the test file checked now you can move to your actual data analysis
####################################

####################################
###get data and combine ped and map file to create pedmap file
####################################

Zambia_pedmap_samples238_final <- list()

Zambia_pedmap_samples238_final[[1]]= read.delim("zambia_ped_file_238samples.txt",  header=FALSE,stringsAsFactors=FALSE) # For large data it take time

Zambia_pedmap_samples238_final[[2]] <- read.delim("zambia_map_file_238samples.txt", header=FALSE,
                                                  stringsAsFactors=FALSE) %>% mutate(V2=as.character(V2)) 

#str(Zambia_pedmap_samples238_final)

#Get metadata and create province_colours_sorted.txt or put ped and map file in provincial sequential order and use pid id (1st column) from ped file to color isolates per province (I am using the second option)

prov_colour <- read.delim("zambia_province_colours_sorted.txt", col.names=c( "Province", "Colour"), header=FALSE, stringsAsFactors=FALSE) 

# reformat and filter genotypes
samples238_final_genotypes<-getGenotypes(ped.map = Zambia_pedmap_samples238_final,
                                         reference.ped.map = NULL,
                                         maf = 0.02, # you can modify the threshold based on you data type
                                         isolate.max.missing = 0.35,# you can modify the threshold based on you data type
                                         snp.max.missing = 0.2,# you can modify the threshold based on you data type
                                         chromosomes = NULL,
                                         input.map.distance = "cM",
                                         reference.map.distance = "cM")#  centimorgan 


####################################
### Estimate Parameters
####################################

##Next we estimate the model parameters. These parameters can be useful as an initial measure of relatedness.
# Estimates the number of meiosis and the probabilities of sharing 0, 1 and 2 alleles IBD between all pairwise combinations of isolates

samples238_final_parameters <- getIBDparameters(ped.genotypes = samples238_final_genotypes,
                                                number.cores = 4) 

head(samples238_final_parameters)

write.csv(samples238_final_parameters, "samples238_final_parameters.csv") # modify this file accordingly to reatian highly_related_pairs (pairs sharing>0.75)

#Assess probability of allele sharing

figpdf = paste('Probability of sharing 1 allele IBD.pdf', sep="")
pdf(file = figpdf)
hist(samples238_final_parameters$ibd1,breaks = 50, col="grey", xlab= " one allele sharing", main= "Probability of sharing 1 allele IBD between 
                   all pairwise combinations of isolates")
dev.off()

figpdf = paste('Probability of sharing zero allele IBD.pdf', sep="")
pdf(file = figpdf)
hist(samples238_final_parameters$ibd0,breaks = 50, col="grey", xlab= " two alleles sharing", main= "Probability of sharing 1 allele IBD between 
                   all pairwise combinations of isolates")
dev.off()


####################################

### Detect IBD Segments

####################################

#Following parameter estimation we can detect IBD segments
#detects genomic regions shared IBD between all pairwise combinations of isolates.
#We have set thresholds on the minimum number of SNPs and length of IBD segments reported in order to reduce false positive IBD calls that are most likely due to population linkage disequilibrium from extremely distant relatedness.

# infer IBD
samples238_final_ibd_segments<-getIBDsegments(ped.genotypes= samples238_final_genotypes,                    
                                              parameters = samples238_final_parameters, 
                                              number.cores = 4, #The number of cores used for parallel execution
                                              minimum.snps = 20, #the minimum number of SNPs in an IBD segment for it to be reported
                                              minimum.length.bp = 50000,# The minimum length of a reported IBD segment
                                              error = 0.001) #The genotyping error rate



#bp ibd summary

head(samples238_final_ibd_segments)

summary(samples238_final_ibd_segments$length_bp)

write.table(samples238_final_ibd_segments, file="samples238_final_ibd_segments_final.txt", quote=FALSE, row.names = F, sep="/t")

ibdbp<- density(samples238_final_ibd_segments$length_bp)

plot(ibdbp, xlab="Length of segement (bp)", main = "Distribution of length of IBD segments across genome")

polygon(ibdbp,  col="dodgerblue", border="black",lty = 1, lwd = 1,angle = c(-45, 45), bg= "red", ylab="Frequency",main="Paiwise allele sharing")

#SNPs 
summary(samples238_final_ibd_segments$number_snps)

bdsnp<- density(samples238_final_ibd_segments$number_snps)

plot(bdsnp, xlab="Number of SNPs", main = "Distribution of SNPs across IBD segments")

polygon(bdsnp,  col="orange", border="black",lty = 1, lwd = 1,angle = c(-45, 45), bg= "red", ylab="Frequency",main="Paiwise allele sharing")

# Centimorgan

summary(samples238_final_ibd_segments$length_M*100)

Ibdcm<- density(samples238_final_ibd_segments$length_M*100)

plot(Ibdcm, xlab="Length of segement (cM)", main = "Distribution of length of IBD segments across genome")

polygon(Ibdcm,  col="darkred", border="black",lty = 1, lwd = 1,angle = c(-45, 45), bg= "red", ylab="Frequency",main="Paiwise allele sharing")


#We can get a brief summary of the inferIBD segments as follows

getIBDsummary(ped.genotypes = samples238_final_genotypes, 
              ibd.segments = samples238_final_ibd_segments) 


#IBD segments are found using the Viterbi algorithm, which finds the single most likely sequence of IBD states that could have generated the observed genotypic data. An alternative method to this is to calculate the posterior probability of IBD sharing, which calculates the probability of sharing 0, 1 or 2 alleles IBD at each SNP, given the genotypic data. Thus, in addition to the Viterbi algorithm, we provide a function to generate the average posterior probability of IBD sharing for each isolate pair, which is calculated as; 

#avePostPr = /frac{PostPr(IBD = 1)}{2} + PostPr(IBD = 2).$$
# **Please note, the output from the following function has number of rows equal to the number of SNPs and number of columns equal to the number of pairwise comparisons. As such, it will be very large when there are many isolates in the analysis and is likely to take some time to run. Additionally, the output may exceed R's memory limit and, in some instances, can cause R to unexpectedly crash.**


# calculate the posterior probability of IBD sharing
samples238_final_posterior_probability <- getIBDposterior(ped.genotypes=samples238_final_genotypes,
                                                          parameters = samples238_final_parameters, 
                                                          number.cores =4, 
                                                          error = 0.001)


head(samples238_final_posterior_probability[,1:10])

#write.table(samples238_final_posterior_probability, file="samples238_final_posterior_probability.txt", quote=FALSE, row.names = F, sep="/t")

####################################

### Exploring IBD Results

####################################

#isoRelate provides a number of graphical functions to help with exploring the vast number of IBD segments detected.
#We begin by plotting the IBD segments as colored blocks across the genome to get an idea of their size and distribution.

# plot IBD segments

pdf("Distribution of IBD segments across Zambia_final.pdf",width=15,height=20,paper='special')
plotIBDsegments(ped.genotypes = samples238_final_genotypes, 
                ibd.segments = samples238_final_ibd_segments, 
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
                add.rug =F,
                plot.title = "Distribution of IBD segments across Zambia", 
                add.legend = F, 
                segment.color = NULL) 

dev.off()

pdf("Distribution of IBD segments across and SNP density Zambia_final.pdf",width=15,height=20,paper='special')
plotIBDsegments(ped.genotypes = samples238_final_genotypes, 
                ibd.segments = samples238_final_ibd_segments, 
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

dev.off()

#Isolate pairs who are highly related will have many IBD segments spanning a large amount of the genome. 
#Identical isolates will have IBD blocks spanning the entire genome.

         


# generate a binary IBD matrix

samples238_final_IBD_matrix<- getIBDmatrix(ped.genotypes = samples238_final_genotypes, 
                                           ibd.segments = samples238_final_ibd_segments)



# calculate the proportion of pairs IBD at each SNP for all samples together 

samples238_final_proportion <- getIBDproportion(ped.genotypes = samples238_final_genotypes,
                                                ibd.matrix = samples238_final_IBD_matrix, 
                                                groups = NULL)


head(samples238_final_proportion)




# plot the proportion of pairs IBD

summary(samples238_final_proportion$prop_ibd)

write.csv( samples238_final_proportion, "samples238_final_proportion_across_provinces.csv")


dd<-density(samples238_final_proportion$prop_ibd*100)
plot(dd,lty = 1, lwd = 1,angle = c(-45, 45), bg= "red", xlab= " Proportion of pairs", main="The proportion of pairs who are IBD at each SNP")

polygon(dd,  col="blue", border="black",lty = 1, lwd = 1,angle = c(-45, 45), bg= "red", ylab="Frequency",main="Paiwise allele sharing")


# plot the proportion of pairs IBD

pdf("Proportion of pairs who are IBD at each SNP across across all isolates in Zambia_final.pdf", width=15,height=20,paper='special')
plotIBDproportions(ibd.proportions =samples238_final_proportion, 
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

# creating a stratification data set 
samples238_groups <- samples238_final_genotypes[[1]][,1:2]
samples238_groups[1:14,"pid"] <- "Copperbelt"
samples238_groups[15:30,"pid"] <- "Eastern"
samples238_groups[31:156,"pid"] <- "luapula"
samples238_groups[157:173,"pid"] <- "Muchinga"
samples238_groups[174:189,"pid"] <- "North-Western"
samples238_groups[190:210,"pid"] <- "Northern"
samples238_groups[211:238,"pid"] <- "Western"


head(samples238_groups)

tail(samples238_groups)

# My color codes

colors <-c("red","orange", "tomato4","mediumvioletred","darkgreen", "darkblue", "dodgerblue")

# calculate the proportion of pairs IBD at each SNP
samples238_groups_proportion_final <- getIBDproportion(ped.genotypes =samples238_final_genotypes,
                                                       ibd.matrix =  samples238_final_IBD_matrix, 
                                                       groups = samples238_groups)


# plot the proportion of pairs IBD for group within and between

pdf("Proportion of pairs IBD within and between parasite population in Zambia.pdf", width=20,height=15,paper='special')
plotIBDproportions(ibd.proportions=samples238_groups_proportion_final,
                   interval = NULL, 
                   annotation.genes = NULL,
                   annotation.genes.color = NULL,
                   highlight.genes = NULL,
                   highlight.genes.labels = TRUE,
                   highlight.genes.color = NULL,
                   highlight.genes.alpha = 0.1,
                   add.rug = FALSE, 
                   plot.title = "Proportion of pairs IBD within and between population in Zambia", 
                   add.legend = FALSE,
                   line.color = NULL, 
                   facet.label = TRUE, 
                   facet.scales = "fixed", 
                   subpop.facet = T)

dev.off()

# plot the proportion of pairs IBD within province only

write.csv( samples238_groups_proportion_final, "samples238_groups_proportion_final.csv") # Modify this data and only save within province output with the following format -#must be a data.frame with 7 columns: chr, snp_id, pos_M, pos_bp, pop, subpop and prop_ibd

samples238_groups_proportion_mod <-read.csv('samples238_groups_proportion_final_per_province_only.csv', header=T,  sep=",") 



pdf("Proportion of pairs who are IBD at each SNP within province final.pdf", width=20,height=15,paper='special')
plotIBDproportions(ibd.proportions =  samples238_groups_proportion_mod,
                   interval = NULL,
                   annotation.genes = NULL,
                   annotation.genes.color = NULL,
                   highlight.genes = NULL,
                   highlight.genes.labels = FALSE,
                   highlight.genes.color = NULL,
                   highlight.genes.alpha = 0.1,
                   line.color = colors,
                   add.rug = F,
                   plot.title = "Proportion of pairs who are IBD at each SNP within province",
                   add.legend = T,
                   facet.label = TRUE,
                   facet.scales = "fixed",
                   subpop.facet = T)
dev.off()

###############################################

### calculate the significance of IBD sharing       

###############################################

# To identify genomic loci with significant amounts of excess IBD we can apply a transformation to the binary IBD matrix to 
# account for variations in isolate relatedness as well as SNP allele frequencies, then calculate a summary statistic 
#  that can be used to assess significance.

# calculate the significance of IBD sharing
samples238_final_iR <- getIBDiR(ped.genotypes =samples238_final_genotypes,
                                ibd.matrix =  samples238_final_IBD_matrix, 
                                groups = NULL)

# plot the iR statistics
pdf("Significant IBD sharing across parasite population in Zambiamaf02.pdf", width=20,height=15,paper='special')
plotIBDiR(ibd.iR = samples238_final_iR, 
          interval = NULL, 
          annotation.genes = NULL,
          annotation.genes.color = NULL,
          highlight.genes = NULL,
          highlight.genes.labels =T,
          highlight.genes.color = NULL,
          highlight.genes.alpha = 0.4,
          point.size = 0.5,
          point.color = NULL,
          add.rug = F, 
          plot.title = "Significant IBD sharing across parasite population in Zambia", 
          add.legend = F,
          facet.label = TRUE, 
          facet.scales = "fixed")                 

dev.off()       


# calculate the significance of IBD sharing for groups 

samples238_final_iRgroup <- getIBDiR(ped.genotypes =samples238_final_genotypes,
                                     ibd.matrix =  samples238_final_IBD_matrix, 
                                     groups =   samples238_groups )

# plot the iR statistics
pdf("Significant IBD sharing within and between populations in Zambia.pdf", width=20,height=15,paper='special')
plotIBDiR(ibd.iR = samples238_final_iRgroup , 
          interval = NULL, 
          annotation.genes = NULL,
          annotation.genes.color = NULL,
          highlight.genes = NULL,
          highlight.genes.labels =T,
          highlight.genes.color = NULL,
          highlight.genes.alpha = 0.4,
          point.size = 1,
          point.color = NULL,
          add.rug = T, 
          plot.title = "Significant IBD sharing within and between populations in Zambia", 
          add.legend = T,
          facet.label = TRUE, 
          facet.scales = "fixed")              

dev.off()          


## write the file in .csv and format accordingly 
write.csv(samples238_final_iR, "samples238_final_iR_across_province.csv") # modify accordigly and plot manhton plot


# Make the Manhattan plot on the Results dataset
#https://www.r-graph-gallery.com/101_Manhattan_plot.html
library(qqman)
library(calibrate)
library("CMplot")
#install.packages("CMplot")

samples238_iR1<- read.csv("samples238_iR_manhattan_formated.csv", header=TRUE)  # formatted as manhattan plot format 
#(column= chromosome number 1, 2, 3, etc,column 2= bp number, column3 = snp_id like and Pf3D7_11_v3:592805 and last column= P.value. NB you have to change log10(p.value) result from isoRelate to p.value =1/10^log10(p.value))

manhattan(samples238_iR1, col = c("black", "royalblue"), annotateTop= T, suggestiveline = F, genomewideline = F, ylim= c(0,20), chr="CHR", bp="BP", snp="SNP", p="P" , main= "Manhattan plot of the significant IBD sharing")

##qqplot It is a good practice to draw a qqplot from the output of a GWAS. It allows to compare the distribution of the pvalue with an expected distribution by chance. Its realization is straightforward thanks to the qq function:

qq(samples238_iR1$P, main="Distribution of the P.value with an expected distribution by chance")

##annotate You probably want to know the name of the SNP of interest: the ones with a high pvalue. You can automatically annotate them using the annotatePval argument:

manhattan(samples238_iR1, ylim= c(0,20), annotatePval = 0.01, main= "Manhattan plot of the significant IBD sharing")


#circular version with CMplot https://github.com/YinLiLin/CMplot

samples238_iR_circular<- read.csv("samples238_iR_COMPLOT_formated.csv", header = T) #format the data according to cmplot format

CMplot(samples238_iR_circular, plot.type="c")

#Circular-Manhattan plot


CMplot(samples238_iR_circular,type="p",plot.type="c",r=0.4,cir.legend=TRUE,
       outward=FALSE,cir.legend.col="black",cir.chr.h=1.3,chr.den.col="black",file="jpg",
       memo="",dpi=300,file.output=TRUE,verbose=TRUE,width=10,height=10)


CMplot(samples238_iR_circular,type="p",plot.type="d",bin.size=1e6,chr.den.col=c("darkgreen", "yellow", "red"),file="jpg",memo="",dpi=300,
       file.output=TRUE,verbose=TRUE,width=9,height=6)


CMplot(samples238_iR_circular,type="p",plot.type="c",r=0.4,col=c("grey30","grey60"),
       threshold=c(1e-6,1e-4),cir.chr.h=1.5,amplify=TRUE,threshold.lty=c(1,2),threshold.col=c("red",
                                                                                                   "blue"),signal.line=1,signal.col=c("red","green"),chr.den.col=c("darkgreen","yellow","red"),
                                                                                                   bin.size=1e6,outward=FALSE,file="jpg",memo="",dpi=300,file.output=TRUE,verbose=TRUE,width=10,height=10)

CMplot(samples238_iR_circular,type="p",plot.type="c",LOG10=FALSE,outward=TRUE,col=matrix(c("#4DAF4A",NA,NA,"dodgerblue4",
                                                                                                    "deepskyblue",NA,"dodgerblue1", "olivedrab3", "darkgoldenrod1"), nrow=3, byrow=TRUE),
                                                                                                    threshold=NULL,r=1.2,cir.chr.h=1.5,cir.legend.cex=0.5,
       cir.band=1,file="jpg", memo="",dpi=300,chr.den.col="black",file.output=TRUE,verbose=TRUE,
       width=10,height=10)

# PLOT

CMplot(samples238_iR_circular,type="p",plot.type="m",LOG10=TRUE,threshold=NULL,file="jpg",memo="",dpi=300,
       file.output=TRUE,verbose=TRUE,width=14,height=6,chr.labels.angle=45)

#SNP-density plot

CMplot(samples238_iR_circular,plot.type="d",bin.size=25e3,chr.den.col=c("darkgreen", "yellow", "orange", "red"),file="jpg",memo="",dpi=300,
       file.output=TRUE,verbose=TRUE,width=9,height=6)          

# plot the significant IBD sharing only within populations in Zambia

colors <-c("red","orange", "tomato4","mediumvioletred","darkgreen", "darkblue", "dodgerblue")

write.csv(samples238_final_iRgroup, "samples238_final_iRgroup.csv") # write as csv


samples238_iRgroup_mod <-read.csv('samples238_iRgroup_final_mod.csv', header=T,  sep=",") #must be a data.frame with 7 columns: chr, snp_id, pos_M, pos_bp, pop, subpop and prop_ibd


# plot the iR statistics only within population

pdf("Significant IBD sharing within parasite population in Zambia-final.pdf", width=20,height=15,paper='special')
plotIBDiR(ibd.iR = samples238_iRgroup_mod , 
          interval = NULL, 
          annotation.genes = NULL,
          annotation.genes.color = NULL,
          highlight.genes = NULL,
          highlight.genes.labels =T,
          highlight.genes.color = NULL,
          highlight.genes.alpha = 0.4,
          point.size = 0.2,
          point.color = colors,
          add.rug = T, 
          plot.title = "Significant IBD sharing within parasite population in Zambia", 
          add.legend = T,
          facet.label = TRUE, 
          facet.scales = "fixed")  


dev.of



# =========================================
# =========================================
#################################   
#plot the network of clusters - use posterior_probability value to plot IBD network
##################################################

#We can investigate the specific isolate pairs that are contributing to the excess IBD sharing on chromosome 13 by creating an IBD network.
# Such a network will identify clusters of isolates sharing a common haplotype over this region. We include our example grouping of isolates
# when plotting the network.


IBD_metrics <- list()
IBD_metrics[["fraction_IBD"]] <- (samples238_final_posterior_probability[,5:ncol(samples238_final_posterior_probability)]*2) %>% 
  colMeans() %>% as.data.frame() %>% 
  tibble::rownames_to_column() %>% 
  tidyr::extract(rowname, c("iid1", "iid2"),
                 regex="[A-Za-z0-9_()//.-]+/([A-Za-z0-9_()//.-]+)/[A-Za-z0-9_()-]+/([A-Za-z0-9_()-]+)") %>%
  mutate(iida=pmin(iid1, iid2), iidb=pmax(iid1, iid2)) %>%
  dplyr::select(-iid1, -iid2)
colnames(IBD_metrics[["fraction_IBD"]]) <- c("fraction_IBD", "iida", "iidb")
IBD_metrics[["fraction_IBD"]]$fraction_IBD <- IBD_metrics[["fraction_IBD"]]$fraction_IBD*100

write.csv(IBD_metrics, "IBD_metrics_final.csv")

# ======================== metadata Network analysis-grouped provinces into bigger region to minimize sample size variation =============================

metadata <- samples238_final_genotypes[["pedigree"]][, c("iid", "fid")] %>% 
  tidyr::extract(fid, c("Region", "Pop"), convert=TRUE, remove=FALSE,
                 regex="([A-Za-z-_]+)([0-9]+)") %>%
  mutate(Region=gsub("-|_", "", Region) %>% stringr::str_to_title(),
         Region=ifelse(Region=="Regionone", "Luapula-Northern",Region),
         Region=ifelse(Region=="Regiontwo", "Eastern-Muchinga",Region),
         Region=ifelse(Region=="Regionthree", "Copperbelt-NorthWestern",Region),
         Region=ifelse(Region=="Regionfour", "Western",Region),
         Pop=ifelse(Pop=="101", "Copperbelt",Pop),
         Pop=ifelse(Pop=="102", "Eastern",Pop),
         Pop=ifelse(Pop=="103", "Luapula",Pop),
         Pop=ifelse(Pop=="104", "Muchinga",Pop),
         Pop=ifelse(Pop=="105", "North-Western",Pop),
         Pop=ifelse(Pop=="106", "Northern",Pop),
         Pop=ifelse(Pop=="107", "Western",Pop),
         Province=ifelse(Region=="Luapula-Northern","Regionone", "Others")) %>%
  dplyr::rename(Sample=iid, Study=fid) 

rownames(metadata) <- metadata$Sample

population_levels <- c("Study", "Region", "Pop", "Province")

pop_colours <- pop_shapes <- list()

pop_colours[["Pop"]] <- data.frame(row.names=c("Copperbelt", "Eastern", "Luapula", "Muchinga", "North-Western", "Northern","Western"),
                                   Color=c("red","orange", "tomato4","mediumvioletred","darkgreen", "darkblue", "dodgerblue"))

pop_colours[["Study"]] <- read.delim("zambia_province_colours_yearformated_sorted.txt", col.names=c("Pop", "Color"), header=FALSE, stringsAsFactors=FALSE) %>% tibble::column_to_rownames(var="Pop")


pop_shapes[["Region"]] <- data.frame(row.names=c("Luapula-Northern", "Eastern-Muchinga", "Copperbelt-NorthWestern", "Western"),
                                     Shape=c(16, 16, 16, 16))


# ======================== generate IBD networks =================

# generate layout and igraph object for pairs of isolates with IBD sharing over
# an arbitrary threshold
set_pairwise_threshold <- function(threshold, IBD_metrics, metric, metadata, population_levels, province=NULL) {
  pair_net <- list()
  if (is.null(province)) {
    ibd_sum <- IBD_metrics[[metric]]
  } else {
    keep_samples <- (metadata %>% subset(Province==province))$Sample
    ibd_sum <- IBD_metrics[[metric]] %>% 
      subset(iida %in% keep_samples & iidb %in% keep_samples)
  }
  ibd_sum$length <- ibd_sum[, metric]
  pairs_ibd <- ibd_sum %>% subset(length>=threshold) %>%
    dplyr::select(iida, iidb, length)
  
  # represent all isolates on network, irrespective of IBD sharing
  nodes <- metadata[as.character(unique(c(ibd_sum$iida, ibd_sum$iidb))), 
                    c("Sample", population_levels)] %>%
    subset(!is.na(Sample))
  
  # generate graph with weighted edges
  pair_net[["graph"]] <- igraph::graph.data.frame(pairs_ibd, vertices=nodes, directed=F)
  
  #zambia.network <- ggnetwork::ggnetwork(pair_net[["graph"]])
  #zambia.network.i.network <- zambia.network[!duplicated(zambia.network[,c(1,2,4)]),]
  #pair_net[["layout"]] <- as.matrix(zambia.network.i.network[,1:2], ncol=2)
  
  # isolates with no IBD sharing relegated to the periphery
  pair_net[["layout"]] <- layout_with_fr(pair_net[["graph"]], weights = NULL)
  colnames(pair_net[["layout"]]) <- c("x", "y")
  
  E(pair_net[["graph"]])$weight <- pairs_ibd$length
  return(pair_net)
}

# given a predefined graph/layout for pairs of isolates with IBD sharing above
# a particular threshold, generate IBD network by a particular population level
plot_pairwise_network <- function(pair_net, pop_colours, pop_shapes, title="") {
  my_nodes <- pair_net[["layout"]] %>% as.data.frame %>% 
    mutate(Sample=as_ids(V(pair_net[["graph"]])),
           Pop=as.character(vertex.attributes(pair_net[["graph"]])[["Pop"]]),
           Region=vertex.attributes(pair_net[["graph"]])[["Region"]])
  
  my_edges <- igraph::as_data_frame(pair_net[["graph"]]) %>%
    merge(my_nodes, by.x="from", by.y="Sample") %>%
    merge(my_nodes, by.x="to", by.y="Sample") %>%
    transmute(x=x.x, xend=x.y, y=y.x, yend=y.y, weight=weight)
  
  keep_Regions <- my_nodes$Region %>% unique %>% as.character %>% sort
  keep_Pops <- my_nodes$Pop %>% unique %>% as.character %>% sort
  
  my_net_test <- ggplot() +
    geom_point(data=my_nodes, aes(x=x, y=y, color=Pop, pch=Region), size=3) +
    geom_segment(data=my_edges, aes(x=x, xend=xend, y=y, yend=yend),
                 colour="black", alpha=0.5, size=0.5) +
    scale_shape_manual(values=(pop_shapes[["Region"]][keep_Regions,"Shape"])) +
    scale_color_manual(values=as.character(pop_colours[["Pop"]][keep_Pops,"Color"])) +
    ggtitle(title) +
    #scale_alpha_continuous(range = c(0, 1), limits=c(0, 100)) +
    theme_void() + theme(plot.title = element_text(face="bold", hjust=0.5))
  
  return(my_net_test)
}

# EXAMPLE USAGE

network_cutoffs <- list(ibd_1= 1, ibd_5= 5, ibd_10= 10, ibd_25=25, ibd_35=35,
                        ibd_45=45, ibd_50=50, ibd_75=75, ibd_90=90)

poster_networks <- lapply(network_cutoffs, function(x) {
  set_pairwise_threshold(x, IBD_metrics, "fraction_IBD", metadata, population_levels)})

poster_network_plots <- lapply(rev(names(network_cutoffs)), function(x) {
  plot_pairwise_network(poster_networks[[x]], pop_colours,
                        pop_shapes, title=paste0("Isolates sharing IBD% >", network_cutoffs[[x]]))})

View(poster_network_plots) # save the figures .pdf format


# more data analysis 


# histogram of IBD to select highly related parasite pairs 

IBD_metrics_final<-read.csv('IBD_metrics_final.csv', header=T,  sep=",")

hist(IBD_metrics_final$IBD_percent, breaks = 100, type = "percent")
library(lattice)
library(openintro)
hoziHist(IBD_metrics_final$IBD_percent, xlim = c(5,10),
         main = "Pairwise Relatedness of P.faciparum across Zambia",  #Create a title for the chart
         xlab = "Shared IBD%",  #Label the x-axis. Include units!
         col = "darkslategray3")

library(ggplot2)
ggplot(IBD_metrics_final,aes(IBD_percent,y = ..density..))+geom_histogram(bins=100,col = "black", fill="black") +coord_flip() + scale_y_continuous(labels = scales::percent) + 
  labs(x="IBD (%)", y="Frequency (%)") +
  ggtitle("Pairwise Relatedness of P.faciparum across Zambia")

ggplot(IBD_metrics_final,aes(IBD_percent,y = ..density..))+geom_histogram(bins=200, col = "black",fill="black")+ scale_y_continuous(labels = scales::percent) + 
  labs(x="IBD (%)", y="Frequency (%)") +
  theme_bw()+
  theme_classic()+
  xlim(5, 100)+
  ggtitle(" >5% Relatedness of P.faciparum across Zambia")

# ================ |-------------| ================
# ================ |-------------| ================




# ================ |-------------| ================
# ================ |-------------| ================

##################################################
# Spatial distribution  IBD >= 0.25 sharing across Zambia
##################################################

# get plotting data

samples238<- read.csv("metadata_FWS_283samples.csv", header=TRUE, sep=",")

### column names: 

### critical ones sample id, latitude and logitude plus whatever category you want to color by.

#===================================================================#
#           load packages, download maps (GADM)                     #
#===================================================================#
# read in shapefile and .RData formats from Global Administrative Areas (gadm.org), 
library(ggplot2)
library(maptools)
library(ggmap) 
library(mapproj) 
library(raster) 
library(maps)
library(RColorBrewer)

# global administrative district maps-downlod Geopackage and rename as GADM
adm1 <- getData('GADM', country='ZMB', level=0)
adm2 <- getData('GADM', country='ZMB', level=1)
adm3 <- getData('GADM', country='ZMB', level=2)

fadm1 = fortify(adm1)
fadm2 = fortify(adm2)
fadm3 = fortify(adm3)

#===================================================================#
#                       Prepare names of regions                    #
#===================================================================#
### Go to http://gadm.org/country, select Country: "Zambia" and File format: "R (SpatialPolygonsDataFrame)" use Download GADM data (version 2.8) from old file format down RS all 0,1 and 2 levels
### save in the working directory and load as below
GADM <- readRDS("ZMB_adm1.rds")
Zambia.adm1.spdf<-get("GADM") # complex list and download all version gadm36 file 0, 1 and 2...


### Get centroids of spatialPolygonDataFrame and convert to dataframe ----------- 
Zambia.adm1.centroids.df <- data.frame(long=coordinates(Zambia.adm1.spdf)[, 1], lat=coordinates(Zambia.adm1.spdf)[, 2])

### Get names and id numbers corresponding to administrative areas-----------
Zambia.adm1.centroids.df[, 'ID_1'] <- Zambia.adm1.spdf@data[,'ID_1']
Zambia.adm1.centroids.df[, 'NAME_1'] <- Zambia.adm1.spdf@data[,'NAME_1']
save(Zambia.adm1.centroids.df, file = "Zambia.adm1.centroids.df.RData")

load("Zambia.adm1.centroids.df.RData") # Zambia.adm1.centroids.df

#===================================================================#
# Plot samples on Zambia map, color by provinces
#===================================================================#

p <- ggplot(samples283, aes(x=long,y=lat)) +
  geom_polygon(data=fadm2, aes(x=long, y=lat, group = group), fill = "white", colour="grey40", size = 0.2) +
  geom_jitter(data=samples283, position=position_jitter(width=0.2, height=0.2), size = 0.5, color="black",
              # aes(x=lat, y=long),alpha = 6/10)+ #, color=variety, shape=factor(status)
              # aes(x=lat, y=long, alpha = 6/10, color=variety, shape=factor(status))) + ### if you want to color and shape
              aes(x=long, y=lat),alpha=1) +
  #scale_color_brewer(palette = "Set1",type = "qual") + 
  #scale_color_manual(values= c("red","orange", "tomato4","mediumvioletred","darkgreen", "darkblue", "dodgerblue")) +
  geom_text(data = Zambia.adm1.centroids.df, aes(label = NAME_1, x = long, y = lat, group = NAME_1), size = 2)+
  coord_equal(ratio=1) +
  labs(title = " 283 DBS samples spatial distribution", x = "Longitude", y = "Latitude") +theme_bw()
#theme_grey(base_size = 10)
png(filename = "Zambia_283_DBS_Sample_Distribution2.png",width = 8,height = 8,res = 600, units = "in")
p
dev.off()

# Define colors based your data  and plot is proportional to 

##################################################
# Highly related parasites 
################################################## 

# Load Library
library("maps")
library(igraph)
library("geosphere")
library(ggmap)
library(lubridate)
library(data.table)
library(ggrepel)
library(rgdal)

# The data we will use contains PWIBD and Geoinfo among them. 
# The Geoinfo file includes info about latitude and longitude.
# If we did not have those, we could use 'geocode()' from 'ggmap'
# to get latitude and longitude for an address.

Geoinfo<- read.csv("high5precent_NODES1.csv", header=TRUE, sep=",") # # Sample id and Lat and Long data

PWIBD <- read.csv("high5precen_data_EDGES1.csv", header=TRUE, as.is=F) # format your data Source (sample1) vs Target (sample2) and Freq as IBD proportion

head(PWIBD)                                             

head(Geoinfo)

# Plot a map of the zambia:
#ggplot(data=adm1,aes(x=long,y=lat,  border="gray10", fill=TRUE, bg="gray30", group=group)) + geom_path()  #+theme_bw()

plot(adm1)
# Plot a map of the united states:
map(adm2, col="grey20", fill=TRUE, bg="black", lwd=0.5)


# Add a point on the map for each airport:
points(x=Geoinfo$long, y=Geoinfo$lat, pch=19, 
       cex=0.2, col="blue")


# Select at least 1 connections in the data.
tab <- table(PWIBD$Source)
big.id <- names(tab)[tab>0]
Geoinfo <- Geoinfo[Geoinfo$ID %in% big.id,]
PWIBD  <- PWIBD[PWIBD$Source %in% big.id & 
                  PWIBD$Target %in% big.id, ]


# Generate edge colors: lighter color means higher flight volume.
col.1 <- adjustcolor("deepskyblue", alpha=1)
col.2 <- adjustcolor("deepskyblue4", alpha=1)
edge.pal <- colorRampPalette(c(col.1, col.2), alpha = TRUE)
edge.col <- edge.pal(100)

# For each connectivity, we will generate the coordinates of an arc that connects
# its star and end point, using gcIntermediate() from package 'geosphere'.
# Then we will plot that arc over the map using lines().
for(i in 1:nrow(PWIBD))  {
  node1 <- Geoinfo[Geoinfo$ID == PWIBD[i,]$Source,]
  node2 <- Geoinfo[Geoinfo$ID == PWIBD[i,]$Target,]
  
  arc <- gcIntermediate( c(node1[1,]$long, node1[1,]$lat), 
                         c(node2[1,]$long, node2[1,]$lat), 
                         n=1000, addStartEnd=TRUE )
  edge.ind <- round(100*PWIBD[i,]$Freq / max(PWIBD$Freq))
  
  lines(arc, col=edge.col[edge.ind], lwd=edge.ind/30)
}




# ================ |-------------| ================




