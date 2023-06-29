
#############################################
#############################################
#############################################
# FigS10 - This script generates a genetic distance and heatmap using fasta file as input file
# Written by Abebe Fola 06/03/2022
#############################################
#############################################

# Set working directory
setwd ("C:/Users/afola/Desktop/Abefola_github/Zambia_and_bordering_countries_Pf_popgene_project/FigS10")

#load the library
library(dplyr)
library(plotly)
library(ape)
library(scatterplot3d)
library(PopGenome)
library(factoextra)

# Input file - VCF file was changed to fasta file using PGDSpider software -http://www.cmpg.unibe.ch/software/PGDSpider/PGDSpider%20manual_vers%202-1-1-5.pdf

GENOME_class <- readData("FASTA")

fasta_238samples="238samples_fasta_final.fa"

# Reads fasta file
dna <- ape::read.dna(fasta_238samples, format="fasta") %>% as.character()
# Calculates distance
dist <- dist.dna(as.DNAbin(dna), model="raw", pairwise.deletion=TRUE)
# Stores into a matrix
hist(dist, col="gray", xlab="Genetic distance", ylab="Frequency", main="Pairwise Genetic Distance between Samples" )

mat <- dist %>% as.matrix()

heatmap.2(mat,Rowv=TRUE,Colv=FALSE,
          
          dendrogram= "row",
          
          labRow=rownames(dist2),
          
          distfun = dist,
          
          hclustfun = hclust,
          
          key=TRUE,
          
          keysize=1,
          
          lhei=c(0.5,4.5),
          
          trace="none",
          
          density.info="none",
          
          margins=c(3,8),
          
          col=rev(redgreen(256)),
          
          main="Pairwise genetic distance")

