#############################################
#############################################
#############################################
# FigS14 - IBD Segements by Centimorgan
# See the Fig4 code how to generate input data used for this analysis 
# Write by Abebe Fola
# Date 06-06-2023
#############################################
#############################################
#############################################


# load libraries 

library(isoRelate)
library(network)
library(dplyr)
library(igraph)
library(ggnetwork)
library(ggplot2)
library(ape)
library(stats)
library(tibble)

# Set working directory

setwd ("C:/Users/afola/Desktop/Abefola_github/Zambia_and_bordering_countries_Pf_popgene_project/FigS14")

# Input data 
IBDC_CM<- read.csv("samples238_final_ibd_segments_final.csv")



# Plot IBD per Centimoran 

summary(IBDC_CM$length_M*100)

Ibdcm<- density(IBDC_CM$length_M*100)

plot(Ibdcm, xlab="Length of segement (cM)", main = "Distribution of length of IBD segments across genome")

polygon(Ibdcm,  col="orange", border="black",lty = 1, lwd = 1,angle = c(-45, 45), bg= "red", ylab="Frequency",main="Paiwise allele sharing")
abline(h=10, col="red")

