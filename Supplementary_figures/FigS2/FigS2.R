####################################################
####################################################
# R script to calcuate number of PET-PCR positive samples seqeunced per Region/province across Zambia
#   Written by Abebe Fola       
####################################################
####################################################
#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install("gdsfmt")
#BiocManager::install("SNPRelate")
#BiocManager::install("SeqArray")
#BiocManager::install("SeqVarTools")
#BiocManager::install("Rgraphviz")
#BiocManager::install("graph")
#BiocManager::install("tess3r")
#BiocManager::install("BiocParallel")
#devtools::install_github("bahlolab/moimix") # R version >R 3.6.2 Required

#Load library 

library(tidyr)
library(ggplot2)
library(dplyr)


## This is the version info of the R packages used by Fola et al. 
#moimix_0.0.2.9001    SeqVarTools_1.24.1   knitr_1.30           graph_1.64.0         flexmix_2.3-17       lattice_0.20-38     
# adegenet_2.1.3       ade4_1.7-16          GenomicRanges_1.38.0 GenomeInfoDb_1.22.1  IRanges_2.20.2       S4Vectors_0.24.4    
# BiocGenerics_0.32.0  SeqArray_1.26.2      vcfR_1.12.0          doMC_1.3.7           iterators_1.0.13     foreach_1.5.1       
# tidyr_1.1.3          isoRelate_0.1.0      devtools_2.3.2       usethis_2.0.0        ggnetwork_0.5.8      ggplot2_3.3.3       
# igraph_1.2.6         dplyr_1.0.5          SNPRelate_1.20.1     gdsfmt_1.22.0     

#) Set your working directory

setwd("C:/Users/afola/Desktop/Abefola_github/Zambia_2018_MIS_Pf_WGS_project/FigS2")


Complete_459samples<- read.csv("Complete_459samples_metadata.csv", header = TRUE)

Complete_459samples %>% 
  group_by(province) %>%
  tally()  %>% 
  arrange(n)


# Number of samples in each Region:



barplotreg <- ggplot(Complete_459samples, aes(province))
barplotreg  + 
  geom_bar() +
  labs(x="Province", y="Number of DBS Samples") +
  theme(axis.text.x = element_text(angle = 90))+
  theme() +
  ggtitle("Samples Distribution per Province") 

ggsave("FigS2.svg", dpi=600, width=7.5, height=7)
ggsave("FigS2.pdf", dpi=600, width=7.5, height=7)

dev.off()
sessionInfo