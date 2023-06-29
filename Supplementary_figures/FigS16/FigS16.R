
####################################################
# FigS16 iHS -whole-genome 
#LINK: https://cran.r-project.org/web/packages/rehh/vignettes/rehh.html
# Written by Abebe Fola 
# 06/26/2026
####################################################

# Load the required Libraries

library(rehh)
library(vcfR)
library(data.table)
library(R.utils)
library(ggplot2)        

# Set working directory 

setwd("C:/Users/afola/Desktop/Zambia_projects/Zambia_2018_WGS_project_ms/Zambia_2018_MIS_Data_Analysis/IHS")

# load vcf file

Pf238_samples<- ("special_tes_allchrom_biallelic_SNPs_maxmix10minQ30minDP10Maf02samplemis90.fixed.recode.vcf.g") # NB three samples where removed from downstream analysis
vcf_header <- seqVCF_Header(Pf238_samples)


chr1_ehh<-rehh::data2haplohh("special_tes_allchrom_biallelic_SNPs_maxmix10minQ30minDP10Maf02samplemis90.fixed.recode.vcf.gz",
                               min_perc_geno.mrk =90,
                               vcf_reader = "data.table", 
                               verbose = T , 
                               polarize_vcf = F,
                               remove_multiple_markers=T,
                               chr.name = "Pf3D7_01_v3")



# Computing genowide scan 

chr1.scan <- scan_hh(chr1_ehh,discard_integration_at_border = F)

write.csv (chr1.scan, "chr1.scan.csv")

# Computing iHS

chr1.ihs<-ihh2ihs(chr1.scan, freqbin = 1, verbose = F, standardize = T)

write.table(x=chr1.ihs$ihs, file = "chr1_ihs.tsv")


# plot iHS 
rehh::manhattanplot(chr1.ihs,
                    pval =  T,
                    threshold = 5,
                    main = "iHS (Zambia_MIS_WGS_Pf - Chr1)")
#Or

rehh::manhattanplot(chr13.ihs,
                    pval =  T,
                    threshold = 5,cex = 1, cex.axis = 0.9, 
                    suggestiveline = F, genomewideline = F, chrlabs = c(1:20,   "P", "Q"),
                    main = "iHS (Zambia_MIS_WGS_Pf - Chr1)")


###########################
###########################
# Chrom 1-14 (each chrom analysis)
# Do the above step for each chrom and combine all chromosomes data set for down stream analysis

###########################
###########################


# All chromosomes analysis

# Combine all ihs.cvs

allchrom.ihs <- read.csv("All_chrom_msmt_tes_iHS.csv")


# Distribution of standardized values: the function distribplot()
##
#Under the assumption that most sites evolve neutrally, the standardized iHS values should follow a normal distribution with the sites under selection as outliers.
###
distribplot(allchrom.ihs$IHS, xlab = "iHS")

#qqplot - use qqplot = T, to determine if data set comes from a population
#   with common distributions 

distribplot(allchrom.ihs$IHS, 
            xlab = "iHS", 
            qqplot = TRUE)

#manhattanplot

rehh::manhattanplot(allchrom.ihs,
                    pval =  T,
                    threshold = 5,
                    main = "iHS (Zambia-MIS - allchroms)")
#Or

rehh::manhattanplot(allchrom.ihs, inset = 0, ylab = "|iHS|",
                    pval =  T,
                    threshold = 5,cex = 1, cex.axis = 0.9, 
                    suggestiveline = F, genomewideline = F, chrlabs = c(1:20,   "P", "Q"),
                    main = "iHS (Zambia-MIS - allchroms)")







