
#############################################
#############################################
#############################################
# FigS7 - SNP density plot
# Written by Abebe Fola 06/03/2022
#############################################
#############################################
#############################################


#Plot SNP-density from vcf files 

# load the following libraries 
library(qqman)
library(calibrate)
library(CMplot)

# Set working directory 

setwd("C:/Users/afola/Desktop/Abefola_github/Zambia_2018_MIS_Pf_WGS_project/FigS7")  

# Load SNP density formatted data 

samples238_snp_density<- read.csv("snp_density_formated_data.csv", header = T) # Look the attached file. This one of intermediate file for selection analysis. It contains four columns Snp_coordinate, Chro_no, SNP_position and P_value (not required for SNP density plot)

#plot SNP-density plot

CMplot(samples238_snp_density,plot.type="d",bin.size=25000,chr.den.col=c("darkgreen", "yellow", "red"),file="jpg",memo="",dpi=300,bin.range = c(1, 60), 
       file.output=TRUE,verbose=TRUE,width=9,height=6)          


