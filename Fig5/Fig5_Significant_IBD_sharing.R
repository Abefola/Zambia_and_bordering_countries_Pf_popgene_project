

###########
# By Abebe Fola
# 03-11-2020
############

###############
###############
# Fig5 - Manhattan plot 
#https://www.r-graph-gallery.com/101_Manhattan_plot.html
###############
###############


# Download required libraries 

library(qqman)
library(calibrate)
library(CMplot)

# Input data (this is output data from IBD analysis - Significant IBD sharing).  Check the details in Fig4_IBD_analysis.r code

samples238_iR<- read.csv("samples238_iR_manhattan_formated.csv", header=TRUE)  # formatted as manhattan plot format 
#(column= chromosome number 1, 2, 3, etc,column 2= bp number, column3 = snp_id like and Pf3D7_11_v3:592805 and last column= P.value. NB you have to change log10(p.value) result from isoRelate to p.value =1/10^log10(p.value))


#plot manhattan plot
manhattan(samples238_iR, col = c("black", "royalblue"), annotateTop= T, suggestiveline = F, genomewideline = F, ylim= c(0,20), chr="CHR", bp="BP", snp="SNP", p="P" , main= "Manhattan plot of the significant IBD sharing")

##qqplot 
#It is a good practice to draw a qqplot from the output. It allows to compare the distribution of the pvalue with an expected distribution by chance. Its realization is straightforward thanks to the qq function:

qq(samples238_iR$P, main="Distribution of the P.value with an expected distribution by chance")

##annotate 
#You probably want to know the name of the SNP of interest: the ones with a high pvalue. You can automatically annotate them using the annotatePval argument:

manhattan(samples238_iR, ylim= c(0,20), annotatePval = 0.01, main= "Manhattan plot of the significant IBD sharing")
