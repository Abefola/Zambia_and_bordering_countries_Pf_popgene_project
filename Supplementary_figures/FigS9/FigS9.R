
#################
#################
# 11-12-2022
# Written by Abebe Fola 
#################
#################

#################
#################
# FigS9 - PET-pcr vs polygenomic infections
#################
#################

# Set working directory

setwd("C:/Users/afola/Desktop/Abefola_github/Zambia_2018_MIS_Pf_WGS_project/FigS9") # change this to your working directory

# Install Library

library(ggplot2)
library(tidyverse)
library(dplyr)
library(tidyr)


# Load the data 

FigS9_file <- read.csv("polyvspetprc.csv", header=TRUE, sep="," )

cbp2 <- c( "dodgerblue4","grey25")

# plot
ggplot(FigS9_file,aes(x=Province, y=prev, fill = category))+
  geom_bar(stat = "identity", position = "dodge")+
  labs(x="Region", y="Prevalence (%)") +
  scale_fill_discrete(name="Region") +
  theme(axis.text.x=element_text(size=rel(1.2), angle=270))+
  #theme_bw() +
  theme(text = element_text(size = 16)) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  #scale_fill_jco()+
scale_fill_manual(values=cbp2)+
ggtitle("prevalence and polygenomic")+
theme_bw()

# Save plot 
ggsave("FigS9.svg", dpi=600, width=7.5, height=7)
ggsave("FigS9.pdf", dpi=600, width=7.5, height=7)


