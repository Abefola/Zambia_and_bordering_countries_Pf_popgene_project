

####################################################
####################################################
# Fig1B - R script to determine rate of within host multiplicity infection using FWS statistics
# Written by Abebe Fola 05/26/2022    
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
library(moimix)
library(flexmix)
library(lattice)
library(tidyr)
library(ggplot2)
library(cowplot)
library(knitr)
library(kableExtra)
library(dplyr)
library(readr)
library(SeqVarTools)
library(gdsfmt)
library(SeqArray)
library(graph)


#) Set your working directory

setwd("C:/Users/afola/Desktop/Abefola_github/Zambia_2018_MIS_Pf_WGS_project/Fig1")

#) import metadata and please change the file name as needed
# File format: csv 
# 1st column is sample name
# 2nd column is country name 

pop<-read.csv('meta_info_sample_id_province.csv', header=T,  sep=",") 

#) Please change the colors if needed (seven colors I used for provinces)

col.list <- c("red","orange", "tomato4","mediumvioletred","darkgreen", "darkblue", "dodgerblue")

# Import vcf file and change the name of vcf file to the one you want to use

vcf_file_samples238='Pf_capture238_final_MAFfilter_002.recode.vcf' # from Pf_capture241_final_MAFfilter_002.recode.vcf files, We removed 3 samples one from southern, one from Central and one from Lusaka provinces
# bcftools view -S Pf_capture238_samples.tsv Pf_capture241_final_MAFfilter_002.recode.vcf > Pf_capture238_final_MAFfilter_002.recode.vcf


vcf_header <- seqVCF_Header(vcf_file_samples238)

# recode header format for AD 
vcf_header$format$Number[vcf_header$ID == "AD"] <- "."

# info columns to retain
info.import <- c("AC", "AF", "AN", "BaseQRankSum", "DP", "DS",
                 "ExcessHet", "FS", "InbreedingCoeff", "MQ",
                 "MQRankSum", "QD", "ReadPosRankSum", "SOR", "ANN")

# format columns to retain
format.import <- c("AD", "DP", "GQ", "GT", "PL", "RGQ", "SB")

# convert VCF to GDS

GDS_Pf238_samples<-seqVCF2GDS(vcf_file_samples238,
           "GDS_Pf238_samples.gds",
           header=vcf_header, info.import=info.import,
           fmt.import=format.import)

# check summary 
seqSummary(GDS_Pf238_samples)

# save gds file
isolates238samples <- seqOpen(GDS_Pf238_samples)

# check summary 
seqSummary(isolates238samples)

#save sample identifiers

sample.id <- seqGetData(isolates238samples, "sample.id")

# get SNP coordinators 
coords <- getCoordinates(isolates238samples)

head(coords)

# Save SNP coordinates as .csv file and use this file as input file for snp density plot (modify accordingly)
 write.csv(coords, "list_of_snps_for_snp_density_plot.csv")

####################################################
# Calculate FWS value
####################################################

 fws_all_238samples <- getFws(isolates238samples)

 fws_all_238samples<-read.csv('fws_all_238samples.csv', header=T,  sep=",") 

 # Merge fws result with metadata

 fws_metadata <- left_join(pop,fws_all_238samples, by ="sample_id" )

 fws_metadata <- write.csv(fws_metadata, "fws_metadata.csv")# remove 1 sample central, 1 from South and 1 Lusaka
 
 fws_metadata_mod <- read.csv("fws_metadatamod.csv") 

 colors <-c(rep("red", 14), rep("orange", 16), rep("tomato4", 125), rep("mediumvioletred", 17), rep("darkgreen", 16), rep("darkblue", 20),rep("dodgerblue", 30))
 
# Plot  FWS Distribution per province across Zambia
figpdf = paste('Fig1B.pdf', sep="")
pdf(file = figpdf)
plot(fws_metadata_mod$fws, pch=19, cex=1.2, col=colors, ylim=c(0,1), ylab=" FWS", xlab= "Sample", main = " Per Province FWS Distribution")
abline(h=0.95, col="black", lwd=2, lty=2)
legend(x="bottomleft", legend=c("Copperbelt", "Eastern", "Luapula", "Muchinga","North-Western", "Northern", "Western"),
       col=c("red","orange", "tomato4","mediumvioletred","darkgreen", "darkblue", "dodgerblue"), title(xlab = " Per Province FWS Distribution"),
       cex=1, pch=c(19))
dev.off()


#OR

ggplot(fws_metadata_mod, 
       aes( y = fws,
            x = province, color=province)) +
  geom_boxplot(color="black")+
  geom_jitter(position=position_jitter(0.0), col=colors)
  
##################################################
##################################################
# Fig1A - Select only polyclonal samples (fws<=0.95 as polyclonal) and plot per province and visualize on Zambia map
##################################################
##################################################

Onlypolyclonal<- fws_metadata_mod %>% 
  group_by(province) %>% 
  summarise(mono = sum(fws >= 0.95, na.rm = TRUE),
            poly = sum(fws < 0.95, na.rm = TRUE),
            N = mono + poly,
            polyprev = poly / N*100) %>% 
  mutate(ci.lower = polyprev - 1.96 * sqrt(polyprev * (100 - polyprev) / N),
         ci.upper = polyprev + 1.96 * sqrt(polyprev * (100 - polyprev) / N))

ggplot(data = Onlypolyclonal,
       aes (fill=province,
         y = polyprev,
         x = province,
         ymin = ci.lower,
         ymax = ci.upper)) +
  geom_bar(stat = "identity", position = "dodge",col= col.list ) +
  geom_errorbar(position = "dodge") +
  scale_fill_manual(values=col.list)+
labs(x="Province", y="Polyclonal (%)") +
theme(axis.text.x=element_text(size=rel(1.2), angle=90))+
  ggtitle("Polygenonomic infection prevalence ")

# Save the data
write.csv(Onlypolyclonal, "polyclonalprev.csv")

ggsave("polyclonalprev.svg", dpi=600, width=7.5, height=7)
ggsave("polyclonalprev.pdf", dpi=600, width=7.5, height=7)

########################
#Using polyclonalprev.csv data determine spatial distribustion per regions

#######################
#===================================================================#
#           load packages, download maps (GADM)                     #
#===================================================================#
# read in shapefile and .RData formats from Global Administrative Areas (gadm.org), 
library(ggplot2)
library(maptools)
library(rgdal)
library(ggmap) 
library(mapproj) 
library(raster) 
library(maps)
library(sp)
library(RColorBrewer)
library(tidyverse)
library(rnaturalearth) 
library(rnaturalearthdata)
library(sf) #see ubuntu issues here: https://rtask.thinkr.fr/installation-of-r-4-0-on-ubuntu-20-04-lts-and-tips-for-spatial-packages/
library(ggspatial)
library(colorspace)
library(dplyr)

prev <- read.csv("polyclonalprev.csv")

#get admin shape file from naturalearth
admin10 <- ne_download(scale="large", type = "admin_1_states_provinces",
                       category = "cultural", returnclass = "sf")
rivers10 <- ne_download(scale = 10, type = 'rivers_lake_centerlines', 
                        category = 'physical', returnclass = "sf")
lakes10 <- ne_download(scale = "medium", type = 'lakes', 
                       category = 'physical', returnclass = "sf")
oceans10 <- ne_download(scale = "medium", type = "coastline",
                        category = 'physical', returnclass = "sf")
sov110 <- ne_download(scale="medium", type = "sovereignty",
                      category = "cultural", returnclass = "sf")


# global administrative district maps
getwd()
Zamb<-st_read("ZMB_adm0.shp") # Source https://www.diva-gis.org/datadown
Zamb_regions<-st_read("ZMB_adm1.shp")
View(Zamb_regions)
Zamb_districts<-st_read("ZMB_adm2.shp") # name NAME_2 for combining
#View(Zamb_districts)
#Set coordinate reference system for Zambia
Zamb_crs<-st_crs(Zamb)
africa <- ne_countries(scale="medium", type = "sovereignty", continent = "Africa", returnclass = "sf", )


#get admin shape file from naturalearth

polyclonalspatial <- right_join(Zamb_regions, prev, by = "NAME_1")

#incidence

samplessites <-read.csv("metadata_238samples.csv", header=T)
cor.reg <- c("red","orange", "tomato4","mediumvioletred","darkgreen", "darkblue", "dodgerblue") 


pop1 = factor(samplessites$province)

ggplot() + 
  geom_sf(data=africa, fill="white")+
  geom_sf(data = Zamb, fill="white", lwd=3) +
  geom_sf(data = Zamb_regions, color="grey10", cex =10, fill="white") +
  geom_sf(data = polyclonalspatial,aes(geometry=geometry, fill=polyprev
  )) +
  scale_fill_distiller(type="seq", palette = "Greys", direction = 1, name=" Polyclonal \nPrev (%)")+
  geom_sf(data=lakes10, color="grey40", fill ="aliceblue", size= 0.8)  +
  annotation_scale(location = "bl", width_hint = 0.25) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.5, "in"), pad_y = unit(0.2, "in"),
                         style = north_arrow_fancy_orienteering)+
  
  #coord_sf(xlim = c(29, 41), ylim = c(-12, 0), expand = FALSE)+xlab("Longitude")+ ylab("Latitude")+
  # coord_sf(xlim = c(22.3, 33.3), ylim = c(-8,-18), expand = TRUE) +
  coord_sf(xlim = c(22.5, 33.3), ylim = c(-8,-19), expand = TRUE) +
  theme(panel.background = element_rect(fill = 'aliceblue', colour = 'grey97')) +
  theme(axis.ticks = element_blank(),
        
        axis.text.x = element_blank(),
        
        axis.text.y = element_blank()) +
  xlab("") +
  
  ylab("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_jitter(data=samplessites, position=position_jitter(width=0.0, height=0.0), size = 2.5, pch=19, col=cor.reg[as.integer(pop1)], 
              # aes(x=lat, y=long),alpha = 6/10)+ #, color=variety, sha5e=factor(status)
              # aes(x=lat, y=long, alpha = 6/10, color=variety, shape=factor(status))) + ### if you want to color and shape
              aes(x=cluster_longitude, y=cluster_latitude, col= "province"))+
  ggtitle(expression(paste("Regional ", italic("Polyclonal"), "infections prevalence")))
  
# Save the plot
ggsave("Fig1A.pdf", dpi=600, width=7.5, height=7)



##################################################
##################################################
# Fig1C - Correlation b/n Polyclonal infection vs PET-PCR parasite prevalence at cluster level
# For this analysis clusters having at least 3 samples were included 
##################################################
##################################################


# Load the file

Polygenomic_pcr<-read.csv("polyclonal_lm_3samplesatleast.csv", header=T,  sep=",")

#attach(Polygenomic_pcr)

# correlation between Polyclonal infection, Prevalence_at_cluster_level

cor3<-ggplot() + geom_point(data=Polygenomic_pcr, aes(x=PET_PCR_Prevalence, col=province, y=Polygenomic_percent, size= sample_size),  alpha=1)+
  theme_bw()+
  scale_size(range = c(2, 6),
             breaks = c(3, 6, 9, 12),
             labels = c("3-6", "6-9", "9-12", ">12"),
             name= "Sample size",
             guide = "legend")+
  scale_y_continuous(limit =c(20, 100))+ scale_x_continuous(limit =c(10, 70))+labs( x="PET-PCR Prevalence (%)", y="Polygenomic Infection (%)")

cor4<- cor3 + scale_color_manual(values=c("red","orange", "tomato4","mediumvioletred","darkgreen", "darkblue", "dodgerblue"), name= "Province") +
  stat_smooth(method="lm", col="black", size=0.3,  aes(x=PET_PCR_Prevalence, y=Polygenomic_percent)) +
  ggtitle("Correlation between PET-PCR vs Polygenomic Infection ")

cor4
  
###########################
## Adjusted sample size dots and legends 
###########################

cor4+scale_size_area(max_size=8)

# Save the plot
ggsave("Fig1C.pdf", dpi=600, width=7.5, height=7)


# See the details software versions
sessionInfo()
