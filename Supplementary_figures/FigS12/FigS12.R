
#############################################
#############################################
#############################################
# FigS12 - Mantel test (correlation between FST and geographic distance)
# Written by Abebe Fola 06/03/2022
#############################################
#############################################
#############################################

# Load library
library(vegan)
library(ape)
library(ggpubr)


# Set working directory 

setwd("C:/Users/afola/Desktop/Abefola_github/Zambia_and_bordering_countries_Pf_popgene_project/FigS12")


# input files 
geo<- read.csv("istance_matrix_mantelinputfile.csv") # pairwise fst value and geographic distance in km

geo <- geo %>% as.matrix()

gen<- read.csv("ST_mantelinputfile.csv") # pairwise fst value and geographic distance in km

gen <- gen %>% as.matrix()

# Mante test 
mantel.test(geo, gen, graph = TRUE,
            main = "Mantel test: a random example with 6 X 6 matrices
representing asymmetric relationships",
            xlab = "z-statistic", ylab = "Density",
            sub = "The vertical line shows the observed z-statistic")

# Plot the results 

plot(geo,gen,pch=16,cex=1,col="grey",bty="l", ylab="Fst", xlab="Spatial distance in Km",
     main="Correlation between FST and Geographic distance in Km")

abline(lm(gen~geo),col="blue")


ggscatter(data1, x = "Distance_km", y = "FST",
          color = "black", shape = 21, size = 3, # Points color, shape and size
          add = "reg.line",  # Add regressin line
          add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
          conf.int = F, # Add confidence interval
          ylab="FST", xlab="Spatial distance in km",
          cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
          main="Correlation between Genetic distance (FST) and Geographic distance in Km",
          cor.coeff.args = list(method = "spearman", label.x = 0.5)
)


ggscatter(data1, x = "Distance_km", y = "FST",
          add = "reg.line", conf.int = F, 
          cor.coef = TRUE, cor.method = "spearman",
          ylab="FST", xlab="Spatial distance in km",
          main="Correlation between Genetic distance (FST) and Geographic distance in Km")


mantel(geo, gen, method="pearson", permutations=999)


lines(smooth.spline(geo, gen), col="grey", cel=2)


lines(lowess(geo, gen, f=0.33))

abline(lm(geo~gen), col="red")

legend(x="topleft", legend=c("Mantel statistic r = 0.328", "P.value = 0.123"))

ggscatter(dat1, x = "Km", y = "Ms_Fst",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman",
          ylab="Genetic distance", xlab="Spatial distance", col="grey", 
          main="Correlation between Genetic distance (Fst) and Geographic distance in Km")