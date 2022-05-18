###################################################
###############  Chemical Analysis  ###############
###################################################

## Files you'll need...
##### A .csv file with retention time as the columns and sample names as the rows with peak areas in the cells
##### A .csv metadata file with any metadata information you want to list for each sample (species, sex, location, etc...)


## Start with a fresh slate to work on
cat("\014") # Clear console
rm(list=ls()) # Remove all variables to start with a clean slate

## Load required packages
library(vegan)
library(ggplot2)

## Read in the .csv files that you will be working with
CHEM <-read.csv("Methyl Ether Matrix.csv", sep=",", row.names=1)
METADATA <- read.csv("Metadata.csv")

################################################
########  Convert Data to Proportions  #########
################################################
## We first need to convert the peak areas for each sample into "percentage" data, meaning that each peak area is proportional to the other peak areas within a sample
## Make sure to take out any internal standards
## You might also want to prune your data before this step so that you get rid of any contamination or unwanted peak areas


## Check the dimension of the file you just loaded (you will be using this for the next step)
#below function spits out results showing your total number of rows and columns
dim(CHEM) #should be 246 rows, 124 cols
dims <- dim(CHEM) 
## Calculate the relative proportions of the areas of each peak
rsums <- rowSums (CHEM[,1:dims[2]]) #sums up values across all columns of each row
#for the above code "rowSUMS (CHEM[,A:B])", you're indicating to sum each row from column A through column B
CHEMPP <- CHEM[,1:dims[2]]/rsums #calculates relative proportions + puts into a new matrix
dim(CHEMPP) #gives you number of rows, columns in the new matrix

rowSums (CHEMPP[,1:dims[2]]) #verifies that all porportions add up to 1 
sum(rowSums (CHEMPP[,1:dims[2]])) #should be 246 = number of files
#if everything doesn't add up to 1, this probably means that you need to change the A:B situation

## Export the CHEMPP results as a .csv file
write.csv(CHEMPP, file="CHEMPP.csv",row.names=TRUE)

################################################
##########  NMDS Plots Using ggplot2  ##########
################################################

CHEM_NMDS <- metaMDS(CHEMPP, distance="bray", k=2, trymax=100)

MDS1 = CHEM_NMDS$points[,1]
MDS2 = CHEM_NMDS$points[,2]
NMDS = data.frame(MDS1 = MDS1, MDS2 = MDS2, Species = METADATA$Species, Location = METADATA$Location, Sex = METADATA$Sex)

ggplot(NMDS, aes(x=MDS1, y=MDS2, col=Species)) +
  geom_point() +
  stat_ellipse(level=0.65) +
  theme_bw() +
  labs(title = "NMDS Plot")

ggplot(NMDS, aes(x=MDS1, y=MDS2, col=Location)) +
  geom_point() +
  stat_ellipse(level=0.65) +
  theme_bw() +
  labs(title = "NMDS Plot")

ggplot(NMDS, aes(x=MDS1, y=MDS2, col=Sex)) +
  geom_point() +
  stat_ellipse(level=0.65) +
  theme_bw() +
  labs(title = "NMDS Plot")

################################################
##########  Statistical Significance  ##########
################################################

#ANOSIM is a method that tests whether two or more groups of samples are significantly different.

CHEM_DIST <- as.matrix((vegdist(CHEMPP, "bray")))

anosim_results = anosim(CHEM_DIST, METADATA$Species)
anosim_results # take a look at results
summary(anosim_results)
plot(anosim_results)

anosim_results = anosim(CHEM_DIST, METADATA$Location)


#Adonis is a nonparametric statistical method that takes a distance matrix and a category to determine sample grouping from.

adonis_results = adonis(CHEM_DIST ~ Sex, METADATA)
adonis_results # take a look at results; there are no summary() or plot() methods included

###################################################
###################################################
## Notes!
# Vegan package likes data to be formated so that samples = rows, and chemicals = columns


###################################################
###################################################
## Useful links!
# https://jonlefcheck.net/2012/10/24/nmds-tutorial-in-r/
# https://chrischizinski.github.io/rstats/vegan-ggplot2/
# https://ourcodingclub.github.io/2018/05/04/ordination.html
# http://geoffreyzahn.com/nmds-example/
# http://tutorials.iq.harvard.edu/R/Rgraphics/Rgraphics.html




