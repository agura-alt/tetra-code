library(ape)
library(phangorn)
library(vegan)
library(ggplot2)
library(dplyr)
library(phytools)
library(rgl)
library(stringr)
library(evoldiver)

# Now in color
CHEM <- read.csv("Methyl Ether Matrix.csv")

CHEM[CHEM > 1] <- 1 # CONVERTS CHEM INTO A PRESENCE ABSENCE MATRIX

CHEM

write.csv(CHEM, file="CHEM.csv")

CHEM <- read.csv("CHEM.csv")

## Check the dimension of the file you just loaded (you will be using this for the next step)
#below function spits out results showing your total number of rows and columns
dims <- dim(CHEM) #246 files, 124 MEs + 2 filename column X
CHEM
dims

CHEM <- CHEM[,3:dims[2]] #remove extraneous rownames columns
dims <- dim(CHEM) #reset dims: 246 files, 124 MEs
CHEM
dims

## Calculate the relative proportions of the areas of each peak
rsums <- rowSums(CHEM[,1:dims[2]]) #sums up values across columns of each row

#for the above code "rowSUMS (CHEM[,A:B])", you're indicating to sum each row from column A through column B
CHEMPP <- CHEM[,1:dims[2]]/rsums #calculates relative proportions + puts into a new matrix
CHEMPP

#gives you number of rows, columns in the new matrix
dimspp <- dim(CHEMPP) #246 files, 124 MEs
dimspp

rowSums(CHEMPP[,1:dimspp[2]]) #verifies that all porportions add up to 1 
#if everything doesn't add up to 1, this probably means that you need to change the A:B situation

## Clean out the ".D.CSV_Area" stuff so you have clean sample names
## Export the CHEMPP results as a .csv file
write.csv(CHEMPP, file="CHEMPP.csv",row.names=TRUE)

CHEMPP <- read.csv("CHEMPP.csv", sep=",", row.names=1)

CHEM_NMDS <- metaMDS(CHEMPP, distance="bray", k=2, trymax=100) #runs multiple NMDS to find a stable solution

METADATA <- read.csv("Metadata.csv") # Reads in Metadata
length(METADATA$Sex)
length(METADATA$Species)


MDS1 <- CHEM_NMDS$points[,1]
MDS2 <- CHEM_NMDS$points[,2]
NMDS <- data.frame(MDS1 = MDS1, MDS2 = MDS2, Species = METADATA$Species)
NMDS$Species <- as.factor(NMDS$Species)

NMDS

ggplot(NMDS, aes(x=MDS1, y=MDS2, col=Species)) +
  geom_point() +
  stat_ellipse(level=0.65) +
  theme_bw() +
  labs(title = "NMDS Plot (Presence/Abundance)")

cent <- aggregate(cbind(MDS1,MDS2)~Species,NMDS,mean)
row.names(cent) <- cent$Species
cent <- cent[-1]

species_tree <- read.tree(file = "phylogeny_pruned_without_palu")

#map full species names to shortened names
name_dict <- c("Quasimodo"="Quasi", "Brevignatha"="Brevi", "Waikamoi"="Waika", "Kamakou"="Kama",
               "Restricta"="Rest", "Pallescens"="Pall", "Eurychasma"="Eury", "Stelarobusta"="Stella",
               "Trituberculata"="Tri", "Filiciphilia"="Fili", "ElongateForest"="EF", "Acuta"="Acuta",
               "Versicolor"="Vers")

shorten_names <- function(x) {
  name_dict[x]
}

species_tree$tip.label <- shorten_names(species_tree$tip.label)

#reorder centroids
cent <- slice(cent,match(species_tree$tip.label,row.names(cent)))
dim(cent)

plot(species_tree)

edges <- species_tree$edge
edge_lengths <- (species_tree$edge.length)*7 #scale edge lengths


colour <- c(hsv(edge_lengths[3],1,1),hsv(edge_lengths[7],1,1),hsv(edge_lengths[8],1,1),
           hsv(edge_lengths[9],1,1),hsv(edge_lengths[10],1,1),hsv(edge_lengths[12],1,1),
           hsv(edge_lengths[15],1,1),hsv(edge_lengths[17],1,1),hsv(edge_lengths[18],1,1),
           hsv(edge_lengths[20],1,1),hsv(edge_lengths[22],1,1),hsv(edge_lengths[23],1,1),
           hsv(edge_lengths[24],1,1))



yy <- chronoPTS2D(x=cent$MDS1,y=cent$MDS2,tree = species_tree, radius = 0.05, tip.label = "tip",
                 tip.adj = 1.7, tip.cex = 0.8,node.cex = 0, col = colour, clade = colour)

play3d(spin3d(axis = c(0,0,1),rpm = 30), duration = 7)
#spin3d(yy)
