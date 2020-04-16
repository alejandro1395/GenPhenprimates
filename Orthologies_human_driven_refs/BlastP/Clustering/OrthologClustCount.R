##########################################
########DESCRIPTIVE PLOTS FROM ORTHOLOG###
############QUALITIES#####################
##########################################

library(phytools)
library(ggplot2)
library(ggnewscale)
library(ggtree)
packageVersion("phytools")
library(ggstance)
library(reshape2)

#READ THE SPECIES TREE
PrimatesTree <- read.nexus("~/Downloads/consensusTree_10kTrees_Primates_Version3_withCladeCredibilityValues.nex")
Primates_ggtree <- ggtree(PrimatesTree)

#READ THE TFILE WITH COUNTS

total_orthos <- read.csv("../summary_counts", header = FALSE, sep = "\t")
summary(total_orthos)

#Get vectors from count files
species_total_raw <- total_orthos$V1[seq(2,length(total_orthos$V1), by=2)]
species_total_raw <- gsub(" processed", "", species_total_raw)
species_total <- c()
species_total_raw <- strsplit(as.character(species_total_raw), "/")
for (element in 1:length(species_total_raw)){
  species_total[element] <- species_total_raw[[element]][1]
}
northos_total <- total_orthos$V1[seq(1,length(total_orthos$V1), by=2)]
sum(as.numeric(as.character(northos_total)))


#CORRECT MISSING names from tree species
species_total[species_total == "Cercocebus_atys"] <- "Cercocebus_torquatus_atys"
species_total[species_total == "Hoolock_hoolock"] <- "Bunopithecus_hoolock"
species_total[species_total == "Prolemur_simus"] <- "Hapalemur_simus"
species_total[species_total == "Sapajus_apella"] <- "Cebus_apella"
species_total[species_total == "Pan_troglodytes"] <- "Pan_troglodytes_troglodytes"
species_total[species_total == "Carlito_syrichta"] <- "Tarsius_syrichta"


#build dataframes
total_dataset <- data.frame(cbind(as.character(species_total), 
                                  as.numeric(as.character(northos_total))))
colnames(total_dataset) <- c("species", "northos_total")
sum(as.numeric(as.character(total_dataset$northos_total)))

total_dataset <- total_dataset[-which(total_dataset$species == "Rhinopithecus_strykeri"),]


#CREATE PLOT FOR TOTAL
Primates_ggtree  <- Primates_ggtree + geom_tiplab(align = TRUE) 
#reorder databy aligning to tiplabels
tree_subset <- Primates_ggtree$data[1:length(total_dataset$species),]
ordered_species <- tree_subset$label[order(tree_subset$y, decreasing = TRUE)]

#plot
total_dataset <- total_dataset[order(match(total_dataset$species, ordered_species), decreasing = TRUE),]
Primates_ggtree + geom_facet(panel = 'Total nº of ortholoygy clusters', data = total_dataset, 
                             geom = ggstance::geom_barh,
                            aes(x = as.numeric(as.character(total_dataset$northos_total))), 
                            stat='identity', col="black", fill ="gold") + theme_tree2() + xlim_tree(1) 
 
 Primates_ggtree <- Primates_ggtree + geom_tiplab(align=TRUE)
 View(Primates_ggtree$data)
  theme_tree2() + xlim_tree(1) 

#include stacked for reduced redundancy

p2 <- facet_plot(Primates_ggtree, panel = 'Reduction with 99% redund', data = melted_total, 
           geom = geom_barh, 
           mapping = aes(x = as.numeric(as.character(value)), fill = as.factor(variable)), 
           stat='identity', col="black", position = "identity") + theme_tree2() +
  scale_fill_manual(values = c("gold1", "red3")) + geom_tiplab() + xlim_tree(1) + xlim_expand(c(0, 0.1), 'Tree')


url <- paste0("https://raw.githubusercontent.com/TreeViz/", "metastyle/master/design/viz_targets_exercise/") 
info <- read.csv(paste0(url, "tip_data.csv"))
info
