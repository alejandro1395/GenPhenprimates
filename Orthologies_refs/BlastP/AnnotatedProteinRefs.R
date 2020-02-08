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
PrimatesTree <- read.nexus("consensusTree_10kTrees_Primates_Version3_withCladeCredibilityValues.nex")
Primates_ggtree <- ggtree(PrimatesTree)

#READ THE TFILE WITH COUNTS

total_npeps <- read.csv("TotalPepPrimates.csv", header = FALSE, sep = "\t")
redundantspecies99_npeps <- read.csv("nr99PepSinglePrimates.csv", header = FALSE, sep = "\t")
summary(total_npeps)
summary(redundantspecies99_npeps)

#Get vectors from count files
species_total <- total_npeps$V1[seq(2,length(total_npeps$V1), by=2)]
species_total <- gsub(" processed", "", species_total)
npeps_total <- total_npeps$V1[seq(1,length(total_npeps$V1), by=2)]
#same but for 99 species
species_99 <- redundantspecies99_npeps$V1[seq(2,length(redundantspecies99_npeps$V1), by=2)]
species_99total <- c()
species99_list <- strsplit(as.character(species_99), "/")
for (element in 1:length(species99_list)){
  species_99total[element] <- species99_list[[element]][1]
}
npeps_99 <- redundantspecies99_npeps$V1[seq(1,length(redundantspecies99_npeps$V1), by=2)]

#CORRECT MISSING names from tree species
species_total[species_total == "Cercocebus_atys"] <- "Cercocebus_torquatus_atys"
species_total[species_total == "Hoolock_hoolock"] <- "Bunopithecus_hoolock"
species_total[species_total == "Prolemur_simus"] <- "Hapalemur_simus"
species_total[species_total == "Sapajus_apella"] <- "Cebus_apella"
species_total[species_total == "Pan_troglodytes"] <- "Pan_troglodytes_troglodytes"


#build dataframes
total_dataset <- data.frame(cbind(as.character(species_total), as.numeric(as.character(npeps_total)), as.numeric(as.character(npeps_99))))
colnames(total_dataset) <- c("species", "npeps_total", "npeps_99")
sum(as.numeric(as.character(total_dataset$npeps_total)))

total_dataset <- total_dataset[-which(total_dataset$species == "Rhinopithecus_strykeri"),]

#melt dataframe
melted_total <- melt(data = total_dataset, id.vars = "species", measure.vars = c("npeps_total", "npeps_99"))

#CREATE PLOT FOR TOTAL

p1 <- facet_plot(Primates_ggtree, panel = 'Total nº of annotated peps', data = total_dataset, 
           geom = geom_barh, 
           mapping = aes(x = as.numeric(as.character(total_dataset$npeps_total))), 
           stat='identity', col="black", fill ="gold") + theme_tree2() + geom_tiplab() + xlim_tree(1) + xlim_expand(c(0, 0.1), 'Tree')


#include stacked for reduced redundancy

p2 <- facet_plot(Primates_ggtree, panel = 'Reduction with 99% redund', data = melted_total, 
           geom = geom_barh, 
           mapping = aes(x = as.numeric(as.character(value)), fill = as.factor(variable)), 
           stat='identity', col="black", position = "identity") + theme_tree2() +
  scale_fill_manual(values = c("gold1", "red3")) + geom_tiplab() + xlim_tree(1) + xlim_expand(c(0, 0.1), 'Tree')


#SUM TOTALS


