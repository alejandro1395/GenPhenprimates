###############################################
#########MAP CLUSTERS ONTO PHYLOGENY###########
###############################################

#while no phylogeny available, map onto vector of sorted species

##########################################
#####DESCRIPTIVE PLOTS FROM JOE TRAITS####
##########################################
##########################################

library(ggplot2)
library(gridExtra)
library(tibble)
library(data.table)
library(ggrepel)
library(dplyr)
library(tidyr)
library(tidyverse)
library(reshape2)
library(dendextend)
library(ComplexHeatmap)
library("factoextra")
library("circlize")
library("RColorBrewer")
library(ggnewscale)
library(scales)
library(cowplot)
library(plyr)
library(ggdendro)
library(viridis)
library(seabron)




#READ ALL TRAITS FILES (1st STEP)
Galan_Acedo <- read.csv("../data/Phenome/Primate_Traits/INPUT_DATA/Galan-Acedo_2019_V4_MergedValues.csv", header=T)
Brains <- read.table("../data/Phenome/Primate_Traits/INPUT_DATA/DeCasien2019_BrainRegionsSpecies.txt", header=T)
LifeHistory <- read.csv("../data/Phenome/Primate_Traits/INPUT_DATA/VanSchaik_Isler2012_LifeHistoryData.csv", header=T)
MatingSystem <- read.csv("../data/Phenome/Primate_Traits/INPUT_DATA/DeCasien2017_MatingSystem.csv", header=T)
Social <- read.csv("../data/Phenome/Primate_Traits/INPUT_DATA/SocialSytemDecasiens2017_Lupold2019_Schultz2011.csv", header=T)
Dispersal <- read.csv("../data/Phenome/Primate_Traits/INPUT_DATA/Schultz2011DispersalGroups.csv", na.strings=c("","NA"), header=T)
SexDimorph <- read.csv("../data/Phenome/Primate_Traits/INPUT_DATA/Lupold2019_SexualOrnamentationData.csv", header=T)
CraniaMale <- read.csv("../data/Phenome/Primate_Traits/INPUT_DATA/PRIMO_CranialMale.csv", header=T)
CraniaFemale  <- read.csv("../data/Phenome/Primate_Traits/INPUT_DATA/PRIMO_CranialFemale.csv", header=T)
MandibleFemale <- read.csv("../data/Phenome/Primate_Traits/INPUT_DATA/PRIMO_MandibleFemale.csv", na.strings=c("","NA"), header=T)
MandibleMale  <- read.csv("../data/Phenome/Primate_Traits/INPUT_DATA/PRIMO_MandibleMale.csv", na.strings=c("","NA"), header=T)
PAD <- read.csv("../data/Phenome/Primate_Traits/INPUT_DATA/PAD_measurements_subset.tsv", header=T, sep = "\t")
#PAD <- read.table("~/PRIMATES/PRIMATE_TRAITS/INPUT_DATA/PAD_measurements_subset.csv", header=T)
Parasites <-read.table("../data/Phenome/Primate_Traits/INPUT_DATA/Parasites_Primates_single_Host_obs.txt", header=T, sep="\t")
SeqSpec <- read.table("../data/Phenome/Primate_Traits/INPUT_DATA/SequencedSpeciesNoSp.txt", header=T)
GroupNames <- read.csv("../data/Phenome/Primate_Traits/INPUT_DATA/GroupNames.csv", header=T)
NameTable <- read.csv("../data/Phenome/Primate_Traits/INPUT_DATA/Old2NewSpeciesNames2.csv", header=T)


##### #READ ALL TRAITS FILES (2nd STEP)
##### Update old species names from datasets to match the "SpeciesBROAD" nomenclature we are using. I am basing the conversion on the pairings in the "NameTable" dataframe

# Brains
NewNames <- setNames(as.character(NameTable$NewName), NameTable$OldName)
Brains$Taxon <- recode(Brains$Taxon,!!!NewNames)
BrainRegions <- Brains
# Mating System
NameDiffs <- setdiff(MatingSystem$Taxon,(intersect(Galan_Acedo$SpeciesBROAD, MatingSystem$Taxon)))
NewNames <- setNames(as.character(NameTable$NewName), NameTable$OldName)
MatingSystem$Taxon <- recode(MatingSystem$Taxon,!!!NewNames)
# Dispersal
NameDiffs <- setdiff(Dispersal$Taxon,(intersect(Galan_Acedo$SpeciesBROAD, Dispersal$Taxon)))
NewNames <- setNames(as.character(NameTable$NewName), NameTable$OldName)
Dispersal$Taxon <- recode(Dispersal$Taxon,!!!NewNames)
# SexualOrnamentation
NameDiffs <- setdiff(SexDimorph$Species,(intersect(Galan_Acedo$SpeciesBROAD, SexDimorph$Species)))
NewNames <- setNames(as.character(NameTable$NewName), NameTable$OldName)
SexDimorph$Species <- recode(SexDimorph$Species,!!!NewNames)
# SocialSystem
#NameTable <- read.csv("~/PRIMATES/INPUT_DATA/Old2NewSpeciesNames2.csv", header=T)
NameDiffs <- setdiff(Social$SpeciesBROAD,(intersect(Galan_Acedo$SpeciesBROAD, Social$SpeciesBROAD)))
NewNames <- setNames(as.character(NameTable$NewName), NameTable$OldName)
Social$SpeciesBROAD <- recode(Social$SpeciesBROAD,!!!NewNames)
# PAD
#PAD <- PAD[-196,] # Remove duplicate Saimiri sciureus
NameDiffs <- setdiff(PAD$Species,(intersect(Galan_Acedo$SpeciesBROAD, PAD$Species)))
NewNames <- setNames(as.character(NameTable$NewName), NameTable$OldName)
PAD$Species <- recode(PAD$Species,!!!NewNames)
# Parasites
NameDiffs <- setdiff(Parasites$Species,(intersect(Galan_Acedo$SpeciesBROAD, Parasites$Species)))
NewNames <- setNames(as.character(NameTable$NewName), NameTable$OldName)
Parasites$Species <- recode(Parasites$Species,!!!NewNames)


##########GENERATE MERGED ALL TRAITS 3RDT step
MergedTraits <- GroupNames
MergedTraits <- merge(MergedTraits, Galan_Acedo,  by.x="SpeciesBROAD", by.y="SpeciesBROAD", all = T)
MergedTraits <- merge(MergedTraits, MatingSystem, by.x="SpeciesBROAD", by.y="Taxon", all = T)
MergedTraits <- merge(MergedTraits, Social, by.x="SpeciesBROAD", by.y="SpeciesBROAD", all = T)
MergedTraits <- merge(MergedTraits, Dispersal, by.x="SpeciesBROAD", by.y="Taxon", all = T)
MergedTraits <- merge(MergedTraits, LifeHistory, by.x="SpeciesBROAD", by.y="Species", all = T)
MergedTraits <- merge(MergedTraits, SexDimorph, by.x="SpeciesBROAD", by.y="Species", all = T)
MergedTraits <- merge(MergedTraits, BrainRegions, by.x="SpeciesBROAD", by.y="Taxon", all = T)
MergedTraits <- merge(MergedTraits, CraniaMale, by.x="SpeciesBROAD", by.y="SpeciesBROAD", all = T)
MergedTraits <- merge(MergedTraits, CraniaFemale, by.x="SpeciesBROAD", by.y="SpeciesBROAD", all = T)
MergedTraits <- merge(MergedTraits, MandibleFemale, by.x="SpeciesBROAD", by.y="SpeciesBROAD", all = T)
MergedTraits <- merge(MergedTraits, MandibleMale, by.x="SpeciesBROAD", by.y="SpeciesBROAD", all = T)
MergedTraits <- merge(MergedTraits, PAD, by.x="SpeciesBROAD", by.y="Species", all = T)
MergedTraits <- merge(MergedTraits, Parasites, by.x="SpeciesBROAD", by.y="Species", all = T)
MergedTraits <- merge(MergedTraits, SeqSpec, by.x="SpeciesBROAD", by.y="SequencedSpecies", all = T)
MergedTraits <- MergedTraits[order(MergedTraits$GroupName),]


#ADDITIONAL FILTERING FOR MISSING SEQUENCED SPECIES - 4TH STEP
# Remove species which are sequenced but not in trait database
#NoTraitSpecies <- setdiff(SeqSpec$SequencedSpecies,(intersect(Galan_Acedo$SpeciesBROAD, SeqSpec$SequencedSpecies)))
#NoTraitSpecies <- c(NoTraitSpecies, "Chiropotes_israelita")
#SeqSpec <- SeqSpec[ ! SeqSpec$SequencedSpecies %in% NoTraitSpecies,]




######REMOVE SPECIES WITH MISSING GROUP INFO
InformativeMergedTraits <- subset(MergedTraits, !is.na(GroupName))
write.table(InformativeMergedTraits, "MergedTraits.tsv", row.names = FALSE, sep = "\t")
ALLTraits <- InformativeMergedTraits



#before calling MCA algorithm for factor analysis of traitsxspecies,
#let's transform the data into presence/absence of traits

ALLTraitsPresence <- ALLTraits %>%
  mutate_all(funs(ifelse(is.na(.), 0, 1)))

#Load FactorMineR package
library(FactoMineR)

as.character(ALLTraits$SpeciesBROAD)
ncol(subset_traits)
subset_traits <- ALLTraitsPresence[,-which(names(ALLTraitsPresence) %in% c("SpeciesBROAD", "GroupName"))]
i=0
while(i < ncol(subset_traits)){
  i=i+1  
  subset_traits[,i] = as.factor(subset_traits[,i])
}
row.names(subset_traits) <- as.character(ALLTraits$SpeciesBROAD)
summary(subset_traits)
res.mca = MCA(subset_traits)





#LET'S STUDY PATTERNS OF CLUSTERED SPECIES onto phylogeny

sugg_spc_res.hcpc <- HCPC(res.mca, nb.clust = 10, graph = FALSE)
View(head(sugg_spc_res.hcpc$data.clust))
species_clust_num <- as.data.frame(cbind(rownames(sugg_spc_res.hcpc$data.clust), sugg_spc_res.hcpc$data.clust$clust))
colnames(species_clust_num) <- c("species", "clusters")
species_clust_num$species <- factor(species_clust_num$species, levels = species_clust_num$species)

#spineplot(factor(species_clust_num$species)~factor(species_clust_num$clusters), 
 #         col = factor(species_clust_num$clusters))

#IDEA1

palette <- brewer.pal(n = 10, name = "Paired")

ggplot(species_clust_num, aes(x = 1,y = species, fill = factor(clusters)), 
stat="identity") + 
  geom_tile() + scale_fill_manual(values=palette) + theme(axis.text.y = element_text(size = 6)) +
  labs(fill = "Species clusters")
                                                          

