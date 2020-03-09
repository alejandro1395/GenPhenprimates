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










###########################
###########################
###1) STUDY OF CATEGORIES##
###########################

#First let's do the first Component analysis of our matrix data to know
#the variables closer together in the axis
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





#LET'S STUDY PATTERNS OF CLUSTERED SPECIES and CLUSTERED TRAITS together HCPC values

sugg_spc_res.hcpc <- HCPC(res.mca, nb.clust = 10, graph = FALSE)

#plot(sugg_spc_res.hcpc, choice = "map")



#convert to dataframe
#best tested method is melt
#TRAITS IN SPECIES CLUSTERS
#sugg_spc_res.hcpc$desc.axes
#plot(sugg_spc_res.hcpc)

traits_by_spc_clust_raw <- melt(sugg_spc_res.hcpc$desc.var$category)
traits_by_spc_clust_raw <- traits_by_spc_clust_raw[which(traits_by_spc_clust_raw$Var2 == "p.value") ,]
variables <- unlist(strsplit(as.vector(traits_by_spc_clust_raw$Var1), "="))
Traits <- variables[seq(1, length(variables), by=2)]
Categories <- variables[seq(2, length(variables), by=2)]
traits_by_spc_clust <- as.data.frame(cbind(traits_by_spc_clust_raw$L1, 
                                           as.character(Traits),
                                           traits_by_spc_clust_raw$value))
variables <- unlist(strsplit(as.vector(traits_by_spc_clust_raw$Var1), "="))
colnames(traits_by_spc_clust) <- c("Cluster", "Trait", "p-value")
traits_by_spc_clust <- traits_by_spc_clust[duplicated(traits_by_spc_clust), ]
#View(traits_by_spc_clust)_raw)




#SPECIES IN SPECIES CLUSTERS

#The easiest idea here is to build a heatmap based on group belonging
#in this case, cluster belonging (and then compare with a clop)
species_by_spc_clust_raw <- as.data.frame(as.matrix(unlist(sugg_spc_res.hcpc$desc.ind$para)))
clusters_spcs <- as.data.frame(unlist(strsplit(row.names(species_by_spc_clust_raw), split='.', 
                                               fixed=TRUE)))
colnames(clusters_spcs) = "all_together"
clusters <- clusters_spcs$all_together[seq(1, length(clusters_spcs$all_together), by=2)]
species <- clusters_spcs$all_together[seq(2, length(clusters_spcs$all_together), by=2)]
clusters_spcs <- as.data.frame(cbind(clusters, as.character(species), species_by_spc_clust_raw$V1))
colnames(clusters_spcs) = c("Cluster", "Species", "coord")
#View(clusters_spcs)



#SPECIES IN TRAITS CLUSTERS

transposed_traits <- as.data.frame(t(subset_traits))
inv.res.mca = MCA(transposed_traits)
sugg_trait_res.hcpc <- HCPC(inv.res.mca, nb.clust = 10, graph = FALSE)

species_by_trait_clust_raw <- melt(sugg_trait_res.hcpc$desc.var$category)
species_by_trait_clust_raw <- species_by_trait_clust_raw[which(species_by_trait_clust_raw$Var2 == "p.value") ,]
variables <- unlist(strsplit(as.vector(species_by_trait_clust_raw$Var1), "="))
Species <- variables[seq(1, length(variables), by=2)]
Categories <- variables[seq(2, length(variables), by=2)]
species_by_trait_clust<- as.data.frame(cbind(species_by_trait_clust_raw$L1, 
                                              as.character(Species),
                                              species_by_trait_clust_raw$value))
colnames(species_by_trait_clust) <- c("Cluster", "Species", "p-value")
species_by_trait_clust <- species_by_trait_clust[duplicated(species_by_trait_clust), ]
#View(species_by_trait_clust)_raw)





#TRAITS IN SPECIES CLUSTERS

#The easiest idea here is to build a heatmap based on group belonging
#in this case, cluster belonging (and then compare with a clop)
traits_by_trait_clust_raw <- as.data.frame(as.matrix(unlist(sugg_trait_res.hcpc$desc.ind$para)))
clusters_traits_raw_list <- strsplit(row.names(traits_by_trait_clust_raw), split='.', 
                                     fixed=TRUE)
clusters <- unlist(lapply(clusters_traits_raw_list, `[[`, 1))#get cluster numbers
traits <- gsub("^[0-9]+.","",row.names(traits_by_trait_clust_raw))

clusters_traits <- as.data.frame(cbind(clusters, as.character(traits), traits_by_trait_clust_raw$V1))
colnames(clusters_traits) = c("Cluster", "Trait", "coord")
View(clusters_traits)







##############################
##############################
###PREPARE DATA FOR MERGED PLOT######
##############################
##############################

species_clusters_assigned <- unique(as.data.frame(cbind(sugg_spc_res.hcpc$data.clust$clust, rownames(sugg_spc_res.hcpc$data.clust))))
colnames(species_clusters_assigned) <- c("Cluster", "Species")
traits_clusters_assigned <- unique(as.data.frame(cbind(sugg_trait_res.hcpc$data.clust$clust, rownames(sugg_trait_res.hcpc$data.clust))))
colnames(traits_clusters_assigned) <- c("Cluster", "Trait")

spc_data <- as.matrix(sugg_spc_res.hcpc$data.clust)
trait_data <- as.matrix(sugg_trait_res.hcpc$data.clust)
Species_names <- unique(rownames(spc_data)[match(labels(sugg_spc_res.hcpc$call$t$tree), rownames(spc_data))])
Trait_names <-  unique(rownames(trait_data)[match(labels(sugg_trait_res.hcpc$call$t$tree), rownames(trait_data))])
total_df <- expand.grid(Species_names, Trait_names)
colnames(total_df) <- c("Species", "Trait")



#fill 2 matrices with log values both for species clusters and for trait clusters
plot_spc_clustmatrix <-  matrix(nrow = length(Species_names), 
                                ncol=length(Trait_names),
                                byrow = TRUE, 
                                dimnames = list(Species_names,Trait_names))
plot_traits_clustmatrix <-  matrix(nrow = length(Species_names), 
                                ncol=length(Trait_names),
                                byrow = TRUE, 
                                dimnames = list(Species_names,Trait_names))
spcs_clusts <- c()
trait_clusts <- c()


#Matrix for species
for(j in 1:ncol(plot_spc_clustmatrix)){
  for(i in 1:nrow(plot_spc_clustmatrix)){
    #look for the significance association values in individual dfs
    current_clast <- as.numeric(as.character(species_clusters_assigned$Cluster[species_clusters_assigned$Species==Species_names[i]]))
    if (Species_names[i] %in% clusters_spcs$Species) {
      species_effect <- as.numeric(as.character(clusters_spcs$coord[as.character(clusters_spcs$Species)==Species_names[i]]))
    }
    else{
      species_effect = 1
    }
    if (Trait_names[j] %in% traits_by_spc_clust$Trait[as.numeric(as.character(traits_by_spc_clust$Cluster))==current_clast]) {
      current_clast <- as.numeric(as.character(species_clusters_assigned$Cluster[species_clusters_assigned$Species==Species_names[i]]))
      trait_effect <- as.numeric(as.character(traits_by_spc_clust$`p-value`[traits_by_spc_clust$Trait==Trait_names[j] &
                                                                              as.numeric(as.character(traits_by_spc_clust$Cluster))==current_clast]))
    }
    else{
      trait_effect = 1
    }
    plot_spc_clustmatrix[i, j] = log(species_effect) + log(trait_effect)
  }}

#Matrix for traits
for(j in 1:ncol(plot_traits_clustmatrix)){
  for(i in 1:nrow(plot_traits_clustmatrix)){
    #look for the significance association values in individual dfs
    current_clast_trait <- as.numeric(as.character(traits_clusters_assigned$Cluster[as.character(traits_clusters_assigned$Trait)==as.character(Trait_names[j])]))
    if (Trait_names[j] %in% clusters_traits$Trait) {
      trait_effect <- as.numeric(as.character(clusters_traits$coord[as.character(clusters_traits$Trait)==Trait_names[j]]))
    }
    else{
      trait_effect <- 1
    }
    if (Species_names[i] %in% species_by_trait_clust$Species[as.numeric(as.character(species_by_trait_clust$Cluster))==current_clast_trait]) {
      current_clast_trait <- as.numeric(as.character(traits_clusters_assigned$Cluster[as.character(traits_clusters_assigned$Trait)==Trait_names[j]]))
      species_effect <- as.numeric(as.character(species_by_trait_clust$`p-value`[as.character(species_by_trait_clust$Species)==as.character(Species_names[i]) &
                                                                                   as.numeric(as.character(species_by_trait_clust$Cluster))==current_clast_trait]))
    }
    else{
      species_effect <- 1
    }
    plot_traits_clustmatrix[i, j] = log(species_effect) + log(trait_effect)
  }}


#convert matrixes filled in dataframe for ggplot2 tile heatmap overlaying colors
#FINAL DRAW OF PLOT#######################
df_species <- as.data.frame(melt(-plot_spc_clustmatrix))
df_traits <- as.data.frame(melt(-plot_traits_clustmatrix))
df_species[sapply(df_species, is.infinite)] <- 1000
df_traits[sapply(df_traits, is.infinite)] <- 1000
df_species$value2 <- df_traits$value






##############################
##############################
###MAKE MERGED PLOT######
##############################
##############################
#make cuts
df_species$value_cut <- cut(df_species$value, breaks = c(-100000, 0 , 100000))
df_species$value_cut2 <- cut(df_species$value2, breaks = c(-100000, 0 , 100000))
df_species$cuts <- paste(df_species$value_cut, df_species$value_cut2, sep = "-")

levels_comb <- expand.grid(lev1 = levels(df_species$value_cut), lev2 = levels(df_species$value_cut))
levels_comb$cuts <- paste(levels_comb$lev1, levels_comb$lev2, sep = "-")
levels_comb$filling <- c("#e8e8e8","#73ae80","#6c83b5","#2a5a5b")
colnames(levels_comb) <- c("Species_cluster", "Traits_cluster", "cuts", "value")
data_m <- left_join(df_species, levels_comb, by = "cuts")
colnames(data_m) <- c("species", "traits", "value.x", "value2", "value_cut", "value_cut2",
                      "cuts", "Species_cluster", "Traits_cluster", "value.y")


#NOW LET'S MAKE THE PLOT BY PARTS
#A) PREPARING DENDOGRAM DATA

#ADD TREES OF CLUSTERS

colorsRainbow <- rainbow(10, alpha = 1)
dendogram_species <- as.dendrogram(sugg_spc_res.hcpc$call$t$tree) %>% 
  set("branches_k_color", k = 10) %>% 
  set("branches_lwd", 0.7)
ggdend_species <- as.ggdend(dendogram_species)
dendogram_traits <- as.dendrogram(sugg_trait_res.hcpc$call$t$tree) %>% 
  set("branches_k_color", k = 10) %>% 
  set("branches_lwd", 0.7)
ggdend_traits <- as.ggdend(dendogram_traits)



# Setup the data, so that the layout is inverted (this is more 
# "clear" than simply using coord_flip())
segment_data_species <- with(
  segment(ggdend_species), 
  data.frame(x = y, y = x, xend = yend, yend = xend, col = col))
#same for traits
segment_data_traits <- with(
  segment(ggdend_traits), 
  data.frame(x = y, y = x, xend = yend, yend = xend, col = col))
# Use the dendrogram label data to position the gene (species)
gene_pos_table <- with(
  ggdend_species$labels, 
  data.frame(y_center = x, species = as.character(label), height = 1))

# Use the dendrogram label data to position the samples (traits)
sample_pos_table <- with(
  ggdend_traits$labels, 
  data.frame(x_center = x, traits = as.character(label), width = 1))

# Neglecting the gap parameters
heatmap_data <- data_m %>% 
  left_join(gene_pos_table) %>%
  left_join(sample_pos_table)

# Limits for the vertical axes
species_axis_limits <- with(
  gene_pos_table, 
  c(min(y_center - 0.5 * height), max(y_center + 0.5 * height))
) + 
  0.1 * c(-1, 1) # extra spacing: 0.1

# Limits for the horizontal axes
traits_axis_limits <- with(
  sample_pos_table,
    c(min(x_center - 0.5 * width), max(x_center + 0.5 * width))) + 
  0.1 * c(-1, 1) # extra spacing: 0.1

# Heatmap plot
plt_hmap <- ggplot(heatmap_data, 
                   aes(x = x_center, y = y_center, fill = value.y, 
                       height = height, width = width)) + 
  geom_tile() +
  scale_x_continuous(breaks = sample_pos_table[, "x_center"], 
                     labels = sample_pos_table$traits, 
                     limits = traits_axis_limits, 
                     expand = c(0, 0)) + 
  # For the y axis, alternatively set the labels as: gene_position_table$gene
  scale_y_continuous(breaks = gene_pos_table[, "y_center"], 
                     labels = rep("", nrow(gene_pos_table)),
                     limits = species_axis_limits, 
                     expand = c(0, 0)) + 
  scale_fill_identity() +
  theme_minimal() +
  theme(axis.text.x = element_text(size = rel(1), hjust = 1, angle = 45), 
        # margin: top, right, bottom, and left
        plot.margin = unit(c(1, 0.2, 0.2, -0.7), "cm"), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

# Dendrogram plot species
plt_species_dendr <- ggplot(segment_data_species) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend, col = col)) + 
  scale_x_reverse() + 
  scale_y_continuous(breaks = gene_pos_table$y_center, 
                     labels = gene_pos_table$species, 
                     limits = species_axis_limits, 
                     expand = c(0, 0)) + 
  labs(x = "Distance", y = "", colour = "", size = "") +
  theme_minimal() + 
  theme(panel.grid.minor = element_blank(),
        legend.position = "None")

# Dendrogram plot traits
plt_traits_dendr <- ggplot(segment_data_traits) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend, col=col)) + 
  scale_x_continuous(position = "right") + 
  scale_y_continuous(breaks = sample_pos_table$x_center, 
                     labels = sample_pos_table$traits, 
                     limits = traits_axis_limits, 
                     position = "top", 
                     expand = c(0, 0)) + coord_flip() +
  theme_minimal() + 
  theme(panel.grid.minor = element_blank(), 
        axis.text.x = element_text(size = rel(0.8), angle = 45, hjust = 0),
        legend.position = "None")

library(cowplot)


final_plot <- plot_grid(
  NULL,    NULL,           NULL,             NULL, 
  NULL,    NULL,           plt_traits_dendr, NULL, 
  NULL,    plt_species_dendr, plt_hmap,         NULL, 
  NULL,    NULL,           NULL,             NULL, 
  nrow = 4, ncol = 4, align = "hv", 
  rel_heights = c(0.5, 1, 2, 0.5), 
  rel_widths = c(0.5, 1, 2, 0.5)
)









