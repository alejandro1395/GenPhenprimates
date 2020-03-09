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
plot.MCA(res.mca)

# Summary of the variables
summary(res.mca)

#interpretation of the first component outcome
eig.val <- get_eigenvalue(res.mca)
fviz_screeplot(res.mca, addlabels = TRUE, ylim = c(0, 60)) #variance explained by each dimension

#first overlook of the components underlying both traits and species
fviz_mca_biplot(res.mca, 
                repel = FALSE, # Avoid text overlapping (slow if many point)
                ggtheme = theme_minimal(),
                col.ind = "chocolate4", col.var = "darkblue")

#TOP10 results only in contribution to axis
fviz_mca_biplot(res.mca, 
                repel = TRUE, # Avoid text overlapping (slow if many point),
                axes = c(1, 2),
                select.ind = list(contrib = 10), 
                select.var = list(contrib = 10),
                ggtheme = theme_minimal(),
                col.ind = "chocolate4", col.var = "darkblue")

#contribution of each variable trait to the dimensions
var <- get_mca_var(res.mca)
var$contrib
fviz_contrib(res.mca, choice = "var", axes = 1, top = 370)
fviz_contrib(res.mca, choice = "var", axes = 2, top = 230)

#now looking at the individuals level (species)
spcs <- get_mca_ind(res.mca)
fviz_mca_ind(res.mca, col.ind = "contrib",
             repel = TRUE,
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), # Avoid text overlapping (slow if many points)
             ggtheme = theme_minimal())

#SOME EXAMPLES
fviz_mca_ind(res.mca, 
             label = "none", # hide individual labels
             habillage = "Acantocheilonema.gracile", # color by groups 
             palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE, ellipse.type = "confidence",
             ggtheme = theme_minimal())
fviz_mca_ind(res.mca, 
             label = "none", # hide individual labels
             habillage = "Trypanosoma.brucei", # color by groups 
             palette = c("#00AFBB", "#E7B800"),
             addEllipses = TRUE, ellipse.type = "confidence",
             ggtheme = theme_minimal())
parasites <- colnames(ALLTraitsPresence[,380:735])
colourCount <- 2*length(parasites)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))
fviz_ellipses(res.mca, parasites,
              geom = "point",
              palette = getPalette(colourCount))
fviz_ellipses(res.mca, parasites,
             geom = "point")

#All dimensions analysis
res.desc <- dimdesc(res.mca, axes = c(1,2,3,4,5))
var <- get_mca_var(res.mca)
dimensions_cont <- as.data.frame(var$contrib)
rownames(dimensions_cont) <- rownames(var$contrib)

melted_dimensions_cont <- melt(setDT(dimensions_cont, keep.rownames = TRUE), "rn")

#TOP HITS FOR EACH DIMENSION
dim1 <- melted_dimensions_cont[which(melted_dimensions_cont$variable=="Dim 1"),] %>%
  group_by(variable) %>%
  top_n(n = 10, wt = as.numeric(as.character(value))) %>% slice(., 1:10)
dim2 <- melted_dimensions_cont[which(melted_dimensions_cont$variable=="Dim 2"),] %>%
  group_by(variable) %>%
  top_n(n = 10, wt = as.numeric(as.character(value))) %>% slice(., 1:10)
dim3 <- melted_dimensions_cont[which(melted_dimensions_cont$variable=="Dim 3"),] %>%
  group_by(variable) %>%
  top_n(n = 10, wt = as.numeric(as.character(value))) %>% slice(., 1:10)
dim4 <- melted_dimensions_cont[which(melted_dimensions_cont$variable=="Dim 4"),] %>%
  group_by(variable) %>%
  top_n(n = 10, wt = as.numeric(as.character(value))) %>% slice(., 1:10)
dim5 <- melted_dimensions_cont[which(melted_dimensions_cont$variable=="Dim 5"),] %>%
  group_by(variable) %>%
  top_n(n = 10, wt = as.numeric(as.character(value))) %>% slice(., 1:10)

#concat all dimensions
besthits <- rbind(dim1, dim2, dim3, dim4, dim5)

ggplot(besthits, aes(x = as.numeric(as.character(value)), y = as.numeric(as.character(value)))) + 
  geom_jitter(aes(colour = variable)) +
  geom_label_repel(aes(label = rn),
                   box.padding   = 0.5, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  facet_wrap(~variable, scales = "free") +
  theme_bw() + theme(legend.position="none")


var <- get_mca_var(res.mca)
var$contrib

#Factor loadings with FactorMineR

#sweep(res.mca$var$coord,2,sqrt(res.mca$eig[1:ncol(res.mca$var$coord),1]),FUN="/")

#Hierarchical clustering of those variables
res.hcpc <- HCPC(res.mca, graph = FALSE)
#descrip by categories
View(res.hcpc$desc.var$test.chi2)
View(res.hcpc$desc.var$category)
#descrip by components/dimensions
res.hcpc$desc.axes
res.hcpc$desc.ind
fviz_dend(res.hcpc, 
          cex = 1,                     # Label size
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco",           # Rectangle color
          labels_track_height = 0.1     # Augment the room for labels
)
fviz_cluster(res.hcpc,
             repel = TRUE,            # Avoid label overlapping
             show.clust.cent = TRUE, # Show cluster centers        # Color palette see ?ggpubr::ggpar
             ggtheme = theme_minimal(),
             main = "Factor map"
)




#LET'S STUDY PATTERNS OF CLUSTERED SPECIES HCPC values

sugg_spc_res.hcpc <- HCPC(res.mca, nb.clust = 10, graph = FALSE)
fviz_dend(sugg_spc_res.hcpc, 
          cex = 1,                     # Label size
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco",           # Rectangle color
          labels_track_height = 0.1     # Augment the room for labels
)

plot(sugg_spc_res.hcpc, choice = "map")



#convert to dataframe
#best tested method is melt
#TRAITS IN SPECIES CLUSTERS
sugg_spc_res.hcpc$desc.axes
plot(sugg_spc_res.hcpc)

traits_by_spc_clust_raw <- melt(head(sugg_spc_res.hcpc$desc.var$category))
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
View(clusters_spcs)





#FILL MATRIX VALUES AND PLOT (Here we reorder based on species cluster and combine information
#for each cluster coming from both species association to the cluster and trait association to
#the cluster
species_clusters_assigned <- unique(as.data.frame(cbind(sugg_spc_res.hcpc$data.clust$clust, rownames(sugg_spc_res.hcpc$data.clust))))
colnames(species_clusters_assigned) <- c("Cluster", "Species")


row_dend = sugg_spc_res.hcpc$call$t$tree
spc_data <- as.matrix(sugg_spc_res.hcpc$data.clust)
sorted_spc_data <- spc_data[match(labels(sugg_spc_res.hcpc$call$t$tree), rownames(spc_data)), ]
Species_names <- unique(rownames(sorted_spc_data))
Trait_names <- unique(colnames(sorted_spc_data))
plot_spc_clustmatrix <-  matrix(nrow = length(Species_names), 
                                ncol=length(Trait_names),
                                byrow = TRUE, 
                                dimnames = list(Species_names,Trait_names))

for(j in 1:ncol(plot_spc_clustmatrix)){
  for(i in 1:nrow(plot_spc_clustmatrix)){
    print(i)
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

#It takes a few minutes
row_dend = sugg_spc_res.hcpc$call$t$tree

#FULL MAP
plot_spc_clustmatrix <- plot_spc_clustmatrix[order(match(rownames(plot_spc_clustmatrix), row_dend$labels)), , drop = FALSE]
Heatmap(matrix = -plot_spc_clustmatrix, name = "Species clusters",
        row_names_gp = gpar(fontsize = 7),
        cluster_rows = color_branches(row_dend, k = 10),
        cluster_columns = FALSE,
        show_heatmap_legend = TRUE,
        row_dend_width = unit(100, "mm"),
        row_dend_side = "right",
        col = colorRamp2(c(0, 5, 10, 15, 20), brewer.pal(n=5, name="RdYlGn")))

#split clusters
#NOW WE CAN PRINT OUR DESIRED CLUSTER OF SPECIES

selected_species_clust1 <- as.character(species_clusters_assigned$Species[which(species_clusters_assigned$Cluster == 8)])
selected_matrix <- -plot_spc_clustmatrix[selected_species_clust1, ]
Heatmap(matrix = selected_matrix[, apply(selected_matrix,2,function(x) all(x>20))], name = "Species clusters",
        #cluster_rows = color_branches(row_dend, k = 1),
        cluster_columns = FALSE,
        show_heatmap_legend = FALSE,
        row_dend_width = unit(2, "mm"),
        row_dend_side = "right",
        col = colorRamp2(c(0, 5, 10, 15, 20), brewer.pal(n=5, name="RdYlGn")),
        column_names_gp = gpar(fontsize = 7),
        column_names_rot = 45,
        row_names_gp = gpar(fontsize = 5))




