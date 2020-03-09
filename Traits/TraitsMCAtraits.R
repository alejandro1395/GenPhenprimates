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
BiocManager::install("ComplexHeatmap")


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












#######################################################
#######################################################
#######################################################
#2) STUDY OF INDIVIDUALS (IN OUR CASE, PRIMATE TRAITS)
#######################################################
#######################################################



#before calling MCA algorithm for factor analysis of traitsxspecies,
#let's transform the data into presence/absence of traits

ALLTraitsPresence <- ALLTraits %>%
  mutate_all(funs(ifelse(is.na(.), 0, 1)))

#Load FactorMineR package
library(FactoMineR)


#########################
#########################
#1) STUDY OF CATEGORIES
#First let's do the first Component analysis of our matrix data to know
#the variables closer together in the axis

subset_traits <- ALLTraitsPresence[,-which(names(ALLTraitsPresence) %in% c("SpeciesBROAD", "GroupName"))]
i=0
while(i < ncol(subset_traits)){
  i=i+1  
  subset_traits[,i] = as.factor(subset_traits[,i])
}
row.names(subset_traits) <- as.character(ALLTraits$SpeciesBROAD)
summary(subset_traits)
transposed_traits <- as.data.frame(t(subset_traits))
dim(transposed_traits)
#remove variables with only 1 factor level (both for columns and rows)
cond1 <- sapply(transposed_traits, function(col) sum(nlevels(col)) > 1)
transposed_traits <- subset(transposed_traits, select=cond1)
dim(transposed_traits)

##### use this for res.mca
summary(transposed_traits)
inv.res.mca = MCA(transposed_traits)

# Summary of the variables
summary(inv.res.mca)

#interpretation of the first component outcome
library("factoextra")
eig.val <- get_eigenvalue(inv.res.mca)
fviz_screeplot(inv.res.mca, addlabels = TRUE, ylim = c(0, 90)) #variance explained by each dimension
#first overlook of the components underlying both traits and species
fviz_mca_biplot(inv.res.mca, 
                repel = FALSE, # Avoid text overlapping (slow if many point)
                label = "BodyMass_kg",
                ggtheme = theme_minimal(),
                col.ind = "chocolate4", col.var = "darkblue")
#TOP10
fviz_mca_biplot(inv.res.mca, 
                repel = TRUE, # Avoid text overlapping (slow if many point),
                axes = c(1, 2),
                select.ind = list(contrib = 10), 
                select.var = list(contrib = 10),
                ggtheme = theme_minimal(),
                col.ind = "chocolate4", col.var = "darkblue")



#contribution of each variable trait to the dimensions
var <- get_mca_var(inv.res.mca)
var$contrib
fviz_contrib(inv.res.mca, choice = "var", axes = 1, top = 370)
fviz_contrib(inv.res.mca, choice = "var", axes = 2, top = 230)

#now looking at the individuals level (species)
traits <- get_mca_ind(inv.res.mca)
fviz_mca_ind(inv.res.mca, col.ind = "contrib",
             repel = TRUE,
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), # Avoid text overlapping (slow if many points)
             ggtheme = theme_minimal())

#All dimensions analysis
inv.res.desc <- dimdesc(inv.res.mca)
inv.var <- get_mca_var(inv.res.mca)
inv.dimensions_cont <- as.data.frame(inv.var$contrib)
rownames(inv.dimensions_cont) <- as.character(rownames(inv.var$contrib))

melted_dimensions_cont <- melt(setDT(inv.dimensions_cont, keep.rownames = TRUE), "rn")

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
besthits <- rbind(dim1, dim2, dim3)

ggplot(besthits, aes(x = as.numeric(as.character(value)), y = as.numeric(as.character(value)))) + 
  geom_jitter(aes(colour = variable)) +
  geom_label_repel(aes(label = rn),
                   box.padding   = 0.5, 
                   point.padding = 0.5,
                   force = 3,
                   segment.color = 'grey50') +
  facet_wrap(~variable) +
  theme_bw() + theme(legend.position="none")


#contribution of each variable trait to the dimensions
var <- get_mca_var(res.mca)
var$contrib
fviz_contrib(res.mca, choice = "var", axes = 1, top = 900)
fviz_contrib(res.mca, choice = "var", axes = 2, top = 50)

#now looking at the SPECIES level (TRAITS)
spcs <- get_mca_ind(inv.res.mca)
fviz_mca_ind(inv.res.mca, col.ind = "contrib",
             repel = TRUE,
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), # Avoid text overlapping (slow if many points)
             ggtheme = theme_minimal())

#SOME EXAMPLES
fviz_mca_ind(res.mca, 
             label = "none", # hide individual labels
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
res.desc[[1]][1]

#Factor loadings with FactorMineR

sweep(res.mca$var$coord,2,sqrt(res.mca$eig[1:ncol(res.mca$var$coord),1]),FUN="/")

#Hierarchical clustering of those variables
res.hcpc = HCPC(inv.res.mca$var$contrib, graph = FALSE)
#descrip by categories
View(res.hcpc$desc.var$test.chi2)
View(res.hcpc$desc.var$category)
#descrip by components/dimensions
res.hcpc$desc.axes
res.hcpc$desc.ind
fviz_dend(res.hcpc, 
          cex = 2,                     # Label size
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco",           # Rectangle color
          labels_track_height = 1,
          horiz = TRUE,
          repel=TRUE,
          phylo_layout = layout_as_tree# Augment the room for labels
)

fviz_cluster(res.hcpc,
             repel = TRUE,            # Avoid label overlapping
             show.clust.cent = TRUE, # Show cluster centers        # Color palette see ?ggpubr::ggpar
             ggtheme = theme_minimal(),
             main = "Factor map"
)




#LET'S STUDY PATTERNS OF CLUSTERED TRAITS HCPC values

sugg_spc_res.hcpc <- HCPC(inv.res.mca, nb.clust = 10, graph = FALSE)
fviz_dend(sugg_spc_res.hcpc, 
          cex = 1,                     # Label size
          palette = "jco",               # Color palette see ?ggpubr::ggpar
          rect = TRUE, rect_fill = TRUE, # Add rectangle around groups
          rect_border = "jco",           # Rectangle color
          labels_track_height = 0.1     # Augment the room for labels
)

#plot(sugg_spc_res.hcpc, choice = "map")



#convert to dataframe
#best tested method is melt
#SPECIES IN TRAITS CLUSTERS
sugg_spc_res.hcpc$desc.axes
plot(sugg_spc_res.hcpc)

species_by_trait_clust_raw <- melt(sugg_spc_res.hcpc$desc.var$category)
species_by_trait_clust_raw <- species_by_trait_clust_raw[which(species_by_trait_clust_raw$Var2 == "p.value") ,]
variables <- unlist(strsplit(as.vector(species_by_trait_clust_raw$Var1), "="))
Species <- variables[seq(1, length(variables), by=2)]
Categories <- variables[seq(2, length(variables), by=2)]
species_by_trait_clust$Species <- as.data.frame(cbind(species_by_trait_clust_raw$L1, 
                                           as.character(Species),
                                           species_by_trait_clust_raw$value))
colnames(species_by_trait_clust) <- c("Cluster", "Species", "p-value")
species_by_trait_clust <- species_by_trait_clust[duplicated(species_by_trait_clust), ]
#View(species_by_trait_clust)_raw)





#TRAITS IN SPECIES CLUSTERS

#The easiest idea here is to build a heatmap based on group belonging
#in this case, cluster belonging (and then compare with a clop)
traits_by_trait_clust_raw <- as.data.frame(as.matrix(unlist(sugg_spc_res.hcpc$desc.ind$para)))
clusters_traits_raw_list <- strsplit(row.names(traits_by_trait_clust_raw), split='.', 
                                               fixed=TRUE)
clusters <- unlist(lapply(clusters_traits_raw_list, `[[`, 1))#get cluster numbers
traits <- gsub("^[0-9]+.","",row.names(traits_by_trait_clust_raw))

clusters_traits <- as.data.frame(cbind(clusters, as.character(traits), traits_by_trait_clust_raw$V1))
colnames(clusters_traits) = c("Cluster", "Trait", "coord")
View(clusters_traits)





#FILL MATRIX VALUES AND PLOT (Here we reorder based on species cluster and combine information
#for each cluster coming from both species association to the cluster and trait association to
#the cluster
traits_clusters_assigned <- unique(as.data.frame(cbind(sugg_spc_res.hcpc$data.clust$clust, rownames(sugg_spc_res.hcpc$data.clust))))
colnames(traits_clusters_assigned) <- c("Cluster", "Trait")


row_dend <- plot(sugg_spc_res.hcpc$call$t$tree)
trait_data <- as.matrix(sugg_spc_res.hcpc$data.clust)
sorted_trait_data <- trait_data[match(labels(sugg_spc_res.hcpc$call$t$tree), rownames(trait_data)), ]
Trait_names <- unique(rownames(sorted_trait_data))
Species_names <- unique(colnames(sorted_trait_data))
plot_trait_clustmatrix <-  matrix(nrow = length(Species_names), 
                                ncol=length(Trait_names),
                                byrow = TRUE, 
                                dimnames = list(Species_names,Trait_names))

for(j in 1:ncol(plot_trait_clustmatrix)){
  for(i in 1:nrow(plot_trait_clustmatrix)){
    #look for the significance association values in individual dfs
    current_clast <- as.numeric(as.character(traits_clusters_assigned$Cluster[traits_clusters_assigned$Trait==Trait_names[j]]))
    if (Trait_names[j] %in% clusters_traits$Trait) {
      trait_effect <- as.numeric(as.character(clusters_traits$coord[as.character(clusters_traits$Trait)==Trait_names[j]]))
    }
    else{
      trait_effect = 1
    }
    if (Species_names[i] %in% species_by_trait_clust$Species[as.numeric(as.character(species_by_trait_clust$Cluster))==current_clast]) {
      current_clast <- as.numeric(as.character(traits_clusters_assigned$Cluster[traits_clusters_assigned$Trait==Trait_names[j]]))
      species_effect <- as.numeric(as.character(species_by_trait_clust$`p-value`[species_by_trait_clust$Species==Species_names[i] &
                                                                              as.numeric(as.character(species_by_trait_clust$Cluster))==current_clast]))
    }
    else{
      species_effect = 1
    }
    plot_trait_clustmatrix[i, j] = log(species_effect) + log(trait_effect)
  }}

#It takes a few minutes
col_dend = sugg_spc_res.hcpc$call$t$tree

#FULL MAP
plot_trait_clustmatrix <- plot_trait_clustmatrix[, order(match(colnames(plot_trait_clustmatrix), col_dend$labels))]
Heatmap(matrix = -plot_trait_clustmatrix, name = "Traits clusters",
        row_names_gp = gpar(fontsize = 7),
        cluster_columns = color_branches(col_dend, k = 10),
        cluster_rows = FALSE,
        show_heatmap_legend = FALSE,
        column_dend_height = unit(100, "mm"),
        column_dend_side = "bottom",
        col = colorRamp2(c(0, 5, 10, 15, 20), brewer.pal(n=5, name="RdYlGn")))

#split clusters
#NOW WE CAN PRINT OUR DESIRED CLUSTER OF SPECIES

selected_traits_clust1 <- as.character(traits_clusters_assigned$Trait[which(traits_clusters_assigned$Cluster == 9)])
selected_matrix <- -plot_trait_clustmatrix[, selected_traits_clust1]
Heatmap(matrix = selected_matrix[apply(selected_matrix,1,function(x) all(x>10)),], name = "Traits clusters",
        #cluster_rows = color_branches(row_dend, k = 1),
        cluster_columns = FALSE,
        show_heatmap_legend = FALSE,
        column_dend_height = unit(100, "mm"),
        column_dend_side = "bottom",
        col = colorRamp2(c(0, 5, 10, 15, 20), brewer.pal(n=5, name="RdYlGn")),
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 1))




