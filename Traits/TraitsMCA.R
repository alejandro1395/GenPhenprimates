##########################################
#####DESCRIPTIVE PLOTS FROM JOE TRAITS####
##########################################
##########################################

library(phytools)
library(ggplot2)
library(ggnewscale)
library(ggtree)
packageVersion("phytools")
library(dplyr)
library(tidyr)
library(gridExtra)
library(RColorBrewer)

#READ ALL TRAITS FILE
ALLTraits <- read.csv("../Joe_traits/Primate_Traits/OUTPUT/MergedTraitsALL.txt",header=TRUE,sep="\t")
summary(ALLTraits)
glimpse(ALLTraits)

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
to_drop <- as.vector(ifelse(n <- sapply(subset_traits, function(x) length(levels(x))) == 1, "DROP", "NODROP"))
to_drop[to_drop == "DROP"]
ALLTraits$IUCN_ER
res.mca = MCA(subset_traits)
plot.MCA(res.mca)

# Summary of the variables
summary(res.mca)

#interpretation of the first component outcome
library("factoextra")
eig.val <- get_eigenvalue(res.mca)
fviz_screeplot(res.mca, addlabels = TRUE, ylim = c(0, 60)) #variance explained by each dimension
#first overlook of the components underlying both traits and species
fviz_mca_biplot(res.mca, 
                repel = FALSE, # Avoid text overlapping (slow if many point)
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
res.desc <- dimdesc(res.mca)
res.desc

#Factor loadings with FactorMineR

sweep(res.mca$var$coord,2,sqrt(res.mca$eig[1:ncol(res.mca$var$coord),1]),FUN="/")

#Hierarchical clustering of those variables
res.hcpc = HCPC(res.mca)
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


#2) STUDY OF INDIVIDUALS (IN OUR CASE, PRIMATE SPECIES)
