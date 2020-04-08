##############################################
#####VARIABILITY PLOTS FROM JOE TRAITS####
##############################################
##############################################

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
library(ggdendro)
library(viridis)

#READ LAST TRAITS FILE

all_traits <- read.csv("../data/Phenome/AllTraits25-03-20.tsv", header=T, sep = "\t")

#READ LAST METADATA FILE

metadata_traits <- read.csv("../data/Phenome/MetaData-AllTraits - MetaData-AllTraits.tsv", header=T, sep = "\t")

#LET'S GO FOR THE variance classification

############
#Binary data
############
binary_traits<-all_traits[,which(names(all_traits) %in% 
                                  metadata_traits$Trait[as.character(metadata_traits$DataType) == "binary"])]

binary_traits <- mutate_all(binary_traits, function(x) as.numeric(as.factor(x)))
View(sapply(binary_traits, unique))

sds_binary <- sapply(binary_traits, sd, na.rm = TRUE)
vars_binary <- as.data.frame(sapply(binary_traits, var, na.rm = TRUE))
vars_binary$traits <- rownames(vars_binary)
colnames(vars_binary) <- c("variance", "traits")

#filterwrong traits
traits_to_remove <- c("GenitalHair", "HeadFlesh", "Entamoeba.histolytica", "Monopolization")
vars_binary <- vars_binary[-which(vars_binary$traits %in% traits_to_remove), ]


ggplot(data=vars_binary[which(vars_binary$variance>=0.1), ]) + geom_bar(aes(x=reorder(traits, 
                                                  variance, decreasing=TRUE), 
                                        y=variance, fill = "lightblue"), stat="identity") +
  coord_flip() + scale_fill_manual(values = "lightblue") + theme_minimal() + theme(legend.title = element_blank(),
                                                                                   legend.position = "None")

ggplot(data=vars_binary[which(vars_binary$variance<0.026), ]) + geom_bar(aes(x=reorder(traits, 
                                                                                      -variance), 
                                                                            y=variance, fill = "lightblue"), stat="identity") +
  coord_flip() + scale_fill_manual(values = "darkblue") + theme_minimal() + theme(legend.title = element_blank(),
                                                                                   legend.position = "None")





##########
#continous data
##########

continous_traits<-all_traits[,which(names(all_traits) %in% 
                                   metadata_traits$Trait[as.character(metadata_traits$DataType) == "numeric"])]

continous_traits <- mutate_all(continous_traits, function(x) as.numeric(as.factor(x)))
#View(sapply(binary_traits, unique))

sds_cont <- sapply(continous_traits, sd, na.rm = TRUE)
vars_cont <- as.data.frame(sapply(continous_traits, var, na.rm = TRUE))
vars_cont$traits <- rownames(vars_cont)
colnames(vars_cont) <- c("variance", "traits")

#filterwrong traits
traits_to_remove <- c("Bacteria", "Fungus",	"Helminth",	"Protozoa")
vars_cont <- vars_cont[-which(vars_cont$traits %in% traits_to_remove), ]


ggplot(data=vars_cont[which(vars_cont$variance > 300), ]) + geom_bar(aes(x=reorder(traits, 
                                                                                      variance, decreasing=TRUE), 
                                                                            y=variance, fill = "red"), stat="identity") +
  coord_flip() + scale_fill_manual(values = "red") + theme_minimal() + theme(legend.title = element_blank(),
                                                                                   legend.position = "None")

ggplot(data=vars_cont[which(vars_cont$variance < 20), ]) + geom_bar(aes(x=reorder(traits, 
                                                                                       -variance), 
                                                                             y=variance, fill = "lightblue"), stat="identity") +
  coord_flip() + scale_fill_manual(values = "darkred") + theme_minimal() + theme(legend.title = element_blank(),
                                                                                  legend.position = "None")


#CATEGORICAL COUNT VARIABLES

categ_traits<-all_traits[,which(names(all_traits) %in% 
                                      metadata_traits$Trait[as.character(metadata_traits$DataType) == "factor"])]
#transform into counts

list_categ_vars <- sapply(categ_traits, table)


vars <- c()
name_vars <- names(list_categ_vars)
count <- 0
for (i in list_categ_vars){
  count <- count +1
  current_data <- as.numeric(as.vector(i))
  vars[count] <- sum( (current_data - as.numeric(mean(current_data) ))^2 )
}

TSS <- as.data.frame(cbind(name_vars, vars))
colnames(TSS) <- c("traits", "variance")

ggplot(data=TSS) + geom_bar(aes(x=reorder(traits, as.numeric(as.character(variance)), decreasing=TRUE), 
                                                                         y=as.numeric(as.character(variance)), 
                                fill = "red"), 
                            stat="identity") +
  coord_flip() + scale_fill_manual(values = "darkgreen") + theme_minimal() + theme(legend.title = element_blank(),
                                                                             legend.position = "None")





#BETTER TO DO THE COEF OF VARIANCE





