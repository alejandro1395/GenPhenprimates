##############################################
#####SPREADNESS TREE PLOTS FROM JOE TRAITS####
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

all_traits <- read.csv("../data/Phenome/AllTraits28-03-20.tsv", header=T, sep = "\t")

ALLTraitsPresence <- all_traits %>%
  mutate_all(funs(ifelse(is.na(.), 0, 1)))
ALLTraitsPresence <- ALLTraitsPresence[,-c(1:2)]


#Now let's bind the Group info

groups_list <- strsplit(as.character(all_traits$GroupName), "_")
groups <- c()
for(i in 1:length(groups_list)){
  group <- paste(groups_list[[i]][1], "_", groups_list[[i]][2], sep="")
  groups[i] <- group}
groups

ALLTraitsPresenceGroupped <- as.data.frame(cbind(groups, ALLTraitsPresence))
ALLTraitsPresenceGroupped <- melt(ALLTraitsPresenceGroupped)
AllTraitsGroupSpreadness <- aggregate(.~ groups + variable, 
          data = ALLTraitsPresenceGroupped, FUN = sum)

totals <- count(groups)
AllTraitsGroupSpreadness$totals <- rep(totals$freq, length(AllTraitsGroupSpreadness$groups)/length(totals$freq))
head(AllTraitsGroupSpreadness)


#SUBSETTING Traits of interest (covering all families)

best_covered_traits <- list()
for (i in 1:length(unique(AllTraitsGroupSpreadness$variable))){
  data_subset <- subset(AllTraitsGroupSpreadness, 
         AllTraitsGroupSpreadness$variable == unique(AllTraitsGroupSpreadness$variable)[i])
  if (0 %in% as.numeric(as.character(data_subset$value))){
    next()
  }
  else{
    best_covered_traits[[i]] <- data_subset
  }
}

best_cases <- do.call(rbind, best_covered_traits)


#PLOT IT
library(randomcoloR)
n <- length(unique(groups))
palette <- distinctColorPalette(n)

ggplot(data=best_cases) + geom_bar(aes(x=variable, y=value, fill=groups), stat="identity") +
  coord_flip() + scale_fill_manual(values = palette) + theme_minimal() + theme(legend.title = element_blank())

#RANK IT

ggplot(data=best_cases) + geom_bar(aes(x=reorder(variable, value, decreasing=TRUE), y=value, fill=groups), stat="identity") +
  coord_flip() + scale_fill_manual(values = palette) + theme_minimal() + theme(legend.title = element_blank())


#SUBSETTING Traits of interest (covering N-1 families)


second_best_covered_traits <- list()
for (i in 1:length(unique(AllTraitsGroupSpreadness$variable))){
  data_subset <- subset(AllTraitsGroupSpreadness, 
                        AllTraitsGroupSpreadness$variable == unique(AllTraitsGroupSpreadness$variable)[i])
  missing_families <- sum(as.numeric(as.character(data_subset$value)) == 0)
  if (5 > missing_families & missing_families >= 1){
    second_best_covered_traits[[i]] <- data_subset
  }
  else{
    next()
  }
}

second_best_cases <- do.call(rbind, second_best_covered_traits)
ggplot(data=second_best_cases) + geom_bar(aes(x=reorder(variable, value, decreasing=TRUE), y=value, fill=groups), stat="identity") +
  coord_flip() + scale_fill_manual(values = palette) + theme_minimal() + theme(legend.title = element_blank())


#SUBSETTING TRAITS OF interest (covering N-2)

third_best_covered_traits <- list()
for (i in 1:length(unique(AllTraitsGroupSpreadness$variable))){
  data_subset <- subset(AllTraitsGroupSpreadness, 
                        AllTraitsGroupSpreadness$variable == unique(AllTraitsGroupSpreadness$variable)[i])
  missing_families <- sum(as.numeric(as.character(data_subset$value)) == 0)
  if (5 > missing_families & missing_families >= 0){
    third_best_covered_traits[[i]] <- data_subset
  }
  else{
    next()
  }}

  third_best_cases <- do.call(rbind, third_best_covered_traits)
  ggplot(data=third_best_cases) + geom_bar(aes(x=reorder(variable, value, decreasing=TRUE), y=value, fill=groups), stat="identity") +
    coord_flip() + scale_fill_manual(values = palette) + theme_minimal() + theme(legend.title = element_blank())
  
  
  
  
  
  
  
  
  
##############################
#NOW LET'S MERGE THE VARIANCE INFORMATION
##############################
#READ LAST METADATA FILE
  
metadata_traits <- read.csv("../data/Phenome/MetaData-AllTraits - MetaData-AllTraits.tsv", header=T, sep = "\t")

############
#Binary data
############

binary_traits<-all_traits[,which(names(all_traits) %in% 
                                     metadata_traits$Trait[as.character(metadata_traits$DataType) == "binary"])]
  

binary_traits$StableGroups<- as.character(binary_traits$StableGroups)
binary_traits$StableGroups[as.character(binary_traits$StableGroups) == "Yes"] <- 1
binary_traits$StableGroups[as.character(binary_traits$StableGroups) == "No"] <- 0
binary_traits$Monopolization<- as.character(binary_traits$Monopolization)
binary_traits$Monopolization[as.character(binary_traits$Monopolization) == "Yes"] <- 1
binary_traits$Monopolization[as.character(binary_traits$Monopolization) == "No"] <- 0
binary_traits <- mutate_all(binary_traits, function(x) as.numeric(as.character(x)))
BinaryTraitsGroupped <- as.data.frame(cbind(groups, binary_traits))

#split by families
families_binary_list <- split(BinaryTraitsGroupped[,-1], BinaryTraitsGroupped$groups)
final_binary_list <- list()
name_binary_vars <- names(families_binary_list)
count_spc <- 0
for (spc in families_binary_list){
  count_spc <- count_spc + 1
  new_list <- sapply(families_binary_list[[count_spc]], table)
  names_new_list <- names(new_list)
  traits_vec <- c()
  count_tr <- 0
  for (trait in new_list){
    if (length(trait)!=0){
    substract <- sum((as.numeric(as.character(trait))/sum(as.numeric(as.character(trait))))^2)
    disp <- 1-substract
    count_tr<- count_tr + 1
    traits_vec[count_tr] <- disp
    }
    else{
      count_tr<- count_tr + 1
      traits_vec[count_tr] <- NA
    }}
  subset_spc <- as.data.frame(cbind(rep(names(families_binary_list)[count_spc], 
                                        length(names_new_list)), 
                                    names_new_list, traits_vec))
  print(traits_vec)
  final_binary_list[[count_spc]] <- subset_spc
}

disparities_binary <- do.call(rbind, final_binary_list)
colnames(disparities_binary) <- c("Family", "traits", "disparity")

between_disparity <- c()
for (i in 1:length(unique(disparities_binary$traits))){
  data_subset <- subset(disparities_binary, 
                        disparities_binary$traits == unique(disparities_binary$traits)[i])
  mean_of_disparity <- mean(as.numeric(as.character(data_subset$disparity)), na.rm=TRUE)
  between_disparity[i] <- mean_of_disparity}

disparities_binary$between_disparities <- rep(between_disparity, length(unique(disparities_binary$Family)))
disparities_binary <- as.data.frame(disparities_binary)



#######################
#PLOT IN TREE TOGETHER#
#######################

#extrasum for visualization
#View(disparities_binary)
ggplot(data=disparities_binary[which(disparities_binary$traits %in% third_best_cases$variable),]) + 
  geom_bar(aes(x=Family, 
               y=as.numeric(as.character(disparity)), 
               fill=Family), 
           stat="identity",
           colour="black") +
  scale_fill_manual(values = palette) + theme_minimal() + theme(legend.title = element_blank()) +
  facet_wrap(~traits) +
  theme(axis.text.x = element_blank()) +
  geom_hline(aes(yintercept = as.numeric(as.character(between_disparities))), color="black", linetype="dashed")





###################
##################
#CONTINOUS data####
##################

continous_traits<-all_traits[,which(names(all_traits) %in% 
                                   metadata_traits$Trait[as.character(metadata_traits$DataType) == "numeric"])]


continous_traits <- mutate_all(continous_traits, function(x) as.numeric(as.character(x)))
#View(sapply(binary_traits, unique))

#filterwrong traits
traits_to_remove <- c("Bacteria", "Fungus",	"Helminth",	"Protozoa")
ContinousTraitsGroupped <- as.data.frame(cbind(groups, continous_traits))
ContinousTraitsGroupped <- ContinousTraitsGroupped[,-which(colnames(ContinousTraitsGroupped) %in% traits_to_remove)]


ContinousTraitsGroupped <- melt(ContinousTraitsGroupped)
ContinousTraitsMeans <- ContinousTraitsGroupped %>%
  group_by(groups, variable) %>%
  summarize_all(mean, na.rm = TRUE)

ContinousTraitsSds <- ContinousTraitsGroupped %>%
  group_by(groups, variable) %>%
  summarize_all(sd, na.rm = TRUE)

AllStatsContinous <- merge(ContinousTraitsMeans, ContinousTraitsSds, by=c("groups", "variable"))
colnames(AllStatsContinous) <- c("Family", "Trait", "Mean", "Sd")
AllStatsContinous$coef_var <- AllStatsContinous$Sd/AllStatsContinous$Mean*100

#Vector with the same but for each trait in whole phylogeny (with MEANS)

between_coef_var <- c()
for (i in 1:length(unique(AllStatsContinous$Trait))){
  data_subset <- subset(AllStatsContinous, 
                        AllStatsContinous$Trait == unique(AllStatsContinous$Trait)[i])
  mean_of_means <- mean(data_subset$Mean, na.rm=TRUE)
  sd_of_means <- sd(data_subset$Mean, na.rm=TRUE)
  between_coef_var[i] <- sd_of_means/mean_of_means*100}

AllStatsContinous$between_coef_vars <- rep(between_coef_var, length(unique(AllStatsContinous$Family)))
ggplot(data=AllStatsContinous[which(AllStatsContinous$Trait %in% third_best_cases$variable),]) + 
  geom_bar(aes(x=Family, 
               y=as.numeric(as.character(coef_var)), 
               fill=Family), 
           stat="identity",
           colour="black") +
  scale_fill_manual(values = palette) + theme_minimal() + theme(legend.title = element_blank()) +
  facet_wrap(~Trait) +
  theme(axis.text.x = element_blank()) +
  geom_hline(aes(yintercept = as.numeric(as.character(between_coef_vars))), color="black", linetype="dashed")





############
#Categorical data
############


categorical_traits<-all_traits[,which(names(all_traits) %in% 
                                   metadata_traits$Trait[as.character(metadata_traits$DataType) == "factor"])]

CategoricalTraitsGroupped <- as.data.frame(cbind(groups, categorical_traits))


#split by families
families_cat_list <- split(CategoricalTraitsGroupped[,-1], CategoricalTraitsGroupped$groups)
final_cat_list <- list()
name_cat_vars <- names(families_cat_list)
count_spc <- 0
for (spc in families_cat_list){
  count_spc <- count_spc + 1
  new_list <- sapply(families_cat_list[[count_spc]], table)
  names_new_list <- names(new_list)
  traits_vec <- c()
  count_tr <- 0
  for (trait in new_list){
    if (length(trait)!=0){
      substract <- sum((as.numeric(as.character(trait))/sum(as.numeric(as.character(trait))))^2)
      disp <- 1-substract
      count_tr<- count_tr + 1
      traits_vec[count_tr] <- disp
    }
    else{
      count_tr<- count_tr + 1
      traits_vec[count_tr] <- NA
    }}
  subset_spc <- as.data.frame(cbind(rep(names(families_cat_list)[count_spc], 
                                        length(names_new_list)), 
                                    names_new_list, traits_vec))
  print(traits_vec)
  final_cat_list[[count_spc]] <- subset_spc
}

disparities_cat <- do.call(rbind, final_cat_list)
colnames(disparities_cat) <- c("Family", "traits", "disparity")

between_disparity <- c()
for (i in 1:length(unique(disparities_cat$traits))){
  data_subset <- subset(disparities_cat, 
                        disparities_cat$traits == unique(disparities_cat$traits)[i])
  mean_of_disparity <- mean(as.numeric(as.character(data_subset$disparity)), na.rm=TRUE)
  between_disparity[i] <- mean_of_disparity}

disparities_cat$between_disparities <- rep(between_disparity, length(unique(disparities_cat$Family)))
disparities_cat <- as.data.frame(disparities_cat)


#######################
#PLOT IN TREE TOGETHER#
#######################

#extrasum for visualization
#View(disparities_binary)
ggplot(data=disparities_cat[which(disparities_cat$traits %in% third_best_cases$variable),]) + 
  geom_bar(aes(x=Family, 
               y=as.numeric(as.character(disparity)), 
               fill=Family), 
           stat="identity",
           colour="black") +
  scale_fill_manual(values = palette) + theme_minimal() + theme(legend.title = element_blank()) +
  facet_wrap(~traits) +
  theme(axis.text.x = element_blank()) +
  geom_hline(aes(yintercept = as.numeric(as.character(between_disparities))), color="black", linetype="dashed")


















#TRASH SCRIPT BY NOW
#NEXT LINES


BinaryTraitsGroupped <- as.data.frame(cbind(groups, binary_traits))
BinaryTraitsGroupped <- melt(BinaryTraitsGroupped)
BinaryTraitsMeans <- BinaryTraitsGroupped %>%
  group_by(groups, variable) %>%
  summarize_all(mean, na.rm = TRUE)

View(BinaryTraitsSds) <- BinaryTraitsGroupped %>%
  group_by(groups, variable) %>%
  summarize_all(sd, na.rm = TRUE)

AllStatsBinary <- merge(BinaryTraitsMeans, BinaryTraitsSds, by=c("groups", "variable"))
colnames(AllStatsBinary) <- c("Family", "Trait", "Mean", "Sd")
AllStatsBinary$coef_var <- AllStatsBinary$Sd/AllStatsBinary$Mean*100

#Vector with the same but for each trait in whole phylogeny (with MEANS)

between_coef_var <- c()
for (i in 1:length(unique(AllStatsBinary$Trait))){
  data_subset <- subset(AllStatsBinary, 
                        AllStatsBinary$Trait == unique(AllStatsBinary$Trait)[i])
  mean_of_means <- mean( data_subset$Mean, na.rm=TRUE)
  sd_of_means <- sd( data_subset$Mean, na.rm=TRUE)
  print(sd_of_means)
  between_coef_var[i] <- sd_of_means/mean_of_means*100}

AllStatsBinary$between_coef_vars <- rep(between_coef_var, length(unique(AllStatsBinary$Family)))

ggplot(data=AllStatsBinary[which(AllStatsBinary$Trait %in% third_best_cases$variable),]) + 
  geom_bar(aes(x=reorder(Family, coef_var, decreasing=TRUE), y=coef_var, fill=Family), 
                                       stat="identity") +
 scale_fill_manual(values = palette) + theme_minimal() + theme(legend.title = element_blank()) +
  facet_wrap(~Trait, scales="free_y") +
  geom_hline(aes(yintercept = between_coef_vars), col="dodgerblue4")









#View(sapply(binary_traits, unique))
  
sds_binary <- sapply(binary_traits, sd, na.rm = TRUE)
means_binary <- sapply(binary_traits, mean, na.rm = TRUE)
coef_var <- sds_binary/means_binary*100

vars_binary <- as.data.frame(sapply(binary_traits, var, na.rm = TRUE))
vars_binary$traits <- rownames(vars_binary)
colnames(vars_binary) <- c("variance", "traits")
  
  #filterwrong traits
  traits_to_remove <- c("GenitalHair", "HeadFlesh", "Entamoeba.histolytica", "Monopolization")
  vars_binary <- vars_binary[-which(vars_binary$traits %in% traits_to_remove), ]
  
  