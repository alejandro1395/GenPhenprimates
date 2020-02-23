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

#READ ALL TRAITS FILE
ALLTraits <- read.csv("../Joe_traits/Primate_Traits/OUTPUT/MergedTraitsALL.txt",header=TRUE,sep="\t")
summary(ALLTraits)
glimpse(ALLTraits)


#WITHOUT ACCOUNTING FOR PHYLOGENY BY NOW

#a) Check continous variables amount
continuous <-select_if(ALLTraits, is.numeric)
summary(continuous)
ncol(continuous)

# standarize the continous columns because of not shared properties
ALLTraits_rescale <- ALLTraits %>%
  mutate_if(is.numeric, funs(as.numeric(scale(.))))
head(ALLTraits_rescale)

#b) Check factor variables 
factor <- data.frame(select_if(ALLTraits_rescale, is.factor))
ncol(factor)

#There is by now a total of 603 continous+binary variable traits and 21 categorical traits
#Let's plot them
# Create graph for each column
par(mfrow=c(1,1))
View(factor)
graph <- lapply(names(factor),
                function(x) 
                  ggplot(factor, aes(get(x))) +
                  geom_bar(fill="lightblue") +
                  theme(axis.text.x = element_text(angle = 90)) + xlab(colnames(x)))
graph[[18]]
grid.arrange(grobs = graph, ncol = 6)

#Correlation between all variables

library(GGally)
# Convert data to numeric
corr <- data.frame(lapply(recast_data, as.integer))
# Plot the graphggcorr(corr,
method = c("pairwise", "spearman"),
nbreaks = 6,
hjust = 0.8,
label = TRUE,
label_size = 3,
color = "grey50")

