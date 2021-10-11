#Load required packages


library(readxl)
library(vegan)
library(Hotelling)
library(microeco)
library(DataEditR)
library(tidyverse)
library(tidyquant)
library(xlsx) 
library(GUniFrac)
library(tidytree)
library(randomForest)
library(pheatmap)
library(agricolae)



#Microeco: https://chiliubio.github.io/microeco/

#Read data
asv <- read.xlsx2("C:/Users/Gloriana/Dropbox/Publicaciones/Publicaciones en progreso/5_Hongos y murcis/Análisis R/NumASVs_table.xlsx", sheetIndex = 1, row.names = 1)

sample <- read.xlsx2("C:/Users/Gloriana/Dropbox/Publicaciones/Publicaciones en progreso/5_Hongos y murcis/Análisis R/sample_table.xlsx", sheetIndex = 1, row.names = 1)

taxonomy <- read.xlsx2("C:/Users/Gloriana/Dropbox/Publicaciones/Publicaciones en progreso/5_Hongos y murcis/Análisis R/ASV_taxonomy.xlsx", sheetIndex = 1, row.names = 1)


library(dplyr)

asv_clean <- asv %>% mutate_if(is.character, as.numeric)

sapply(asv_clean, class)

# use pipe operator in magrittr package
library(magrittr)

# set.seed is used to fix the random number generation to make the results repeatable
set.seed(123)

# make the plotting background same with the tutorial
library(ggplot2)
theme_set(theme_bw())

# make the taxonomic information unified, important
taxonomy %<>% tidy_taxonomy

dataset <- microtable$new(sample_table = sample, otu_table = asv_clean, tax_table = taxonomy)
class(dataset)

print(dataset)

dataset$tidy_dataset()

print(dataset)

dataset$sample_sums() %>% range


dataset$cal_abund()
# return dataset$taxa_abund
class(dataset$taxa_abund)

# return dataset$taxa_abund
class(dataset$taxa_abund)

dir.create("taxa_abund")
dataset$save_abund(dirpath = "taxa_abund")


#calculate the alpha diversity
# If you want to add Faith's phylogenetic diversity, use PD = TRUE, this will be a little slow
dataset$cal_alphadiv(PD = FALSE)
# return dataset$alpha_diversity
class(dataset$alpha_diversity)


# save dataset$alpha_diversity to a directory
dir.create("alpha_diversity")
dataset$save_alphadiv(dirpath = "alpha_diversity")


#calculate the distance matrix of beta diversity using function cal_betadiv()
# If you do not want to calculate unifrac metrics, use unifrac = FALSE
# require GUniFrac package
dataset$cal_betadiv(unifrac = FALSE)
# return dataset$beta_diversity
class(dataset$beta_diversity)
# save dataset$beta_diversity to a directory
dir.create("beta_diversity")
dataset$save_betadiv(dirpath = "beta_diversity")


# create trans_abund object
# use 10 Phyla with the highest abundance in the dataset.
t1 <- trans_abund$new(dataset = dataset, taxrank = "Class", ntaxa = 10)
t2 <- trans_abund$new(dataset = dataset, taxrank = "Order", ntaxa = 10)
# tx object now include the transformed abundance data t1$abund_data and other elements for the following plotting

t1$plot_bar(others_color = "grey70", facet = "Species", xtext_keep = FALSE, legend_text_italic = FALSE)
ggsave("Figure1a_Abundance per class.png", width = 6, height = 4, dpi = 600)

t2$plot_bar(others_color = "grey70", facet = "Species", xtext_keep = FALSE, legend_text_italic = FALSE)
ggsave("Figure1b_Abundance per order.png", width = 6, height = 4, dpi = 600)
# return a ggplot2 object


# The groupmean parameter can be used to obtain the group-mean barplot.
t1 <- trans_abund$new(dataset = dataset, taxrank = "Class", ntaxa = 10, groupmean = "Species")
t2 <- trans_abund$new(dataset = dataset, taxrank = "Order", ntaxa = 10, groupmean = "Species")

t1$plot_bar(others_color = "grey70", legend_text_italic = FALSE)
ggsave("Figure1c_Relative abundance per class.png", width = 6, height = 4, dpi = 600)
t2$plot_bar(others_color = "grey70", legend_text_italic = FALSE)
ggsave("Figure1d_Relative abundance per order.png", width = 6, height = 4, dpi = 600)

# show 15 taxa at Class level
t1 <- trans_abund$new(dataset = dataset, taxrank = "Order", ntaxa = 15)
t1$plot_box(group = "Species")
ggsave("Figure2a_Differences in relative abundance per group at order level.png", width = 9, height = 6, dpi = 600)

t1 <- trans_abund$new(dataset = dataset, taxrank = "Order", ntaxa = 30)
t1$plot_heatmap(facet = "Species", xtext_keep = FALSE, withmargin = FALSE)
ggsave("Figure2b_Heatmap of differences in relative abundance per group at order level.png", width = 9, height = 6, dpi = 600)



#Venn analysis to analyze the unique and shared OTUs of species
# merge samples as one community for each group
dataset1 <- dataset$merge_samples(use_group = "Species")
# dataset1 is a new microtable object
# create trans_venn object
t1 <- trans_venn$new(dataset1, ratio = "seqratio")
t1$plot_venn()
ggsave("Figure3_Venn diagram.png", width = 9, height = 6, dpi = 600)
# The integer data is OTU number
# The percentage data is the sequence number/total sequence number



dataset1 <- dataset$merge_samples(use_group = "Species")
t1 <- trans_venn$new(dataset1)
# transform venn results to the sample-species table, here do not consider abundance, only use presence/absence information.
t2 <- t1$trans_venn_com(use_OTUs_frequency = TRUE)
# t2 is a new microtable class, each part is considered as a sample
class(t2)


#Alpha diversity 
t1 <- trans_alpha$new(dataset = dataset, group = "Species")
# return t1$alpha_stat
t1$alpha_stat[1:16, ]

t1$cal_diff(method = "t-test")
# return t1$res_alpha_diff
t1$res_alpha_diff[1:5, ]


t1$plot_alpha(measure = "Observed", group = "Species", pair_compare = TRUE)

ggsave("Figure4_Comparison of observed alpha diversity between communities.png", width = 6, height = 4, dpi = 600)



#The distance matrix of beta diversity can be transformed and plotted using trans_beta class. The analysis referred to the beta diversity in this class mainly include ordination, group distance, clustering and manova. 


# we first create an object and select NMDS for ordination
t1 <- trans_beta$new(dataset = dataset, group = "Species", measure = "jaccard", ordination = "PCoA")
# t1$res_ordination is the ordination result list
class(t1$res_ordination)

# calculate and plot sample distances within groups
t1$cal_group_distance()
# return t1$res_group_distance
t1$plot_group_distance(distance_pair_stat = TRUE)
ggsave("Figure5_Jaccard distance among communities.png", width = 6, height = 4, dpi = 600)


t1$plot_ordination(plot_color = "Species", plot_shape = "Species", plot_group_ellipse = TRUE)
ggsave("Figure6_Ordination for distance among communities.png", width = 9, height = 6, dpi = 600)



# manova for all groups
t1$cal_manova(cal_manova_all = TRUE)
t1$res_manova$aov.tab


#LEfSe combines the non-parametric test and linear discriminant analysis
t1 <- trans_diff$new(dataset = dataset, method = "lefse", group = "Species", alpha = 0.01, lefse_subgroup = NULL)
# t1$res_lefse is the LEfSe result
# t1$res_abund is the abundance information
t1$plot_lefse_bar(LDA_score = 4)
ggsave("Figure7_Differential abundance test.png", width = 8, height = 8, dpi = 600)

t1$plot_diff_abund(use_number = 1:10)

t1$res_lefse[1:10, ]

#The third approach is rf, which depends on the random forest[16, 17] and the non-parametric test. The current method can calculate random forest by bootstrapping like the method in LEfSe and only use the significant features. MeanDecreaseGini is selected as the indicator value in the analysis.

# use Genus level for parameter rf_taxa_level, if you want to use all taxa, change to "all"
# nresam = 1 and boots = 1 represent no bootstrapping and use all samples directly
t1 <- trans_diff$new(dataset = dataset, method = "rf", group = "Species", rf_taxa_level = "all")
# t1$res_rf is the result stored in the object
# plot the result
t2 <- t1$plot_diff_abund(use_number = 1:20, only_abund_plot = FALSE)
gridExtra::grid.arrange(t2$p1, t2$p2, ncol=2, nrow = 1, widths = c(2,2))
ggsave("Figure7_Random forest results for contribution of species to differentiating between communities.png", width = 9, height = 6, dpi = 600)
# the middle asterisk represent the significances




#Heat plots for function
#Read data
func_all <- read.xlsx2("C:/Users/Gloriana/Dropbox/Publicaciones/Publicaciones en progreso/5_Hongos y murcis/Análisis R/Function_table-unknown.xlsx", sheetIndex = 1, colClasses = c("character", "character", "character", "integer", rep("numeric", 2)))

func_all




numeric_cols <- sapply(func_all, Hmisc::all.is.numeric)

if (sum(numeric_cols) > 1)  {
  func_all[,numeric_cols] <- data.matrix(func_all[,numeric_cols])
} else {
  func_all[,numeric_cols] <- as.numeric(func_all[,numeric_cols])
}


func_all$Function <- reorder(func_all$Function,func_all$Abundance)


func.heatmap <- ggplot(data = func_all, mapping = aes(x = Sample,
                                                      y = Function,
                                                      fill = Abundance)) +
  geom_tile() +
  xlab(label = "Sample") +
  facet_grid(~ Site, switch = "x", scales = "free_x", space = "free_x") +
  scale_fill_gradient(name = "Abundance",
                      low = "white", high = "orange", limits=c(0,55))
theme_bw()

func.heatmap


#h_total <- func_all %>% 
#group_by(Sample) %>% 
#summarise(Frequency = sum(Frequency)) %>% 
#mutate(Function = 'Total')


v_total <- func_all %>% 
  group_by(Site, Function) %>% 
  summarise(Abundance = sum(Abundance)) %>% 
  mutate(Sample = 'Total')


p <- func.heatmap + 
  geom_point(data = v_total, 
             aes(color = Abundance), 
             size = 10, 
             shape = 19, show.legend = FALSE) +
  scale_color_gradient2(low = "white", 
                        high = "gray",
                        midpoint = 0) +
  geom_text(data = v_total, size = 4, aes(label = round(Abundance,2)))

p




ggsave("Figure8_Heat map with function_abundance.png", width = 12, height = 5, dpi = 600)





#Analysis by tree

#Read data
asv <- read.xlsx2("C:/Users/Gloriana/Dropbox/Publicaciones/Publicaciones en progreso/5_Hongos y murcis/Análisis R/NumASVs_table.xlsx", sheetIndex = 1, row.names = 1)

sample <- read.xlsx2("C:/Users/Gloriana/Dropbox/Publicaciones/Publicaciones en progreso/5_Hongos y murcis/Análisis R/sample_table_tree.xlsx", sheetIndex = 1, row.names = 1)

taxonomy <- read.xlsx2("C:/Users/Gloriana/Dropbox/Publicaciones/Publicaciones en progreso/5_Hongos y murcis/Análisis R/ASV_taxonomy.xlsx", sheetIndex = 1, row.names = 1)


library(dplyr)

asv_clean <- asv %>% mutate_if(is.character, as.numeric)

sapply(asv_clean, class)

# use pipe operator in magrittr package
library(magrittr)

# set.seed is used to fix the random number generation to make the results repeatable
set.seed(123)

# make the plotting background same with the tutorial
library(ggplot2)
theme_set(theme_bw())

# make the taxonomic information unified, important
taxonomy %<>% tidy_taxonomy

dataset <- microtable$new(sample_table = sample, otu_table = asv_clean, tax_table = taxonomy)
class(dataset)

print(dataset)

dataset$tidy_dataset()

print(dataset)

dataset$sample_sums() %>% range


dataset$cal_abund()
# return dataset$taxa_abund
class(dataset$taxa_abund)

# return dataset$taxa_abund
class(dataset$taxa_abund)

dir.create("taxa_abund")
dataset$save_abund(dirpath = "taxa_abund")


#calculate the alpha diversity
# If you want to add Faith's phylogenetic diversity, use PD = TRUE, this will be a little slow
dataset$cal_alphadiv(PD = FALSE)
# return dataset$alpha_diversity
class(dataset$alpha_diversity)


# save dataset$alpha_diversity to a directory
dir.create("alpha_diversity")
dataset$save_alphadiv(dirpath = "alpha_diversity")


# If you do not want to calculate unifrac metrics, use unifrac = FALSE
# require GUniFrac package
dataset$cal_betadiv(unifrac = FALSE)
# return dataset$beta_diversity
class(dataset$beta_diversity)
# save dataset$beta_diversity to a directory
dir.create("beta_diversity")
dataset$save_betadiv(dirpath = "beta_diversity")

# merge samples as one community for each group
dataset1 <- dataset$merge_samples(use_group = "Species")
# dataset1 is a new microtable object
# create trans_venn object
t1 <- trans_venn$new(dataset1, ratio = "seqratio")
t1$plot_venn()
ggsave("Figure9_Figure9_Venn diagram by tree.png", width = 6, height = 4, dpi = 600)


t1$plot_alpha(pair_compare = TRUE, measure = "Observed", shape = "Species")
ggsave("Figure10_Comparison of observed alpha diversity between communities by tree.png", width = 8, height = 4, dpi = 600)

# we first create an trans_beta object
# this operation invoke the distance matrix of bray in dataset$beta_diversity
t1 <- trans_beta$new(dataset = dataset, group = "Species", measure = "jaccard", ordination = "PCoA")

# t1$res_ordination is the ordination result list
class(t1$res_ordination)
# plot the PCoA result
t1$plot_ordination(plot_color = "Species", plot_shape = "Species", plot_group_ellipse = TRUE)
ggsave("Figure11_Ordination for distance among communities by tree.png", width = 6, height = 4, dpi = 600)

# calculate and plot sample distances within groups
t1$cal_group_distance()
# return t1$res_group_distance
t1$plot_group_distance(distance_pair_stat = TRUE)
ggsave("Figure12_Comparison of observed beta diversity between communities by tree.png", width = 8, height = 4, dpi = 600)

t1$cal_manova(cal_manova_all = TRUE)
t1$res_manova$aov.tab
t1$cal_manova(cal_manova_set = "Species")
t1$res_manova$aov.tab
