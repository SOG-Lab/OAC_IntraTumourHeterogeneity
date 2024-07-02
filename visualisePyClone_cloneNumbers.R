## R script to visualise the PyClone results as circle plots ##
cat("*** Project ID: OAC_Heterogeneity_MR_MT ***\n")
cat("Loading: R script to visualise the PyClone results\n")

## Loading external libraries ##
library(ggplot2)
library(ggbeeswarm) ## 
library(ggpubr) ## stat_compare_means

## Loading external scripts ##
## Load script with general project settings and background functionality
source("SCRIPTS/setup.R")

## Environment ##
dir.out <- "RESULTS/PyCloneVI/Plots/"

mode <- "prePost"
data <- read.csv(paste0(dir.out, "CloneNumbers_", mode, ".csv"))
data$type <- sapply(strsplit(data$Sample.ID, "_"), "[", 3)

## Barplot of clone numbers per sample and per tumour
ggplot(data, aes(x = Sample.ID, y = Number.of.clones, fill = type)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(x = "", y = "Clone count") + 
  theme_test() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
ggsave(paste0(dir.out, "Barplot_cloneNumbers_", mode, ".png"), width = 15, height = 5)
ggsave(paste0(dir.out, "Barplot_cloneNumbers_", mode, ".pdf"), width = 15, height = 5)

data <- data[!grepl("tumour", data$Sample.ID),]

if(mode == "prePost"){
  data$Timepoint <- factor(data$Timepoint, levels = c("pre", "post"))
}

## Scatter plot number of clones vs cellularity
ggscatter(data, x = "Number.of.clones", y = "Cellularity",
          add = "reg.line",
          conf.int = TRUE,
          title = "Cellularity",
          xlab = "Number of Clones",
          ylab = "Cellularity [%]") + 
  stat_cor(method = "spearman")
ggsave(paste0(dir.out, "Scatter_ClonesVsCellularity_", mode, ".png"), width = 5, height = 4, type = "cairo")
ggsave(paste0(dir.out, "Scatter_ClonesVsCellularity_", mode, ".pdf"), width = 5, height = 4)

ggplot(data, aes(x = Number.of.clones, y = Cellularity)) +
  geom_point(aes(color = Timepoint)) + 
  theme_minimal() + 
  geom_smooth(method = 'lm', formula = y ~ x) +
  labs(x = "Clone count", y = "Cellularity [%]") 
ggsave(paste0(dir.out, "Dotplot_cloneNumber_Cellularity_", mode, ".png"), width = 6, height = 5)
ggsave(paste0(dir.out, "Dotplot_cloneNumber_Cellularity_", mode, ".pdf"), width = 6, height =5)

if(mode == "prePost"){
  ## Violin plot number of clones vs time point
  ggplot(data, aes(x = Timepoint, y = Number.of.clones)) +
    geom_violin() +
    stat_summary(fun = "mean", geom = "crossbar", color = "black") + 
    stat_compare_means(method = "t.test") + 
    geom_beeswarm(aes(color = Timepoint)) + 
    theme_minimal() + 
    labs(x = "Timepoint", y = "Clone count") 
  ggsave(paste0(dir.out, "Violinplot_cloneNumber_Timepoint_", mode, ".png"), width = 5, height = 5)
  ggsave(paste0(dir.out, "Violinplot_cloneNumber_Timepoint_", mode, ".pdf"), width = 5, height =5)
}

## Violin plot number of clones vs cellularity
data$Number.of.clones <- as.factor(data$Number.of.clones)
ggplot(data, aes(x = Number.of.clones, y = Cellularity)) +
  geom_violin() +
  stat_summary(fun = "mean", geom = "crossbar", color = "black") + 
  stat_compare_means(method = "kruskal.test") + 
  geom_beeswarm(aes(color = Timepoint)) + 
  theme_minimal() + 
  labs(x = "Clone count", y = "Cellularity [%]") 
ggsave(paste0(dir.out, "Violinplot_cloneNumber_Cellularity_", mode, ".png"), width = 6, height = 5)
ggsave(paste0(dir.out, "Violinplot_cloneNumber_Cellularity_", mode, ".pdf"), width = 6, height =5)

