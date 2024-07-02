cat("Project ID: OAC_Heterogeneity_MR_MT\n")
## R script to map neoantigens to clones from PyClone IV ##
cat("Loading: R script to map neoantigens to clones from PyClone IV \n")

## Loading external libraries ##

## Loading external scripts ##
## Load script with general project settings and background functionality
source("SCRIPTS/setup.R")

## Environment ##
dir.input <- "RESULTS/PyCloneVI/PyCloneResults/"
dir.out <- "RESULTS/PyCloneVI/Neoantigens/"
.checkAndCreateDir(dir.out)
## PyClone results files
files <- list.files(dir.input, pattern = "_filtered[.]csv$", full.names = TRUE)
## Load annotation data
cohort <- general.load.data(paste0(dir.data, "Annotation/Cohort.csv")) ## Contains sample IDs and clinical information for each patient
## Patients with pre and post treatment samples
patients <- cohort$Patient[cohort$Post_Samples > 0]

neoantigens <- general.load.data("DATA/Neoantigens/All_MR_Neoantigens.csv", silent = TRUE)
neoantigens$X <- NULL
## Remove neoantigens that don't have a WT score (frameshift mutation neoantigens)
neoantigens <- neoantigens[!is.na(neoantigens$Corresponding.WT.Score),]
## Remove neoantigens on chromosomes X/Y/MT
neoantigens <- neoantigens[neoantigens$Chromosome %in% paste0("chr", 1:22),]
neoantigens$DAI <- neoantigens$Corresponding.WT.Score - neoantigens$Best.MT.Score
## Add mutation ID that matches the one in the pileup file
neoantigens$mutation_id <- paste(neoantigens$Gene.Name, neoantigens$Chromosome, neoantigens$Stop, neoantigens$Stop, neoantigens$Reference, neoantigens$Variant, sep = "_")

neoantigensPerClone <- function(){
  ## Assign neoantigens to pre/post clonal structures for the 10 patients that have pre/post files available
  files.prePost <- files[grepl(paste(patients, collapse = "|"), files) & !grepl("preOnly", files) & !grepl("old", files)]
  i <- 1
  total <- length(files.prePost)
  a <- sapply(files.prePost, function(file){
    patient <- strsplit(basename(file), "_")[[1]][1]
    cat(i, "/", total, "| Assigning neoantigens to clonal structure from", patient, "\n")
    i <<- i + 1
    ## Load and prepare PyClone result
    clones <- .preparePyCloneData(file)
    
    ## Merge neoantigens with cluster IDs
    neos <- merge(neoantigens, clones[, c("Patient", "mutation_id", "cluster_id", "cluster")], by.x = c("Patient_ID", "mutation_id"), by.y = c("Patient", "mutation_id"), all.x = TRUE)
    if("cluster_id.x" %in% names(neos)){
      neos$cluster_id.x[is.na(neos$cluster_id.x)] <- neos$cluster_id.y[is.na(neos$cluster_id.x)]
      neos$cluster.x[is.na(neos$cluster.x)] <- neos$cluster.y[is.na(neos$cluster.x)]
      neos$cluster_id.y <- neos$cluster.y <- NULL
      names(neos) <- gsub("[.]x$", "", names(neos))
    }
    neoantigens <<- neos
  })
  ## pre/post neoantigens
  neos.prepost <- neoantigens[neoantigens$Patient_ID %in% patients,]
  general.save.data(neos.prepost, paste0(dir.out, "Neoantigens_clusters_PrePost"))
  
  ## For plotting remove duplicated entries - keep per patient not per sample
  neos.prepost <- unique(neos.prepost[, names(neos.prepost) %in% c("Patient_ID", "mutation_id", "Best.MT.Score", "Corresponding.WT.Score", "DAI", "cluster_id", "cluster")])
  
  ## Plotting
  #neos.prepost <- .plot.DAI.clusters(neos.prepost, paste0(dir.out, "prePost_DAI"))
  
  ## Assign neoantigens to pre clonal structures for all 29 patients 
  neoantigens$cluster <- neoantigens$cluster_id <- NULL
  files.pre <- files[(!grepl(paste(patients, collapse = "|"), files) | grepl("preOnly", files))]
  i <- 1
  total <- length(files.pre)
  
  a <- sapply(files.pre, function(file){
    patient <- strsplit(basename(file), "_")[[1]][1]
    cat(i, "/", total, "| Assigning neoantigens to clonal structure from", patient, "\n")
    i <<- i + 1
    ## Load and prepare PyClone result
    clones <- .preparePyCloneData(file)
    
    ## Merge neoantigens with cluster IDs
    neos <- merge(neoantigens, clones[, c("Patient", "mutation_id", "cluster_id", "cluster")], by.x = c("Patient_ID", "mutation_id"), by.y = c("Patient", "mutation_id"), all.x = TRUE)
    if("cluster_id.x" %in% names(neos)){
      neos$cluster_id.x[is.na(neos$cluster_id.x)] <- neos$cluster_id.y[is.na(neos$cluster_id.x)]
      neos$cluster.x[is.na(neos$cluster.x)] <- neos$cluster.y[is.na(neos$cluster.x)]
      neos$cluster_id.y <- neos$cluster.y <- NULL
      names(neos) <- gsub("[.]x$", "", names(neos))
    }
    neoantigens <<- neos
  })
  
  ## pre neoantigens
  neos.pre <- neoantigens[neoantigens$Pre.Post == "Pre",]
  general.save.data(neos.pre, paste0(dir.out, "Neoantigens_clusters_PreOnly"))
  
  ## For plotting remove duplicated entries - keep per patient not per sample
  neos.pre <- unique(neos.pre[, names(neos.pre) %in% c("Patient_ID", "mutation_id", "Best.MT.Score", "Corresponding.WT.Score", "DAI", "cluster_id", "cluster")])
  
  ## Plotting
  #neos.pre <- .plot.DAI.clusters(neos.pre, paste0(dir.out, "pre_DAI"))
  
  .plot.DAI.clusters <- function(neos, out){
    w = 10
    h = 10
    if(length(neos$Patient_ID) > 10){
      w = 20
      h = 30
    }
    ## Adjust cluster IDs
    neos$cluster_id <- as.factor(neos$cluster_id)
    neos$cl <- factor(paste0("C", neos$cl), levels = paste0("C", 1:max(neos$cl, na.rm = TRUE)))
    
    ## Plot DAI pre/post per clusters
    ggplot(neos, aes(x = cluster_id, y = DAI)) +
      geom_violin() + 
      geom_beeswarm(aes(color = cluster_id)) +
      stat_summary(fun = "mean", geom = "crossbar", color = "red") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      facet_wrap(. ~ Patient_ID, scale = "free_x")
    ggsave(paste0(out, "_clusterID.png"), type = "cairo", w = 10, h = 5)
    ggsave(paste0(out, "_clusterID.pdf"), w = 10, h = 5)
    
    ggplot(neos, aes(x = cl, y = DAI)) +
      geom_violin() + 
      geom_beeswarm(aes(color = cl)) +
      stat_summary(fun = "mean", geom = "crossbar", color = "red") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      facet_wrap(. ~ Patient_ID, scale = "free_x") 
    ggsave(paste0(out, "_cluster.png"), type = "cairo", w = 10, h = 5)
    ggsave(paste0(out, "_cluster.pdf"), w = 10, h = 5)
    
    ## No scaling
    ggplot(neos, aes(x = cluster_id, y = DAI)) +
      geom_violin() + 
      geom_beeswarm(aes(color = cluster_id)) +
      stat_summary(fun = "mean", geom = "crossbar", color = "red") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      facet_wrap(. ~ Patient_ID, scale = "free")
    ggsave(paste0(out, "_clusterID_noScale.png"), type = "cairo", w = 10, h = 5)
    ggsave(paste0(out, "_clusterID_noScale.pdf"), w = 10, h = 5)
    
    ggplot(neos, aes(x = cl, y = DAI)) +
      geom_violin() + 
      geom_beeswarm(aes(color = cl)) +
      stat_summary(fun = "mean", geom = "crossbar", color = "red") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      facet_wrap(. ~ Patient_ID, scale = "free") 
    ggsave(paste0(out, "_cluster_noScale.png"), type = "cairo", w = 10, h = 5)
    ggsave(paste0(out, "_cluster_noScale.pdf"), w = 10, h = 5)
    
    return(neos)
  }
  neos.pre <- .plot.DAI.clusters(neos.pre, paste0(dir.out, "pre_DAI"))
  neos.prepost <- .plot.DAI.clusters(neos.prepost, paste0(dir.out, "prePost_DAI"))
}

neoantigensPerTimepoint <- function(mode){
  dir.out <- "RESULTS/Neoantigens_Jun23/"
  .checkAndCreateDir(dir.out)
  data.allClonalSubclonal <- general.save.data(paste0("RESULTS/MutationalSignatures_May23/data.allClonalSubclonalMutations_", mode, ".csv"))
  neoantigens$mutationID <- paste(neoantigens$Gene.Name, gsub("chr", "", neoantigens$Chromosome), neoantigens$Stop, neoantigens$Stop, neoantigens$Reference, neoantigens$Variant, sep = "_")
  
  data <- merge(data.allClonalSubclonal, unique(neoantigens[, c("mutationID", "DAI")]), by = "mutationID")
  data$Classification.timepoint[data$Classification == "clonal"] <- "clonal"
  data$Classification.timepoint <- factor(data$Classification.timepoint, levels = c("clonal", "unique.pre", "unique.post"))
  data$type <- sapply(strsplit(data$Sample_ID, "_"), "[", 2)
  data$type <- factor(data$type, levels = c("clonal", "subclonal", "unique.pre", "unique.post"))
  
  my_comparisons <- list(c("clonal", "unique.pre"), c("clonal", "unique.post"), c("unique.pre", "unique.post"))
  ggplot(data, aes(x = type, y = DAI)) +
    geom_violin() + 
    geom_beeswarm() +
    stat_summary(fun = "mean", geom = "crossbar", color = "black") +
    stat_compare_means(comparisons = my_comparisons, method = "t.test") + # Add pairwise comparisons p-value
    theme_minimal()
  ggsave(paste0(dir.out, "Neoantigens_uniquePrePost.png"), type = "cairo", w = 10, h = 5)
  ggsave(paste0(dir.out, "Neoantigens_uniquePrePost.pdf"), w = 10, h = 5)
  
  ggplot(data, aes(x = type, y = DAI)) +
    geom_violin() + 
    stat_summary(fun = "mean", geom = "crossbar", color = "black") +
    stat_compare_means(comparisons = my_comparisons, method = "t.test") + # Add pairwise comparisons p-value
    theme_minimal()
  ggsave(paste0(dir.out, "Neoantigens_uniquePrePost_noDots.png"), type = "cairo", w = 5, h = 5)
  ggsave(paste0(dir.out, "Neoantigens_uniquePrePost_noDots.pdf"), w = 5, h = 5)
  
  ggplot(data, aes(x = type, y = DAI)) +
    geom_violin() + 
    stat_summary(fun = "mean", geom = "crossbar", color = "black") +
    stat_compare_means(comparisons = my_comparisons, method = "kruskal.test") + # Add pairwise comparisons p-value
    theme_minimal() +
    facet_wrap(. ~ Patient, scale = "free") 
  ggsave(paste0(dir.out, "Neoantigens_uniquePrePost_noDots_perPatient.png"), type = "cairo", w = 15, h = 10)
  ggsave(paste0(dir.out, "Neoantigens_uniquePrePost_noDots_perPatient.pdf"), w = 15, h = 10)
  
  
  data.sub <- data[which(data$type %in% c("unique.pre", "unique.post")),]
  data.sub$type <- factor(data.sub$type, levels = c("unique.pre", "unique.post"))

  ggplot(data.sub, aes(x = type, y = DAI)) +
    geom_violin() + 
    geom_beeswarm() +
    stat_summary(fun = "median", geom = "crossbar", color = "black") +
    stat_compare_means() + # Add pairwise comparisons p-value
    theme_minimal()
  ggsave(paste0(dir.out, "Neoantigens_uniquePrePost_noClonal.png"), type = "cairo", w = 5, h = 5)
  ggsave(paste0(dir.out, "Neoantigens_uniquePrePost_noClonal.pdf"), w = 5, h = 5)
  
  ggplot(data.sub, aes(x = type, y = DAI)) +
    geom_violin() + 
 #   geom_beeswarm() +
    scale_y_continuous(trans = "log10") +
    stat_summary(fun = "median", geom = "crossbar", color = "red") +
    stat_summary(fun = "mean", geom = "crossbar", color = "black") +
    stat_compare_means() + # Add pairwise comparisons p-value
    ylab("log10(DAI)") +
    theme_minimal()
  ggsave(paste0(dir.out, "Neoantigens_uniquePrePost_noClonal_log10.png"), type = "cairo", w = 5, h = 5)
  ggsave(paste0(dir.out, "Neoantigens_uniquePrePost_noClonal_log10.pdf"), w = 5, h = 5)

  ggplot(data.sub, aes(x = type, y = DAI)) +
    geom_violin() + 
    #   geom_beeswarm() +
    scale_y_continuous(trans = "log2") +
    stat_summary(fun = "median", geom = "crossbar", color = "red") +
    stat_summary(fun = "mean", geom = "crossbar", color = "black") +
    stat_compare_means() + # Add pairwise comparisons p-value
    ylab("log2(DAI)") +
    theme_minimal()
  ggsave(paste0(dir.out, "Neoantigens_uniquePrePost_noClonal_log2.png"), type = "cairo", w = 5, h = 5)
  ggsave(paste0(dir.out, "Neoantigens_uniquePrePost_noClonal_log2.pdf"), w = 5, h = 5)
  

  ggplot(data.sub, aes(x = type, y = DAI)) +
    geom_violin() + 
    geom_beeswarm() +
    theme_minimal() +
    facet_wrap(. ~ Patient, scale = "free") 
  ggsave(paste0(dir.out, "Neoantigens_uniquePrePost_noClonal_perPatient.png"), type = "cairo", w = 10, h = 5)
  ggsave(paste0(dir.out, "Neoantigens_uniquePrePost_noClonal_perPatient.pdf"), w = 10, h = 5)
  
  p <- ggboxplot(data.sub, x = "type",y = "DAI",
            facet.by = "Patient", scales = "free_y", ggtheme = theme_pubclean())
  p + stat_compare_means(
    aes(label = paste("p=", signif(as.numeric(..p.format..)))),
    method = "kruskal", size = 2)
  ggsave(paste0(dir.out, "Box_Neoantigens_uniquePrePost_noClonal_perPatient.png"), type = "cairo", w = 10, h = 5)
  ggsave(paste0(dir.out, "Box_Neoantigens_uniquePrePost_noClonal_perPatient.pdf"), w = 10, h = 5)
  
  p <- ggboxplot(data.sub, x = "type",y = "DAI", ggtheme = theme_pubclean())
  p + stat_compare_means(
    aes(label = paste("p=", signif(as.numeric(..p.format..)))),
    method = "kruskal")
  ggsave(paste0(dir.out, "Box_Neoantigens_uniquePrePost_noClonal.png"), type = "cairo", w = 10, h = 5)
  ggsave(paste0(dir.out, "Box_Neoantigens_uniquePrePost_noClonal.pdf"), w = 10, h = 5)
  
}




