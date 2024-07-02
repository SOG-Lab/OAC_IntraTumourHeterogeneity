cat("*** Project ID: OAC_Heterogeneity_MR_MT ***\n")
## R script to analyse the single nucleotide mutation signatures using COSMIC v3.3 SBS ##
cat("Loading: R script to analyse the single nucleotide mutation signatures using COSMIC v3.3 SBS\n")

## Loading external libraries ##
library(deconstructSigs) ## Get context of each mutation
library(ggpubr) ## ggarrange() to combine ggplots
library(survival) ## Survival analysis
library(survminer) ## Kaplan-Meier plots
library(ggstatsplot) ## Bar plots with stats
library(rstatix) ## p-value stats per facet
## Install SignatureEstimation package from source (from https://rdrr.io/github/WuyangFF95/SynSigRun/src/R/RunSignatureEstimationAttributionOnly.R)
## remotes::install_url("https://www.ncbi.nlm.nih.gov/CBBresearch/Przytycka/software/signatureestimation/SignatureEstimation.tar.gz")
library(SignatureEstimation) ## Use quadratic programing approach to get mutational signatures

## Loading external scripts ##
## Load script with general project settings and background functionality
source("/SCRIPTS/setup.R")

## Environment ##
mutationSignatures.file <- "DATA/COSMIC_v3.3.1_SBS_GRCh37.txt"
insert <- ""
mode <- "preOnly"
dir.out <- "RESULTS/MutationalSignatures/"
.checkAndCreateDir(dir.out)
.checkAndCreateDir(paste0(dir.out, "Signatures/"))
.checkAndCreateDir(paste0(dir.out, "Plots/"))
.checkAndCreateDir(paste0(dir.out, "KM_Plots/"))
cohortSigs <- paste0("SBS", c(2, 13, 3, 36, 1, 5, 18, 20, 26, 44, 8, 39, 40, 93, "17a", "17b")) ## In plotting order
## Substitution types
substitution.types <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
## Signature contribution cut-off
cut <- 0.05

.simplifySubstitutions <- function(data){
  data.simp <- data.frame()
  a <- sapply(substitution.types, function(typ){
    tmp <- rowSums(data[, grepl(typ, colnames(data))])
    tmp <- data.frame(Sample_ID = names(tmp), typ = as.numeric(tmp))
    names(tmp) <- c("Sample_ID", typ)
    if(nrow(data.simp) == 0){
      data.simp <<- tmp
    }else{
      data.simp <<- merge(data.simp, tmp)
    }
  })
  return(data.simp)
}

## Visualise mutational signatures as barcharts with substitution types and number of mutations
.plotMutSigsBarchartBig <- function(data, out){
  data <- data[data$Percent > 0,]
  
  ## Get signatures in numerical order
#  data$Signature <- factor(data$Signature, levels = names(mutSigs)[names(mutSigs) %in% unique(data$Signature)])
  if(any(c("SBS31", "SBS35") %in% data$Signature)){
    data$Signature <- factor(data$Signature, levels = unique(c(cohortSigs[cohortSigs %in% unique(data$Signature)], c("SBS31", "SBS35"))))
  }else{
    data$Signature <- factor(data$Signature, levels = cohortSigs[cohortSigs %in% unique(data$Signature)])
  }
  
  if(!any(names(data) == "Timepoint")){
    data$Timepoint <- "pre"
  }else{
    if("pre" %in% data$Timepoint){
      color.time <- c("#DAB390")
    }
    if("post" %in% data$Timepoint){
      color.time <- c(color.time, "#9E042C")
    }
    if("clonal" %in% data$Timepoint){
      color.time <- c(color.time, "#595959")
    }
    if("subclonal" %in% data$Timepoint){
      color.time <- c(color.time, "#000000")
      names(color.time) <- c("pre", "post", "clonal", "subclonal")
      data$Timepoint <- factor(data$Timepoint, levels = c("pre", "post", "clonal", "subclonal"))
    }else{
      data$Timepoint <- factor(data$Timepoint, levels = c("pre", "post", "clonal"))
    }
  }
  
  mutbars <- ggplot(unique(data[, c("Sample_ID", "MutCount", "Timepoint")]), aes(x = Sample_ID, y = MutCount, fill = Timepoint)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(x = "", y = "Mutation count") + 
    scale_fill_manual(values = color.time) + 
    theme_test() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  
  tmp <- unique(data[, c("Sample_ID", substitution.types)])
  tmp <- as.data.frame(pivot_longer(tmp, names_to = "Substitution", values_to = "Proportion", cols = substitution.types))
  subs <- ggplot(tmp, aes(x = Sample_ID, y = Proportion, fill = Substitution)) + 
    geom_bar(stat = "identity", color = "black", position = "fill") +
    theme_minimal() +
    labs(x = "", y = "Proportion [%]") +
    scale_fill_manual(values = c("#627E4A", "#EAECD1", "#EFE3A0", "#C3A579", "#913D2C", "#4A1415")) + 
    theme_test() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  sigs <- ggplot(data, aes(x = Sample_ID, y = Percent, fill = Signature)) + 
    geom_bar(stat = "identity", color = "white") +
    theme_minimal() +
    labs(x = "", y = "Proportion [%]") +
    theme_test() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  figure <- ggarrange(mutbars, subs, sigs, ncol = 1, nrow = 3, common.legend = FALSE, legend = "right")
  
  ggsave(paste0(out, ".pdf"), width = 30, height = 15)
  ggsave(paste0(out, ".png"), width = 30, height = 15, type = "cairo")
}


assignMutationalSignatures <- function(mode, insert = ""){
  mutSigs <- general.load.data(mutationSignatures.file, row.names = 1)
  ## Remove sequencing artefect signatures
  mutSigs <- mutSigs[, !names(mutSigs) %in% signatures.seqArtefacts]
  cat("Assigning mutational signatures for", mode, "mode.\n\n")
  ## Load data if available
  if(mode == "preOnly"){
    if(exists("allMutations.preOnly")){
      mutationCollection.file <- allMutations.preOnly
    }else{
      mutationCollection.file <- paste0("DATA/allMutations_collection_preOnly", insert, ".RData")
      if(!file.exists(mutationCollection.file)){
       stop("File doesn't exist. Please run prepareMutationCollections() from the mutationalAnalysis.R script first.")
      }
    }
  }else if(mode == "prePost"){
    if(exists("allMutations.preOnly")){
      mutationCollection.file <- allMutations.prePost
    }else{
      mutationCollection.file <- paste0("DATA/allMutations_collection_prePost", insert, ".RData")
      if(!file.exists(mutationCollection.file)){
        stop("File doesn't exist. Please run prepareMutationCollections() from the mutationalAnalysis.R script first.")
      }
    }
    ## Add platinum-based treatment signatures for pre/post analysis
    cohortSigs <- c(cohortSigs, "SBS31", "SBS35") ## Add platinum-based signatures when comparing pre/post samples
  }
  
  if(!is.list(mutationCollection.file)){
    res <- general.load.data(mutationCollection.file)
  }else{
    res <- mutationCollection.file
  }
  ## Load annotation data
  cohort <- general.load.data(paste0(dir.data, "Annotation/Cohort.csv")) ## Contains sample IDs and clinical information for each patient
  cohort.files <- .prepareCohortFiles(paste0(dir.data, "Annotation/Cohort_sampleDetails.csv")) ## Contains sample information and file storage locations
  
  ## Keep only patients/samples based on mode
  cohort <- cohort[cohort$Patient %in% names(res), ]
  cohort.files <- cohort.files[cohort.files$Patient %in% names(res), ]
  
  ## Store the assigned mutation signatures and single nucleotid substitutions per sample
  signatures.unfiltered <- signatures.SigReassigned <- substitutions <- mutCounts <- data.frame()
  data.allClonalSubclonal <- data.frame()
  ## Assign mutational signatures to each sample per patient
  i <- 1
  total <- length(res)
  a <- sapply(res[cohort$Patient], function(dataList){
    patient <- as.character(dataList$patient)
    cat("\n", as.character(Sys.time()), "-", i, "/", total, "|", patient)
    i <<- i + 1
    ## Get the required data
    data <- dataList$allMutations
    ## Only SNP for this type of analysis
    data <- data[data$Variant_Type == "SNP",]
    if(mode == "preOnly"){
      data <- data[!data$Sample_ID %in% exclude.pre,]
      data$Sample_ID <- ifelse(data$Classification == "clonal", paste0(patient, "_clonal"), paste0(patient, "_unique"))
    }else if(mode == "prePost"){
      data$Sample_ID <- ifelse(data$Classification.timepoint == "unique.pre", paste0(patient, "_unique.pre"), paste0(patient, "_unique.post"))
      data$Sample_ID[which(data$Classification == "clonal")] <- paste0(patient, "_clonal")
      ## Remaining mutations (present in some pre and some post samples)
      data$Sample_ID[is.na(data$Sample_ID)] <- paste0(patient, "_subclonal")
      data <- data[!is.na(data$Sample_ID), ] ## Remove subclonal mutations across pre/post samples 
    }else{
      stop("Please select valid mode (preOnly/prePost).\n")
    }
    ## Remove duplicated entries
    data <- unique(data[order(data$Sample_ID),])
    data$Patient <- patient
    if(nrow(data.allClonalSubclonal) == 0){
      data.allClonalSubclonal <<- data
    }else{
      data.allClonalSubclonal <<- rbind(data.allClonalSubclonal, data)
    }
    ## Get the single nucleotid substitutions based on the mutation position (default is hg19 assembly)
    data.sigs <- mut.to.sigs.input(mut.ref = data, sample.id = "Sample_ID", chr = "Chromosome", pos = "Start_Position", ref = "Reference_Allele", alt = "Tumor_Seq_Allele2")
  
    ## Save substitutions per sample
    if(nrow(substitutions) == 0){
      substitutions <<- data.sigs
    }else{
      substitutions <<- rbind(substitutions, data.sigs)
    }
    
    ## Number of mutations considered for mutational signature analysis
    mutCount <- rowSums(data.sigs, na.rm = TRUE)
    if(nrow(mutCounts) == 0){
      mutCounts <<- data.frame(Sample_ID = names(mutCount), MutCount = as.numeric(mutCount))
    }else{
      mutCounts <<- rbind(mutCounts, data.frame(Sample_ID = names(mutCount), MutCount = as.numeric(mutCount)))
    }
    mutSigs <- mutSigs[names(data.sigs), cohortSigs]
    
    a <- sapply(row.names(data.sigs), function(sample){
      cat(" |", sample)
      ## Assign the mutational signatures # signatures.ref = mutSigs
      tmp <- findSigExposures(M = t(data.sigs[sample,]), P = mutSigs[, cohortSigs], decomposition.method = decomposeQP)$exposures
      if(nrow(signatures.unfiltered) == 0){
        signatures.unfiltered <<- tmp
      }else{
        signatures.unfiltered <<- cbind(signatures.unfiltered, tmp)
      }
      ## To avoid over-fitting, signatures < 10% for each sample were omitted and mutations were reassigned to the remaining signatures.
      selectedSigs <- row.names(tmp[tmp > cut, , drop = FALSE])
      
      if(length(selectedSigs) == 0){
        warning(paste0("No signatures for ", sample, ". Cohort signatures used."))
        selectedSigs <- cohortSigs
        tmp <- findSigExposures(M = t(data.sigs[sample,]), P = mutSigs[, selectedSigs], decomposition.method = decomposeQP)$exposures
        selectedSigs <- row.names(tmp[tmp > cut, , drop = FALSE])
      }
      if(length(selectedSigs) == 1){
        warning(paste0("Only 1 signature passed the threshold for ", sample, ". OAC signatures (SBS17a/b) used."))
        selectedSigs <- unique(c(selectedSigs, "SBS17a", "SBS17b"))
      }
      removedSigs <- names(mutSigs)[!names(mutSigs) %in% selectedSigs]
      tmp.SigReassigned <- findSigExposures(M = t(data.sigs[sample,]), P = mutSigs[, selectedSigs], decomposition.method = decomposeQP)$exposures
      
      tmp.SigReassigned <- rbind(tmp.SigReassigned, matrix(rep(0, length(removedSigs)), dimnames = list(removedSigs, sample)))
      tmp.SigReassigned <- tmp.SigReassigned[names(mutSigs), , drop = FALSE]
      if(nrow(signatures.SigReassigned) == 0){
        signatures.SigReassigned <<- tmp.SigReassigned
      }else{
        signatures.SigReassigned <<- cbind(signatures.SigReassigned, tmp.SigReassigned)
      }
    }) 
  })
  cat("\n")
  ## Save all analysed mutations
  general.save.data(data.allClonalSubclonal, paste0(dir.out, "data.allClonalSubclonalMutations_", mode))
  
  ## Change to long format
  signatures.SigReassigned <- as.data.frame(t(as.data.frame(signatures.SigReassigned)))
  signatures.unfiltered <- as.data.frame(t(as.data.frame(signatures.unfiltered)))
  
  ## Visualisation
  ## Change to long format for plotting
  signatures.unfiltered$Sample_ID <- row.names(signatures.unfiltered)
  signatures.SigReassigned$Sample_ID <- row.names(signatures.SigReassigned)
  
  ## Save data
  general.save.data(signatures.SigReassigned, paste0(dir.out, "Signatures/signatures.SigReassigned_", mode, insert, "_", cut), row.names = TRUE)
  general.save.data(signatures.unfiltered, paste0(dir.out, "Signatures/signatures.unfiltered_", mode, insert), row.names = TRUE)
  
  data.plot <- as.data.frame(pivot_longer(signatures.SigReassigned, names_to = "Signature", values_to = "Percent", cols = starts_with("SBS")))
  data.plot.unfiltered <- signatures.unfiltered
#  data.plot.unfiltered$SBSother <- 0
#  data.plot.unfiltered[data.plot.unfiltered < cut] <- 0
#  data.plot.unfiltered$SBSother <- 1 - rowSums(data.plot.unfiltered[, cohortSigs], na.rm = TRUE)
  data.plot.unfiltered <- as.data.frame(pivot_longer(data.plot.unfiltered, names_to = "Signature", values_to = "Percent", cols = starts_with("SBS")))

  data.plot$Timepoint <- ifelse(grepl("clonal", data.plot$Sample_ID), "clonal", "pre")
  data.plot.unfiltered$Timepoint <- ifelse(grepl("clonal", data.plot.unfiltered$Sample_ID), "clonal", "pre")
  
  subs <- .simplifySubstitutions(substitutions)
  data.plot <- merge(merge(data.plot, subs, all.x = TRUE), mutCounts)
  data.plot <- data.plot[data.plot$Percent > 0,]
  data.plot.unfiltered <- merge(merge(data.plot.unfiltered, subs, all.x = TRUE), mutCounts)
  data.plot.unfiltered <- data.plot.unfiltered[data.plot.unfiltered$Percent > 0,]
  
  if(mode == "prePost"){
    data.plot$Timepoint[grepl("post", data.plot$Sample_ID)] <- "post"
    data.plot.unfiltered$Timepoint[grepl("post", data.plot.unfiltered$Sample_ID)] <- "post"
    if(any(grepl("subclonal", data.plot$Sample_ID))){
      data.plot$Timepoint[grepl("subclonal", data.plot$Sample_ID)] <- "subclonal"     
      data.plot$Sample_ID <- factor(data.plot$Sample_ID, levels = paste0(rep(cohort$Patient, each = 4), c("_clonal", "_subclonal", "_unique.pre", "_unique.post")))
      data.plot.unfiltered$Timepoint[grepl("subclonal", data.plot.unfiltered$Sample_ID)] <- "subclonal" 
      data.plot.unfiltered$Sample_ID <- factor(data.plot.unfiltered$Sample_ID, levels = paste0(rep(cohort$Patient, each = 4), c("_clonal", "_subclonal", "_unique.pre", "_unique.post")))
    }else{
      data.plot$Sample_ID <- factor(data.plot$Sample_ID, levels = paste0(rep(cohort$Patient, each = 3), c("_clonal", "_unique.pre", "_unique.post")))
      data.plot.unfiltered$Sample_ID <- factor(data.plot.unfiltered$Sample_ID, levels = paste0(rep(cohort$Patient, each = 3), c("_clonal", "_unique.pre", "_unique.post")))
    }
  }else{
    patients <- cohort$Patient[order(cohort$OS, max(cohort$OS.endo)-cohort$OS.endo, decreasing = FALSE)]
    data.plot$Sample_ID <- factor(data.plot$Sample_ID, levels = paste0(rep(patients, each = 2), c("_clonal", "_unique")))
    data.plot.unfiltered$Sample_ID <- factor(data.plot.unfiltered$Sample_ID, levels = paste0(rep(patients, each = 2), c("_clonal", "_unique")))
  }
  .plotMutSigsBarchartBig(data.plot, paste0(dir.out, "SignatureEstimation_BIG_barchart_uniClonal_Mutations_", mode, insert, "_", cut))
  .plotMutSigsBarchartBig(data.plot.unfiltered, paste0(dir.out, "SignatureEstimation_BIG_barchart_uniClonal_Mutations_", mode, insert, "_unfiltered"))
  
#  signatures.SigReassigned <- signatures.SigReassigned[, c("Sample_ID", cohortSigs)]
  signatures.SigReassigned <- merge(signatures.SigReassigned, subs)
  signatures.SigReassigned <- merge(signatures.SigReassigned, mutCounts)
  ## Add patient IDs
  signatures.SigReassigned$Patient <- gsub("_clonal", "", gsub("_subclonal", "", gsub("_unique.post", "", gsub("_unique.pre", "", signatures.SigReassigned$Sample_ID))))
  ## Add mutation groups
  signatures.SigReassigned$Category <- "clonal"
  signatures.SigReassigned$Category[grepl("subclonal", signatures.SigReassigned$Sample_ID)] <- "subclonal"
  signatures.SigReassigned$Category[grepl("unique.post", signatures.SigReassigned$Sample_ID)] <- "unique.post"
  signatures.SigReassigned$Category[grepl("unique.pre", signatures.SigReassigned$Sample_ID)] <- "unique.pre"
  ## Add treatment
  signatures.SigReassigned <- merge(signatures.SigReassigned, cohort[, c("Patient", "OS", "OS.endo", "Allocated.treatment")])
  signatures.SigReassigned$Treatment <- gsub("Arm [AB]: ", "", signatures.SigReassigned$Allocated.treatment)
  
  general.save.data(signatures.SigReassigned, paste0(dir.out, "Data_SignatureEstimation_BIG_barchart_uniClonal_Mutations_", mode, insert, "_", cut))
  ## Violin plot of mutation groups
  signatures.SigReassigned <- signatures.SigReassigned[signatures.SigReassigned$Category != "subclonal", ]
  signatures.SigReassigned$Category <- factor(signatures.SigReassigned$Category, levels = c("clonal", "subclonal", "unique.pre", "unique.post"))
  my_comparisons <- list(c("clonal", "unique.pre"), c("clonal", "unique.post"), c("unique.pre", "unique.post") )
  ggplot(signatures.SigReassigned, aes(x = Category, y = MutCount)) +
    geom_violin() + 
    stat_summary(fun = "mean", geom = "crossbar", color = "black") + 
    stat_compare_means(comparisons = my_comparisons, method = "t.test") + # Add pairwise comparisons p-value
    geom_beeswarm()
  ggsave(paste0(dir.out, "Violin_MutationCategories_", mode, insert, ".png"), width = 5, height = 4, type = "cairo")
  ggsave(paste0(dir.out, "Violin_MutationCategories_", mode, insert, ".pdf"), width = 5, height = 4)
  
  ggplot(signatures.SigReassigned, aes(x = Category, y = MutCount)) +
    geom_violin() + 
    geom_line(aes(group = Patient, color = Patient)) +
    geom_beeswarm(aes(color = Patient)) +
    stat_summary(fun = "mean", geom = "crossbar", color = "black") + 
    stat_compare_means(comparisons = my_comparisons, method = "t.test", paired = TRUE) + # Add pairwise comparisons p-value
    theme_test()
  ggsave(paste0(dir.out, "Violin_MutationCategories_paired_", mode, insert, ".png"), width = 5, height = 4, type = "cairo")
  ggsave(paste0(dir.out, "Violin_MutationCategories_paired_", mode, insert, ".pdf"), width = 5, height = 4)
  
}

signatureViolinPlots <- function(mode, sigMode = "reassigned", insert = ""){
  ## Load data
  ## Signature proportions
  if(sigMode == "reassigned"){
    insert <- paste0(insert, "_", cut)
    data <- general.load.data(paste0(dir.out, "Signatures/signatures.SigReassigned_", mode, insert, ".RData"))
  }else if(sigMode == "unfiltered"){
    insert <- paste0(insert, "_", cut)
    data <- general.load.data(paste0(dir.out, "Signatures/signatures.unfiltered_", mode, insert, ".RData"))
    if("Sample_ID" %in% names(data)){
      samples <- data$Sample_ID
      data$Sample_ID <- NULL
    }
    data$SBSother <- 0
    data[data < cut] <- 0
    data$SBSother <- 1 - rowSums(data[, cohortSigs], na.rm = TRUE)
    data$Sample_ID <- row.names(data)    
  }
  dir.plot <- paste0(dir.out, "Plots/")
  
  ## Add combined signatures (e.g. SBS17a/b - OAC or SBS2/13 - APOBEC)
  data$SBSapo <- data$SBS2 + data$SBS13
  data$SBSoac <- data$SBS17a + data$SBS17b
  
  ## Remove signatures not present in cohort
  data <- data[, colSums(data[, grepl("SBS", names(data))]) > 0]
  signatures <- names(data)[grepl("SBS", names(data))]
  
  signatures.plot <- as.data.frame(pivot_longer(data, names_to = "Signature", values_to = "Proportion", cols = signatures))
  signatures.plot$Category <- ifelse(grepl("clonal", signatures.plot$Sample_ID), "clonal", "unique")
  signatures.plot$Signature <- factor(signatures.plot$Signature, levels = signatures)
  
  stat.test <- signatures.plot %>%
    group_by(Signature) %>%
    t_test(Proportion ~ Category, paired = TRUE) 
  stat.test <- stat.test %>% add_xy_position(x = "Category")
  ggp <- ggplot(signatures.plot, aes(x = Category, y = Proportion)) +
    geom_violin() + 
    stat_summary(fun = "mean", geom = "crossbar", color = "black") + 
    geom_beeswarm() +
    facet_grid(. ~ Signature) +
    stat_pvalue_manual(stat.test)
  ggsave(paste0(dir.plot, "Violin_Signatures_", mode, insert, ".png"), width = 18, height = 4, type = "cairo")
  ggsave(paste0(dir.plot, "Violin_Signatures_", mode, insert, ".pdf"), width = 18, height = 4)
  
  sigs <- paste0("SBS", c(2, "apo", 3, 13, 36))
  stat.test <- signatures.plot[signatures.plot$Signature %in% sigs,] %>%
    group_by(Signature) %>%
    t_test(Proportion ~ Category, paired = TRUE) 
  stat.test <- stat.test %>% add_xy_position(x = "Category")
  ggplot(signatures.plot[signatures.plot$Signature %in% sigs,], aes(x = Category, y = Proportion)) +
    geom_violin() + 
    stat_summary(fun = "mean", geom = "crossbar", color = "black") + 
    geom_beeswarm() +
    facet_grid(. ~ Signature) +
    stat_pvalue_manual(stat.test)
  ggsave(paste0(dir.plot, "Violin_Signatures_uniqueHigh_", mode, insert, ".png"), width = 4, height = 4, type = "cairo")
  ggsave(paste0(dir.plot, "Violin_Signatures_uniqueHigh_", mode, insert, ".pdf"), width = 4, height = 4)
  
  sigs <- paste0("SBS", c("17a", "17b", "oac"))
  stat.test <- signatures.plot[signatures.plot$Signature %in% sigs,] %>%
    group_by(Signature) %>%
    t_test(Proportion ~ Category, paired = TRUE) 
  stat.test <- stat.test %>% add_xy_position(x = "Category")
  ggplot(signatures.plot[signatures.plot$Signature %in% sigs,], aes(x = Category, y = Proportion)) +
    geom_violin() + 
    stat_summary(fun = "mean", geom = "crossbar", color = "black") + 
    geom_beeswarm() +
    facet_grid(. ~ Signature) +
    stat_pvalue_manual(stat.test)
  ggsave(paste0(dir.plot, "Violin_Signatures_clonalHigh_", mode, insert, ".png"), width = 3, height = 4, type = "cairo")
  ggsave(paste0(dir.plot, "Violin_Signatures_clonalHigh_", mode, insert, ".pdf"), width = 3, height = 4)
  
  ## p-values
  a <- sapply(signatures, function(sig){
    p <- NA
    ## If there are NAs in the percent then can't do paired
    if(any(is.na(signatures.plot$Proportion))){
      try(p <- t.test(signatures.plot$Proportion[signatures.plot$Category == "clonal" & signatures.plot$Signature == sig], signatures.plot$Proportion[signatures.plot$Category == "unique" & signatures.plot$Signature == sig], paired = FALSE)$p.value, silent = TRUE)
    }else{
      try(p <- t.test(signatures.plot$Proportion[signatures.plot$Category == "clonal" & signatures.plot$Signature == sig], signatures.plot$Proportion[signatures.plot$Category == "unique" & signatures.plot$Signature == sig], paired = TRUE)$p.value, silent = TRUE)
    }
    cat(sig, "- p =", signif(p, digits = 4), "\n")
    ## Combine SBS17a and SBS17b percentages
    if(sig == "SBS17b"){
      clonal.17 <- signatures.plot$Proportion[signatures.plot$Category == "clonal" & signatures.plot$Signature == "SBS17a"] + signatures.plot$Proportion[signatures.plot$Category == "clonal" & signatures.plot$Signature == "SBS17b"]
      unique.17 <- signatures.plot$Proportion[signatures.plot$Category == "unique" & signatures.plot$Signature == "SBS17a"] + signatures.plot$Proportion[signatures.plot$Category == "unique" & signatures.plot$Signature == "SBS17b"]
      if(any(is.na(signatures.plot$Proportion))){
        try(p <- t.test(clonal.17, unique.17, paired = FALSE)$p.value, silent = TRUE)
      }else{
        try(p <- t.test(clonal.17, unique.17, paired = TRUE)$p.value, silent = TRUE)
      }
      cat("SBS17a/b - p =", signif(p, digits = 4), "\n")
    }else if(sig == "SBS2"){
      ## Combine SBS2 and SBS13 percentages (APOBEC signatures)
      clonal.apo <- signatures.plot$Proportion[signatures.plot$Category == "clonal" & signatures.plot$Signature == "SBS2"] + signatures.plot$Proportion[signatures.plot$Category == "clonal" & signatures.plot$Signature == "SBS13"]
      unique.apo <- signatures.plot$Proportion[signatures.plot$Category == "unique" & signatures.plot$Signature == "SBS2"] + signatures.plot$Proportion[signatures.plot$Category == "unique" & signatures.plot$Signature == "SBS13"]
      if(any(is.na(signatures.plot$Proportion))){
        try(p <- t.test(clonal.apo, unique.apo, paired = FALSE)$p.value, silent = TRUE)
      }else{
        try(p <- t.test(clonal.apo, unique.apo, paired = TRUE)$p.value, silent = TRUE)
      }
    }
  })
  
  ## Present/absent tables for fisher's exact test
  ## https://statsandr.com/blog/fisher-s-exact-test-in-r-independence-test-for-a-small-sample/
  data$Category <- ifelse(grepl("clonal", data$Sample_ID), "clonal", "unique")
  signatures <- c(signatures, "SBSapoBRCAMUTYH", "SBSapoBRCA")
  a <- sapply(signatures, function(sig){
    if(sig == "SBSapoBRCAMUTYH"){
      data$Sig <- ifelse(data[, "SBS2"] > 0 | data[, "SBS13"] > 0 | data[, "SBS3"] > 0 | data[, "SBS36"] > 0, "present", "absent")
    }else if(sig == "SBSapoBRCA"){
      data$Sig <- ifelse(data[, "SBS2"] > 0 | data[, "SBS13"] > 0 | data[, "SBS3"] > 0, "present", "absent")
    }else{
      data$Sig <- ifelse(data[, sig] > 0, "present", "absent")
    }
    tab <- table(data$Category, data$Sig)
    p <- signif(fisher.test(tab)$p.value, 2)
    ggbarstats(
      data, Sig, Category,
      perc.k = 2L, ## 2 decimals for percentage labels
      results.subtitle = FALSE,
      subtitle = paste0(
        "Fisher's exact test", ", p-value = ",
        ifelse(p < 0.001, "< 0.001", round(p, 3))
      ),
      legend.title = sig
    )
    ggsave(paste0(dir.plot, "Barplot_Signatures_presentAbsent_", mode, insert, "_", sig, ".png"), width = 3, height = 3, type = "cairo")
    ggsave(paste0(dir.plot, "Barplot_Signatures_presentAbsent_", mode, insert, "_", sig, ".pdf"), width = 3, height = 3)
  })  
}

signaturePresentAbsent_clinical <- function(mode, sigMode = "reassigned", insert = ""){
#if(TRUE){
  dir.plot <- paste0(dir.out, "KM_Plots/")
  .checkAndCreateDir(dir.plot)
  ## Go through signatures present in either the clonal or unique mutations and plot Kaplan-Meier plots
  ## Load data
  ## Signature proportions
  if(sigMode == "reassigned"){
    insert <- paste0(insert, "_", cut)
    signatures.SigReassigned <- general.load.data(paste0(dir.out, "Signatures/signatures.SigReassigned_", mode, insert, ".RData"))
  }else if(sigMode == "unfiltered"){
    insert <- paste0(insert, "_", cut)
    signatures.SigReassigned <- general.load.data(paste0(dir.out, "Signatures/signatures.unfiltered_", mode, insert, ".RData"))
    if("Sample_ID" %in% names(signatures.SigReassigned)){
      samples <- signatures.SigReassigned$Sample_ID
      signatures.SigReassigned$Sample_ID <- NULL
    }
    signatures.SigReassigned$SBSother <- 0
    signatures.SigReassigned[signatures.SigReassigned < cut] <- 0
    signatures.SigReassigned$SBSother <- 1 - rowSums(signatures.SigReassigned[, cohortSigs], na.rm = TRUE)
    signatures.SigReassigned$Sample_ID <- row.names(signatures.SigReassigned)    
  }
  ## Clinial cohort data
  cohort <- general.load.data(paste0(dir.data, "Annotation/Cohort.csv")) ## Contains sample IDs and clinical information for each patient

  ## Add combined signatures to signature list
  signatures.SigReassigned$SBSapo <- signatures.SigReassigned$SBS2 + signatures.SigReassigned$SBS13
  signatures.SigReassigned$SBSoac <- signatures.SigReassigned$SBS17a + signatures.SigReassigned$SBS17b
  signatures.SigReassigned$SBSapoBRCAMUTYH <- signatures.SigReassigned$SBS2 + signatures.SigReassigned$SBS13 + signatures.SigReassigned$SBS3 + signatures.SigReassigned$SBS36
  signatures.SigReassigned$SBSapoBRCA <- signatures.SigReassigned$SBS2 + signatures.SigReassigned$SBS13 + signatures.SigReassigned$SBS3
  
  ## Remove signatures not present in cohort
  signatures.SigReassigned <- signatures.SigReassigned[, colSums(signatures.SigReassigned[, grepl("SBS", names(signatures.SigReassigned))]) > 0]
  signatures <- names(signatures.SigReassigned)[grepl("SBS", names(signatures.SigReassigned))]

  ## Merge signature with patient data
  cohort <- merge(cohort, data.frame(Patient = rep(cohort$Patient, each = 2), Sample_ID = paste0(rep(cohort$Patient, each = 2), c("_clonal", "_unique"))))
  cohort <- merge(cohort, signatures.SigReassigned, all.x = TRUE)
  cohort$Category <- ifelse(grepl("_clonal", cohort$Sample_ID), "clonal", "unique")  
  
  ## Use present/absent signatures to group patients
  i <- 1
  total <- length(signatures)
  cut <- 0.0001 ## 0.01 and 0.15 used currently to play around with
  a <- sapply(signatures, function(sig){
    cat("\n", i, "/", total, "|", sig)
    i <<- i + 1
    ## Signature present/absent in clonal mutations
    tmp <- cohort[cohort$Category == "clonal",]
    tmp$Signature <- tmp[,sig]
    tmp$SigPresent <- ifelse(tmp$Signature >= cut, 1, 0)
    .KaplanMeierAnalysis(tmp, paste0(dir.plot, "KM_OS_SubsetReassigned_PresentAbsent_clonal_", sig, "_", cut, "_", mode, insert), sig)
    ## Signature present/absent in unique mutations
    tmp <- cohort[cohort$Category == "unique",]
    tmp$Signature <- tmp[, sig]
    tmp$SigPresent <- ifelse(tmp$Signature >= cut, 1, 0)
    .KaplanMeierAnalysis(tmp, paste0(dir.plot, "KM_OS_SubsetReassigned_PresentAbsent_unique_", sig, "_", cut, "_", mode, insert), sig)
    ## Signature present/absent in any 
    tmp <- cohort
    tmp$Signature <- tmp[, sig]
    tmp$SigPresent <- ifelse(tmp$Signature >= cut, 1, 0)
    tmp <- unique(tmp[, c("Patient", "SigPresent", "OS.endo", "OS", "PFS.endo", "PFS")])
    .KaplanMeierAnalysis(tmp, paste0(dir.plot, "KM_OS_SubsetReassigned_PresentAbsent_any_", sig, "_", cut, "_", mode, insert), sig)
  })
  cat("\n")
}

.KaplanMeierAnalysis <- function(tmp, out, sig){
  if(grepl("clonal", out)){
    cat(" | clonal - ")
    title <- "clonal"
  }else if(grepl("unique", out)){
    cat(" | unique - ")
    title <- "unique"
  }else{
    cat(" | any - ")
    title <- "any"
  }
  tmp$SigPresent <- as.factor(tmp$SigPresent)
  ## Signature percentage as continuous marker
  if(all(tmp$Signature == 0) | is.null(tmp$Signature)){
    cat("p = NA (OS) ")
  }else{
    cox <- coxph(Surv(OS.endo, OS) ~ Signature, data = tmp) 
    cat("p =", .extractPvalue(cox, "coxph"), "(OS) ")
  }
  ## Only plot Kaplan-Meier if enough patients per group (at least 25% of overall patients)
  ## OS
  fit <- survfit(Surv(OS.endo, OS) ~ SigPresent, data = tmp)
  ggsurv <- ggsurvplot(fit, data = tmp, 
                       break.time.by = 12,
                       pval = TRUE, 
                       risk.table = TRUE)
  ggsurv$plot <- ggsurv$plot + labs(title = paste("Signature", sig, "present/absent in ", title, " mutations"))
  png(paste0(out, ".png"), type = "cairo", width = 5*200, height = 5*200, res = 200)
  print(ggsurv)
  dev.off()
  pdf(paste0(out, ".pdf"), width = 5, height = 5)
  print(ggsurv)
  dev.off()
  ## HR
  if(sig %in% c("SBSapoBRCAMUTYH", "SBSapoBRCA", "SBS3") & grepl("unique", out)){
    print(table(tmp$cTNM.8th.AJCC))
    cox <- coxph(Surv(OS.endo, OS) ~ SigPresent + cTNM.8th.AJCC, data = tmp[tmp$cTNM.8th.AJCC %in% c("IIB", "III"),])
    te <<- tmp
    png(paste0(out, "_HR.png"), type = "cairo", width = 5*200, height = 2*200, res = 200)
    print(ggforest(cox, data = tmp[tmp$cTNM.8th.AJCC %in% c("IIB", "III"),]))
    dev.off()
    pdf(paste0(out, "_HR.pdf"), width = 5, height = 2)
    print(ggforest(cox, data = tmp[tmp$cTNM.8th.AJCC %in% c("IIB", "III"),]))
    dev.off()
    
    tmp$cTNM.8th.AJCC.2groups <- gsub("^III$", "III/IV", gsub("^IVA$", "III/IV",gsub("^IIB$", "I/II", gsub("^IB$", "I/II", tmp$cTNM.8th.AJCC))))
    cox <- coxph(Surv(OS.endo, OS) ~ SigPresent + cTNM.8th.AJCC.2groups, data = tmp)
    print(cox)
    te <<- tmp
    png(paste0(out, "_HR_2groups.png"), type = "cairo", width = 5*200, height = 2*200, res = 200)
    print(ggforest(cox, data = tmp))
    dev.off()
    pdf(paste0(out, "_HR_2groups.pdf"), width = 5, height = 2)
    print(ggforest(cox, data = tmp))
    dev.off()
  }
  ## PFS
  if(all(tmp$Signature == 0) | is.null(tmp$Signature)){
    cat("p = NA (OS)")
  }else{
    cat("p =", .extractPvalue(cox, "coxph"), "(PFS)")
    cox <- coxph(Surv(PFS.endo, PFS) ~ Signature, data = tmp) 
  }
  fit <- survfit(Surv(PFS.endo, PFS) ~ SigPresent, data = tmp)
  ggsurv <- ggsurvplot(fit, data = tmp, 
                       break.time.by = 12,
                       pval = TRUE, 
                       risk.table = TRUE)
  ggsurv$plot <- ggsurv$plot + labs(title = paste("Signature", sig, "present/absent in ", title, " mutations"))
  png(paste0(gsub("OS", "PFS", out), ".png"), type = "cairo", width = 5*200, height = 5*200, res = 200)
  print(ggsurv)
  dev.off()
  pdf(paste0(gsub("OS", "PFS", out), ".pdf"), width = 5, height = 5)
  print(ggsurv)
  dev.off()
}
