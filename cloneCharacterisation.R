cat("Project ID: OAC_Heterogeneity_MR_MT\n")
## R script to characterise clones identified by PyClone ##
cat("Loading: R script to characterise clones identified by PyClone\n")

## Loading external libraries ##
library(deconstructSigs) ## Assign mutation signatures to samples
library(survival) ## Survival analysis (vs clone number)
library(survminer) ## Kaplan-Meier plots for survival visualisation
library(ggbeeswarm) ## Nicer distribution of dots on ggplot
library(ggpubr) ## Add p-value stats to ggplots
## Install SignatureEstimation package from source (from https://rdrr.io/github/WuyangFF95/SynSigRun/src/R/RunSignatureEstimationAttributionOnly.R)
## remotes::install_url("https://www.ncbi.nlm.nih.gov/CBBresearch/Przytycka/software/signatureestimation/SignatureEstimation.tar.gz")
library(SignatureEstimation) ## Use quadratic programing approach to get mutational signatures


## Loading external scripts ##
## Load script with general project settings and background functionality
source("SCRIPTS/setup.R")
## Load script to calculate the platinum enrichment score
source("SCRIPTS/platinumEnrichmentScore.R")

## Environment ##
dir.input <- "RESULTS/PyCloneVI/PyCloneResults/"
dir.out <- "RESULTS/PyCloneVI/CloneCharacterisation/"
.checkAndCreateDir(dir.out)
## PyClone results files
files <- list.files(dir.input, pattern = "_filtered[.]csv$", full.names = TRUE)
## Load annotation data
cohort <- general.load.data(paste0(dir.data, "Annotation/Cohort.csv")) ## Contains sample IDs and clinical information for each patient
## Patients with pre and post treatment samples
patients.prePost <- cohort$Patient[cohort$Post_Samples > 0]
## Mutation signatures
mutationSignatures.file <- "DATA/COSMIC_v3.3.1_SBS_GRCh37.txt"
## Known OAC signatures - COSMIC version 2 (from Noorani et al. 2017; Newell et al. 2019; Naeini et al. 2023)
oac_signatures <- paste0("SBS", c(1, 2, 3, 5, 8, 9, 13, "17a", "17b", 18, 20, 31, 35))
## Same signatures as used for clonal/unique mutation analysis
cohortSigs <- paste0("SBS", c(2, 13, 3, 36, 1, 5, 18, 20, 26, 44, 8, 39, 40, 93, "17a", "17b"))
cut <- 0.05 ## cut-off to select signatures for reassignment step


.preparePyCloneData <- function(file){
  clones <- general.load.data(file, silent = TRUE)
  clones$mutation_id <- gsub("_rRNA", "-rRNA", gsub("_A.4_", "-A.4_", gsub("_A.3_", "-A.3_", gsub("_Z15", "-Z15", gsub("_snr52", "-snr52", gsub("_19_", "-19_", gsub("__", "-", gsub("POST_", "", gsub("_RNA", "-RNA", gsub("_SRP", "-SRP", clones$mutation_id))))))))))
  
  clones$Patient <- gsub("_POST", "", gsub("_T[ULM]", "", clones$sample_id))
  ## Take each mutation only once
  clones <- unique(clones[, c("Patient", "mutation_id", "cluster_id", "cluster")])
  
  ## Split mutation ID in individual parts
  tmp <- strsplit(clones$mutation_id, "_")
  clones$Hugo_Symbol <- sapply(tmp, "[", 1)
  clones$Chromosome <- sapply(tmp, "[", 2)
  clones$Start_Position <- as.numeric(sapply(tmp, "[", 3))
  clones$End_Position <- as.numeric(sapply(tmp, "[", 4))
  clones$Reference_Allele <- sapply(tmp, "[", 5)
  clones$Tumor_Seq_Allele2 <- sapply(tmp, "[", 6)
  clones$Cluster_ID <- sapply(strsplit(clones$cluster, " "), "[", 1)
  
  ## Keep only chr 1-22
  clones <- clones[clones$Chromosome %in% paste0("chr", 1:22),]
  if(!all(grepl("chr", clones$Chromosome))){
    print(table(clones$Chromosome))
    cl <<- clones
    stop()
  }
  
  ## Remove inserts/deletions
  clones <- clones[which(nchar(clones$Reference_Allele) == 1 & clones$Reference_Allele != "-" & nchar(clones$Tumor_Seq_Allele2) == 1 & clones$Tumor_Seq_Allele2 != "-"),]
  return(clones)
}

clones.mutSignatures.v3.3 <- function(signa = cohortSigs){
  dir.out <- paste0(dir.out, "/MutationalSignatures_v3.3_", cut, "/")
  .checkAndCreateDir(dir.out)
  signa <- cohortSigs
  
  mutSigs <- general.load.data(mutationSignatures.file, row.names = 1)
  ## Remove sequencing artefect signatures
  mutSigs <- mutSigs[, !names(mutSigs) %in% signatures.seqArtefacts]
  
  signatures.unfiltered.all <- signatures.SigReassigned.all <- substitutions.all <- mutCounts.all <- data.frame()
  
  ## Assign mutational signatures to each clone per patient
  i <- 1
  total <- length(files)
  a <- sapply(files, function(file){
    ## Patient
    patient <- strsplit(basename(file), "_")[[1]][1]
    if(grepl("2samp", file)){
      patient <- paste0(patient, "_preOnly_2samp")
    }else if(grepl("preOnly", file)){
      patient <- paste0(patient, "_preOnly")
    }
    cat(patient)
    if(patient %in% patients.prePost){
      signa <- unique(c(signa, "SBS31", "SBS35"))
    }
    
    ## Store the assigned mutation signatures and single nucleotid substitutions per clone
    signatures.unfiltered <- signatures.SigReassigned <- substitutions <- data.frame()
    ## Load and prepare PyClone result
    clones <- .preparePyCloneData(file)
    ## Remove filtered out clones (Cluster_ID = NA)
    clones <- clones[!is.na(clones$Cluster_ID),]
    clones$Cluster_ID <- paste0("C", clones$Cluster_ID)
    
    print(table(clones$Cluster_ID))
    
    ## Get the single nucleotid substitutions based on the mutation position (default is hg19 assembly)
    clones.sig <- mut.to.sigs.input(mut.ref = clones, sample.id = "Cluster_ID", chr = "Chromosome", pos = "Start_Position", ref = "Reference_Allele", alt = "Tumor_Seq_Allele2")
#    cat("\n", as.character(Sys.time()), "-", i, "/", total, "|", patient)
    i <<- i + 1

    ## Save substitutions per clone
    substitutions <- clones.sig
    
    ## Get right order
    mutSigs <- mutSigs[names(clones.sig),]
    
    a <- sapply(unique(clones$Cluster_ID), function(clon){
      cat(" | clone", clon)
      ## Assign the mutational signatures # signatures.ref = mutSigs or signatures.nature2013
      tmp <- findSigExposures(M = t(clones.sig[clon,]), P = mutSigs[, signa], decomposition.method = decomposeQP)$exposures
      if(!"SBS31" %in% signa){
        te <- data.frame(clon = rep(0, 2), row.names = c("SBS31", "SBS35"))
        names(te) <- clon
        tmp <- rbind(tmp, te)
      }
      if(nrow(signatures.unfiltered) == 0){
        signatures.unfiltered <<- tmp
      }else{
        signatures.unfiltered <<- cbind(signatures.unfiltered, tmp)
      }
      
      ## To avoid over-fitting, signatures < 10% for each sample were omitted and mutations were reassigned to the remaining signatures.
      selectedSigs <- row.names(tmp[tmp > cut, , drop = FALSE])
      
      if(length(selectedSigs) == 0){
        warning(paste0("No signatures for ", clon, ". Cohort signatures used."))
        selectedSigs <- signa
        tmp <- findSigExposures(M = t(clones.sig[clon,]), P = mutSigs[, selectedSigs], decomposition.method = decomposeQP)$exposures
        selectedSigs <- row.names(tmp[tmp > cut, , drop = FALSE])
      }
      if(length(selectedSigs) == 1){
        selectedSigs <- unique(c(selectedSigs, "SBS17a", "SBS17b"))
      }
      removedSigs <- signa[!signa %in% selectedSigs]
      tmp.SigReassigned <- findSigExposures(M = t(clones.sig[clon,]), P = mutSigs[, selectedSigs], decomposition.method = decomposeQP)$exposures
      tmp.SigReassigned <- rbind(tmp.SigReassigned, matrix(rep(0, length(removedSigs)), dimnames = list(removedSigs, clon)))
      tmp.SigReassigned <- tmp.SigReassigned[signa, , drop = FALSE]
      if(!"SBS31" %in% signa){
        te <- data.frame(clon = rep(0, 2), row.names = c("SBS31", "SBS35"))
        names(te) <- clon
        tmp.SigReassigned <- rbind(tmp.SigReassigned, te)
      }
      ## Save the mutational signature percentages 
      if(nrow(signatures.SigReassigned) == 0){
        signatures.SigReassigned <<- tmp.SigReassigned
      }else{
        signatures.SigReassigned <<- cbind(signatures.SigReassigned, tmp.SigReassigned)
      }
    })
    cat("\nSaving results ...\n")
    ## Save
    general.save.data(signatures.unfiltered, paste0(dir.out, patient, "_unfiltered_clones"), row.names = TRUE, silent = TRUE)
    general.save.data(signatures.SigReassigned, paste0(dir.out, patient, "_Reassigned_clones"), row.names = TRUE, silent = TRUE)
    
    ## Add to overall data frames
    signatures.unfiltered <- as.data.frame(signatures.unfiltered)
    signatures.SigReassigned <- as.data.frame(signatures.SigReassigned)
    names(signatures.unfiltered) <- paste(patient, names(signatures.unfiltered), sep = "_")
    names(signatures.SigReassigned) <- paste(patient, names(signatures.SigReassigned), sep = "_")

    if(nrow(signatures.unfiltered.all) == 0){
      signatures.unfiltered.all <<- signatures.unfiltered
      signatures.SigReassigned.all <<- signatures.SigReassigned
    }else{
      signatures.unfiltered.all <<- cbind(signatures.unfiltered.all, signatures.unfiltered)
      signatures.SigReassigned.all <<- cbind(signatures.SigReassigned.all, signatures.SigReassigned)
    }
    
  })
  cat("\n")
  
  general.save.data(signatures.unfiltered.all, paste0(dir.out, "all_unfiltered_clones"), row.names = TRUE, silent = TRUE)
  general.save.data(signatures.SigReassigned.all, paste0(dir.out, "all_Reassigned_clones"), row.names = TRUE, silent = TRUE)
  
  ## Convert data into long format to be plotted
  signatures.SigReassigned.all$Signature <- row.names(signatures.SigReassigned.all)
  samples <- names(signatures.SigReassigned.all)[names(signatures.SigReassigned.all) != "Signature"]
  sigs.plot <- as.data.frame(pivot_longer(signatures.SigReassigned.all, cols = samples, names_to = "Sample", values_to = "Percent"))
  sigs.plot$Patient <- sapply(strsplit(sigs.plot$Sample, "_"), "[", 1)
  split <- strsplit(sigs.plot$Sample, "_C")
  sigs.plot$Clone <- as.factor(as.numeric(sapply(split, "[", 2)))
  sigs.plot$Sample <- sapply(split, "[", 1)
  sigs.plot$Timepoint[sigs.plot$Sample %in% patients.prePost] <- "prePost"
  sigs.plot$Timepoint <- ifelse(is.na(sigs.plot$Timepoint), "preOnly", "prePost")
  sigs.plot <- sigs.plot[sigs.plot$Percent > 0,]
  if(any(c("SBS31", "SBS35") %in% sigs.plot$Signature)){
    sigs.plot$Signature <- factor(sigs.plot$Signature, levels = c(cohortSigs, "SBS31", "SBS35"))
  }else{
    sigs.plot$Signature <- factor(sigs.plot$Signature, levels = cohortSigs)
  }

  ## Visualise
  ggplot(sigs.plot, aes(x = Clone, y = Percent, fill = Signature)) + 
    geom_bar(stat = "identity", color = "white") +
    facet_wrap(. ~ Sample, scale = "free") +
    theme_minimal() +
    labs(x = "", y = "Proportion [%]") +
    theme_test() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(paste0(dir.out, "Barcharts_allPatients_clones.pdf"), width = 10, height = 15)
  ggsave(paste0(dir.out, "Barcharts_allPatients_clones.png"), width = 10, height = 15, type = "cairo")
  
  
}


clones.calculateEnrichment <- function(){
  res.e <- data.frame()
  sapply(files[grepl(paste(patients.prePost, collapse = "|"), files) & !grepl("preOnly", files) & !grepl("old", files)], function(file){
    cat("Calculating platinum enrichment for", basename(file), "\n")
    ## Load and prepare PyClone result
    clones <- .preparePyCloneData(file)

    ## Add mutation details that are later needed
    clones$mutPattern <- paste(clones$Reference_Allele, clones$Tumor_Seq_Allele2, sep = ">")
    clones <- clones[clones$mutPattern %in% c("C>A", "G>T"),]
    
    ## Find sequence around all mutations
    data <- unique(clones[!is.na(clones$cluster), !names(clones) %in% c("sample_id", "cellular_prevalence", "cellular_prevalence_std")])
    
    i <<- 1
    data$seq <- apply(data, 1, .findSequence, seqLength, nrow(data))
    cat("\n")
    
    ## Add the previous and next base to the data frame
    data$prevBase <- substr(data$seq, start = seqLength, stop = seqLength)
    data$nextBase <- substr(data$seq, start = seqLength+2, stop = seqLength+2)
    
    ## Establish mutation context using the previous and next base
    data$contextPrev <- paste0(data$prevBase, data$Reference_Allele)
    data$contextNext <- paste0(data$Reference_Allele, data$nextBase)
    
    ## Split data up per clone
    data$clusterID_new <- paste0("C", data$Cluster_ID)
    clusters <- unique(data$clusterID_new)
    dataList <- list()
    a <- sapply(clusters, function(clone){
      dataList[[clone]] <<- data[data$clusterID_new == clone,]
    })
    
    ## Calculate enrichment scores
    es <- sapply(dataList, .calculateEnrichmentScore)
    ## Transform into data frame
    df <- data.frame(Clone = names(es), E.score = as.numeric(es), Patient = clones$Patient[1], Mutations = sapply(dataList, nrow))
    df <- merge(df, unique(data[, c("clusterID_new", "cluster")]), by.x = "Clone", by.y = "clusterID_new")
 
    ## Save enrichment scores
    if(nrow(res.e) == 0){
      res.e <<- df
    }else{
      res.e <<- rbind(res.e, df)
    }
    
    ## Plot values per clone
    ggplot(unique(df[, c("Clone", "E.score")]), aes(x = Clone, y = E.score)) +
      geom_point() +
      theme_bw()
    ggsave(paste0(dir.out, "/Clones_PlatinumEnrichment_", df$Patient[1], ".png"), width = 5, height = 5, dpi = 300, type = "cairo")
    ggsave(paste0(dir.out, "/Clones_PlatinumEnrichment_", df$Patient[1], ".pdf"), width = 5, height = 5, dpi = 300)
  })
  
  res.e$Clone <- as.factor(as.numeric(gsub("C", "", res.e$Clone)))
  ## Plot all at once (facet per patient)
  ggplot(unique(res.e[, c("Clone", "E.score", "Patient")]), aes(x = Clone, y = E.score)) +
    geom_point() +
    scale_y_continuous(limits = c(0.6, 2.4), breaks = c(0.6, 1.5, 2.4)) +
    facet_wrap(~ Patient, scales = "free_x") + 
    theme_bw() 
  ggsave(paste0(dir.out, "/Clones_PlatinumEnrichment_allPatients.png"), width = 7, height = 5, dpi = 300, type = "cairo")
  ggsave(paste0(dir.out, "/Clones_PlatinumEnrichment_allPatients.pdf"), width = 7, height = 5, dpi = 300)
  
  ## Save scores
  general.save.data(res.e, paste0(dir.out, "AllScores"))
}


clones.clinicalCharacteristics <- function(){
  cohort <- general.load.data(paste0(dir.data, "Annotation/Cohort.csv")) ## Contains sample IDs and clinical information for each patient
  cohort.files <- .prepareCohortFiles(paste0(dir.data, "Annotation/Cohort_sampleDetails.csv")) ## Contains sample information and file storage locations
  cohort.clones <- general.load.data(paste0(dir.data, "Annotation/CloneNumbers.csv"))
  dir.out <- "RESULTS/Cohort/ClonesVsClinical/"
  .checkAndCreateDir(dir.out)
  
  ## Add clone numbers to cohort
  cohort <- merge(cohort, cohort.clones)
  
  ## Group patients based on mean clone number
  av.clones <- median(cohort$Clones.pre)
  cohort$Clones.pre.group <- ifelse(cohort$Clones.pre > av.clones, yes = "high", no = "low")
  cohort$Clones.pre.group[cohort$Clones.pre == av.clones] <- "average"
  cohort$Clones.pre.group <- factor(cohort$Clones.pre.group, levels = c("low", "average", "high"))
  cohort$Clones.pre.group2 <- gsub("average", "low", cohort$Clones.pre.group)
  my_comparisons <- list( c("low", "average"), c("average", "high"), c("low", "high") )
  
  cohort$ypTNM.8th.AJCC.simp <- gsub("[ABC]", "", cohort$ypTNM.8th.AJCC)
  
  ## Survival vs clone groups
  coxph.p <- signif(.extractPvalue(coxph(Surv(cohort$OS.endo, cohort$OS) ~ as.numeric(cohort$Clones.pre)), mode = "coxph"), digits = 2)
  fit <- survfit(Surv(OS.endo, OS) ~ Clones.pre.group, data = cohort)
  ggsurv <- ggsurvplot(fit, data = cohort, 
                       break.time.by = 12,
                       pval = TRUE, 
                       risk.table = TRUE)
  ggsurv$plot <- ggsurv$plot + labs(title = "OS vs clone numbers (median 5)")
  png(paste0(dir.out, "OS.endo_clones_3groups.png"), type = "cairo", width = 5*200, height = 5*200, res = 200)
  print(ggsurv)
  dev.off()
  pdf(paste0(dir.out, "OS.endo_clones_3groups.pdf"), width = 7, height = 5)
  print(ggsurv)
  dev.off()
  
  fit <- survfit(Surv(OS.endo, OS) ~ Clones.pre.group2, data = cohort)
  ggsurv <- ggsurvplot(fit, data = cohort, 
                       break.time.by = 12,
                       pval = TRUE, 
                       risk.table = TRUE)
  ggsurv$plot <- ggsurv$plot + labs(title = "OS vs clone numbers (median 5)")
  png(paste0(dir.out, "OS.endo_clones_2groups.png"), type = "cairo", width = 5*200, height = 5*200, res = 200)
  print(ggsurv)
  dev.off()
  pdf(paste0(dir.out, "OS.endo_clones_2groups.pdf"), width = 7, height = 5)
  print(ggsurv)
  dev.off()
  cohort$Clones.pre.group2 <- factor(cohort$Clones.pre.group2, levels = c("high", "low"))
  cox <- coxph(Surv(OS.endo, OS) ~ Clones.pre.group2 + cTNM.8th.AJCC, data = cohort[cohort$cTNM.8th.AJCC %in% c("IIB", "III"),])
  png(paste0(dir.out, "OS.endo_clones_2groups_HR.png"), type = "cairo", width = 5*200, height = 2.5*200, res = 200)
  print(ggforest(cox))
  dev.off()
  pdf(paste0(dir.out, "OS.endo_clones_2groups_HR.pdf"), width = 5, height = 2.5)
  print(ggforest(cox))
  dev.off()
  cox <- coxph(Surv(OS.endo, OS) ~ Clones.pre.group2 + cTNM.8th.AJCC + Allocated.treatment, data = cohort[cohort$cTNM.8th.AJCC %in% c("IIB", "III") & cohort$Treatment != "Direct to Surgery" & cohort$Allocated.treatment != "Carboplatin & Paclitaxel & RT",])
  png(paste0(dir.out, "OS.endo_clones_2groups_HR_treatment.png"), type = "cairo", width = 5*200, height = 2.5*200, res = 200)
  print(ggforest(cox))
  dev.off()
  pdf(paste0(dir.out, "OS.endo_clones_2groups_HR_treatment.pdf"), width = 5, height = 2.5)
  print(ggforest(cox))
  dev.off()
  
  ## Survival vs stage
  fit <- survfit(Surv(OS.endo, OS) ~ cTNM.8th.AJCC, data = cohort)
  ggsurv <- ggsurvplot(fit, data = cohort, 
                       break.time.by = 12,
                       pval = TRUE, 
                       risk.table = TRUE)
  ggsurv$plot <- ggsurv$plot + labs(title = "OS vs cTNM stage")
  png(paste0(dir.out, "OS.endo_cTNM.png"), type = "cairo", width = 5*200, height = 5*200, res = 200)
  print(ggsurv)
  dev.off()
  pdf(paste0(dir.out, "OS.endo_cTNM.pdf"), width = 7, height = 5)
  print(ggsurv)
  dev.off()
  cohort$cTNM.8th.AJCC.2groups <- gsub("^III$", "III/IV", gsub("^IVA$", "III/IV",gsub("^IIB$", "I/II", gsub("^IB$", "I/II", cohort$cTNM.8th.AJCC))))
  fit <- survfit(Surv(OS.endo, OS) ~ cTNM.8th.AJCC.2groups, data = cohort)
  ggsurv <- ggsurvplot(fit, data = cohort, 
                       break.time.by = 12,
                       pval = TRUE, 
                       risk.table = TRUE)
  ggsurv$plot <- ggsurv$plot + labs(title = "OS vs cTNM stage")
  png(paste0(dir.out, "OS.endo_cTNM_2groups.png"), type = "cairo", width = 5*200, height = 5*200, res = 200)
  print(ggsurv)
  dev.off()
  pdf(paste0(dir.out, "OS.endo_cTNM_2groups.pdf"), width = 7, height = 5)
  print(ggsurv)
  dev.off()
  cox <- coxph(Surv(OS.endo, OS) ~ Clones.pre.group2 + cTNM.8th.AJCC.2groups, data = cohort)
  png(paste0(dir.out, "OS.endo_clones_2groups_stage_HR.png"), type = "cairo", width = 5*200, height = 2.5*200, res = 200)
  print(ggforest(cox))
  dev.off()
  pdf(paste0(dir.out, "OS.endo_clones_2groups_stage_HR.pdf"), width = 5, height = 2.5)
  print(ggforest(cox))
  dev.off()
  
  
  
  ## PFS
  coxph.p <- signif(.extractPvalue(coxph(Surv(cohort$PFS.endo, cohort$PFS) ~ as.numeric(cohort$Clones.pre)), mode = "coxph"), digits = 2)
  fit <- survfit(Surv(PFS.endo, PFS) ~ Clones.pre.group, data = cohort)
  ggsurvplot(fit, data = cohort,
             break.time.by = 12,
             risk.table = FALSE,
             conf.int = FALSE, pval = TRUE,
             ggtheme = theme_minimal(),
             font.main = 20, font.x = 18, font.y = 18, font.tickslab = 16, font.legend = 14,
             title = "PFS vs clone numbers (median 5)"
  )
  ggsave(paste0(dir.out, "PFS.endo_clones_3groups.pdf"), width = 7, height = 5)
  ggsave(paste0(dir.out, "PFS.endo_clones_3groups.png"), width = 7, height = 5, type = "cairo")

  
  ## Plot clone number vs clinical characteristics
  clinical.num <- c("Age", "endo.Tumour.Length.cm", "path.Tumour.Size.mm", "Total.nodes", "Positive.Nodes", "Positive.node.ratio", "Mutations.pre", "Neoantigens.preOnly", "DAI.pre.max", "DAI.pre.mean", "DAI.pre.median")
  clinical.cat <- c("Treatment", "PET.response", "cTNM.8th.AJCC", "Mandard.Score", "Mandard.simple", "Any.positive.nodes", "ypTNM.8th.AJCC", "ypTNM.8th.AJCC.simp", "Current.Status", "Smoker")

  ## Go through numerical clinical characteristics and plot them + stats vs clone numbers
  a <- sapply(clinical.num, function(clin){
    cat("Plotting", clin, "vs clones\n")
    cohort$clin <- cohort[, clin]
    
    sp <- ggscatter(cohort, x = "Clones.pre", y = "clin",
                    add = "reg.line",  # Add regression line
                    add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                    conf.int = TRUE, # Add confidence interval
                    ylab = clin,
                    xlab = "Clones"
    )
    # Add correlation coefficient
    sp + stat_cor(method = "pearson")
    ggsave(paste0(dir.out, clin, "_numClones.png"), width = 5, height = 5, type = "cairo")
    ggsave(paste0(dir.out, clin, "_numClones.pdf"), width = 5, height = 5)

    ## Group patients based on mean clone number
    ggplot(cohort, aes(x = Clones.pre.group, y = clin)) +
      geom_violin() +
      stat_summary(fun = "mean", geom = "crossbar", width = 0.5, colour = "black") +
      geom_beeswarm() +        
      stat_compare_means(method = "t.test", comparisons = my_comparisons) +
      stat_compare_means() +
      ylab(clin) +
      xlab("Clones") +
      theme_minimal()
    ggsave(paste0(dir.out, clin, "_numClones_group.png"), width = 5, height = 5, type = "cairo")
    ggsave(paste0(dir.out, clin, "_numClones_group.pdf"), width = 5, height = 5)
    
    ggplot(cohort, aes(x = Clones.pre.group2, y = clin)) +
      geom_violin() +
      stat_summary(fun = "mean", geom = "crossbar", width = 0.5, colour = "black") +
      geom_beeswarm() +        
      stat_compare_means(method = "t.test") +
      ylab(clin) +
      xlab("Clones") +
      theme_minimal()
    ggsave(paste0(dir.out, clin, "_numClones_2groups.png"), width = 5, height = 5, type = "cairo")
    ggsave(paste0(dir.out, clin, "_numClones_2groups.pdf"), width = 5, height = 5)
    
  })
  
  ## Go through categorical clinical characteristics and plot them + stats vs clone numbers
  a <- sapply(clinical.cat, function(clin){
    cat("Plotting", clin, "vs clones\n")
    cohort$clin <- cohort[, clin]
    
    ## Continuous clone numbers
    ggplot(cohort, aes(x = clin, y = Clones.pre)) +
      geom_violin() +
      stat_summary(fun = "mean", geom = "crossbar", width = 0.5, colour = "black") +
      geom_beeswarm() +        
      stat_compare_means() +
      xlab(clin) +
      ylab("Clones") +
      theme_minimal()
    ggsave(paste0(dir.out, clin, "_numClones.png"), width = 5, height = 5, type = "cairo")
    ggsave(paste0(dir.out, clin, "_numClones.pdf"), width = 5, height = 5)
  })
}




