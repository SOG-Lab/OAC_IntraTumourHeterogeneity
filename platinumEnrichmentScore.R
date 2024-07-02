cat("Project ID: OAC_Heterogeneity_MR_MT\n")
## R script to calculate the platinum enrichment scores ##
cat("Loading: R script to calculate the platinum enrichment scores\n")

## Loading external librarys ##
library(stringr) ## String manipulation
library(BSgenome.Hsapiens.UCSC.hg19) ## Human genome to extract sequence around mutations
library(ggpubr) ## Add p-values to plot
library(ggbeeswarm) ## Nicer jitter plotting

## Loading external scripts ##
## Load script with general project settings and background functionality
source("SCRIPTS/setup.R")

## Environment ##
seqLength <- 20
samplePvalues <- c()
out.dir <- "RESULTS/PlatinumEnrichmentScore/"
.checkAndCreateDir(out.dir)
vcf.cols <- c("#CHROM", "POS", "REF", "ALT")

## Internally called functions ##

.findSequence <- function(row, length, total){
  cat("\r", i, "/", total)
  i <<- i + 1
  pos <- as.numeric(row["Start_Position"])
  seq <- as.character(BSgenome.Hsapiens.UCSC.hg19[[row["Chromosome"]]][(pos-length):(pos+length)])
  return(seq)
}

.gsubRecursive <- function(seqs, pattern){
  seqs <- unlist(lapply(seqs, function(s){
    while(str_count(s, pattern) > 0 && !is.na(str_count(s, pattern))){
      s <- gsub(pattern, substr(pattern, start = 1, stop = 1), s)
    }
    return(s)
  }))
  return(seqs)
}

.calculateEnrichmentScore <- function(data, mode = NA){
  ## Number of mutated C (G)
  mutationsC <- nrow(data)
  
  ## Split data in C>A mutations and G>T mutations
  data.CA <- data[data$mutPattern == "C>A",]
  data.GT <- data[data$mutPattern == "G>T",]
  
  ## Number of mutated C (G) falling in a CpC (GpG) dinucleotide
  mutationsCpC <- sum(length(which(data$contextPrev == "CC" | data$contextNext == "CC")), length(which(data$contextPrev == "GG" | data$contextNext == "GG")))
  ## Number of C (G) within the 41-base region centered on the mutated C (G)
  contextC <- sum(str_count(data.CA$seq, "C")) + sum(str_count(data.GT$seq, "G"))
  ## Number of CpC (GpG) dinucleotides in 41-base region centered on the mutated C (G)
  contextCpC <- sum(nchar(data.CA$seq) - nchar(.gsubRecursive(data.CA$seq, "CC"))) + sum(nchar(data.GT$seq) - nchar(.gsubRecursive(data.GT$seq, "GG")))
  
  e <- (mutationsCpC*contextC)/(mutationsC*contextCpC)
  
  ## determine enrichment
  if(!is.na(mode)){
    tmp <- matrix(c(mutationsCpC, mutationsC, contextCpC, contextC-contextCpC), nrow = 2, byrow = TRUE)
    samplePvalues <<- c(samplePvalues, chisq.test(tmp)$p.value)
  }
  return(e)
}

.prepareData <- function(data, patient, perSample = FALSE){
  ## Remove inserts/deletions
  dataAll <- data[which(nchar(data$Reference_Allele) == 1 & data$Reference_Allele != "-" & nchar(data$Tumor_Seq_Allele2) == 1 & data$Tumor_Seq_Allele2 != "-"),]
  ## Keep mutations found only in pre- or post-treatment samples (no clonal/subclonal mutations)
  dataAll <- dataAll[which(dataAll$Classification.timepoint %in% c("unique.post", "unique.pre")),]
  
  ## Add mutation details that are later needed
  dataAll$mutPattern <- paste(dataAll$Reference_Allele, dataAll$Tumor_Seq_Allele2, sep = ">")
  dataAll <- dataAll[dataAll$mutPattern %in% c("C>A", "G>T"),]
  d <<- dataAll
  ## Find sequence around all mutations
  tmp <- unique(dataAll[, c("mutationID", "Chromosome", "Start_Position")])
  i <<- 1
  tmp$seq <- apply(tmp, 1, .findSequence, seqLength, nrow(tmp))
  ## Merge sequence information with mutation data frame
  dataAll <- merge(dataAll, tmp)
  dataAll$prevBase <- substr(dataAll$seq, start = seqLength, stop = seqLength)
  dataAll$nextBase <- substr(dataAll$seq, start = seqLength+2, stop = seqLength+2)
  
  ## Establish mutation context using the previous and next base
  dataAll$contextPrev <- paste0(dataAll$prevBase, dataAll$Reference_Allele)
  dataAll$contextNext <- paste0(dataAll$Reference_Allele, dataAll$nextBase)
  
  ## Split data up per sample
  if(perSample){
    samples <- unique(data$Sample_ID)
    dataList <- list()
    a <- sapply(samples, function(samp){
      dataList[[samp]] <<- data[data$Sample_ID == samp,]
    })
  }

  ## Add an entry for all pre-mutations (unique.pre, subclonal.pre, clonal.pre) and post-mutations (unique.post, clonal.post)
  dataList[[paste0(patient, "_pre")]] <- dataAll[dataAll$Classification.timepoint %in% c("unique.pre", "subclonal.pre", "clonal.pre"),]
  dataList[[paste0(patient, "_post")]] <- dataAll[dataAll$Classification.timepoint %in% c("unique.post", "clonal.post"),]
  d <<- dataAll
  return(dataList)
}

## Functions ##
calculatePlatinumEnrichmentScore <- function(){
  ## Load annotation data
  cohort <- general.load.data(paste0(dir.data, "Annotation/Cohort.csv")) ## Contains sample IDs and clinical information for each patient
  cohort.files <- general.load.data(paste0(dir.data, "Annotation/Cohort_sampleDetails.csv")) ## Contains sample information and file storage locations
  cohort.files <- cohort.files[!is.na(cohort.files$Patient),] ## Remove empty rows that sometimes get attached at the bottom when saving from excel
  cohort.files$Sample.type[cohort.files$Sample_ID == "JHH121_NA"] <- "NA" ## Sample type was set to <NA> during import
  ## Remove normal samples
  cohort.files <- cohort.files[!cohort.files$Sample.type %in% c("NA", "BC"), ]
  
  ## Patients with pre and post treatment samples
  patients <- cohort$Patient[cohort$Post_Samples > 0]

  cat("Prepare data ...\n")

  ## Load all mutations
  mutationsList <- general.load.data("DATA/allMutations_collection_prePost.RData")

  j <- 1
  total <- length(patients)
  z <- lapply(patients, function(patient){
    cat(j, "/", total, "|", patient, "... \n")
    j <<- j + 1  
    data <- mutationsList[[patient]]$allMutations
    
    ## Prepare data and split into sample specific list entries
    dataList <- .prepareData(data, patient)
    
    cat("\n")
    ## Calculate enrichment scores
    es <- sapply(dataList, .calculateEnrichmentScore)

    return(list(E.score = es))
  })
  names(z) <- patients
  print(z)
  Escores <- as.data.frame(unlist(sapply(z, "[", 1)))
  names(Escores) <- "E.score"
  Escores <- Escores[!is.na(Escores$E.score), , drop = FALSE]
  row.names(Escores) <- Escores$Sample_ID <- sapply(strsplit(gsub("E.score", "", row.names(Escores)), "[..]"), "[", 3)
  Escores$Patient <- sapply(strsplit(Escores$Sample_ID, "_"), "[", 1)
  Escores$timepoint <- "Pre"
  Escores$timepoint[grepl("POST", Escores$Sample_ID, ignore.case = TRUE)] <- "Post"
  Escores$timepoint <- factor(Escores$timepoint, levels = c("Pre", "Post"))
  general.save.data(Escores, paste0(out.dir, "Escores"))
  
  ## Wide format for easier comparison of pre/post
  tmp <- Escores[!grepl("_T", Escores$Sample_ID), names(Escores) != "Sample_ID"]
  tmp <- as.data.frame(pivot_wider(tmp, names_from = timepoint, values_from = E.score))
  tmp$diff <- tmp$Pre - tmp$Post
  general.save.data(tmp, paste0(out.dir, "Escores_prePost"))



  ## Plot the scores from the sample unique mutations (1 dot per sample)
  data.plot.samp <- Escores[grepl("_T", Escores$Sample_ID),]
  
  if(nrow(data.plot.samp) > 0){
    ggplot(data.plot.samp, aes(x = timepoint, y = E.score, colour = timepoint)) + 
      geom_violin() + 
      geom_jitter() +
      theme_minimal() +
      labs(x = "", y = "Platinum signature C>A (CpC context) enrichment")
    
    ggboxplot(data.plot.samp, x = "timepoint", y = "E.score", color = "timepoint", palette = "jco", add = "jitter", title = "Unique mutations per sample",
              ylab = "Platinum signature C>A (CpC context) enrichment", xlab = "", legend = "right") +
      stat_compare_means(method = "t.test") 
    ggsave(paste0(out.dir, "/PlatinumEnrichment_unique.png"), width = 5, height = 5, dpi = 300, type = "cairo")
    ggsave(paste0(out.dir, "/PlatinumEnrichment_unique.pdf"), width = 5, height = 5, dpi = 300)
    
    data.plot.samp$timepoint <- factor(data.plot.samp$timepoint, levels = c("Pre", "Post"))
    ggplot(data.plot.samp, aes(x = timepoint, y = E.score)) +
      geom_violin() +
      stat_summary(fun = "mean", geom = "crossbar", color = "black") + 
      geom_beeswarm() +
      stat_compare_means(method = "t.test") +
      ylab("Platinum signature C>A (CpC context) enrichment") +
      xlab("")
    ggsave(paste0(out.dir, "/PlatinumEnrichment_unique_violin.png"), width = 5, height = 5, dpi = 300, type = "cairo")
    ggsave(paste0(out.dir, "/PlatinumEnrichment_unique_violin.pdf"), width = 5, height = 5, dpi = 300)
    
    ## Plot the scores from the sample unique mutations (1 dot per sample) - split per patient
    ggplot(data.plot.samp, aes(x = timepoint, y = E.score, color = timepoint))+
      geom_point(size = 4) +
      geom_point(data = data.plot, mapping = aes(x = timepoint, y = E.score), color = "black", fill = NA) + 
      labs(title = "Unique mutations", y = "Platinum signature C>A (CpC context) enrichment", x = "") +
      facet_grid(. ~ Patient) + 
      theme_linedraw()
    ggsave(paste0(out.dir, "/PlatinumEnrichment_unique_patient.png"), width = 10, height = 5, dpi = 300, type = "cairo")
    ggsave(paste0(out.dir, "/PlatinumEnrichment_unique_patient.pdf"), width = 10, height = 5, dpi = 300)
  }

  ## Plot the unique pre vs unique post mutations (1 dot per patient)
  data.plot <- Escores[!grepl("_T", Escores$Sample_ID),]
  
  ggboxplot(data.plot, x = "timepoint", y = "E.score", color = "timepoint", palette = "jco", add = "jitter", title = "Unique mutations per timepoint",
            ylab = "Platinum signature C>A (CpC context) enrichment", xlab = "", legend = "right") +
    geom_line(aes(group = Patient)) + 
    stat_compare_means(method = "t.test", paired = TRUE) 
  ggsave(paste0(out.dir, "/PlatinumEnrichment_prePost.png"), width = 5, height = 5, dpi = 300, type = "cairo")
  ggsave(paste0(out.dir, "/PlatinumEnrichment_prePost.pdf"), width = 5, height = 5, dpi = 300)
  
  ggplot(data.plot, aes(x = timepoint, y = E.score)) +
    geom_violin() + 
    stat_summary(fun = "mean", geom = "crossbar", color = "black") + 
    geom_beeswarm() +
    stat_compare_means(method = "t.test", paired = TRUE)
  ggsave(paste0(out.dir, "/PlatinumEnrichment_prePost_violin.png"), width = 5, height = 5, dpi = 300, type = "cairo")
  ggsave(paste0(out.dir, "/PlatinumEnrichment_prePost_violin.pdf"), width = 5, height = 5, dpi = 300)

}

