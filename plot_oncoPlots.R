cat("*** Project ID: OAC_Heterogeneity_MR_MT ***\n")
## R script to plot oncoplots ##
cat("Loading: R script to plot oncoplots\n")

## Loading external librarys ##
library(maftools) ## oncoplot()
library(biomaRt) ## Convert ENSEMBL transcript IDs to gene IDs
library(ggstatsplot) ## For ggbarstats()

## Loading external scripts ##
## Load script with general project settings and background functionality
source("SCRIPTS/setup.R")
## Load additional functionality for the analysis of mutations
source("SCRIPTS/mutationalAnalysis.R") 

## Environment ##
out.dir <- "RESULTS/Oncoplot/"
.checkAndCreateDir(out.dir)
cancerGenes.file <- "DATA/COSMIC_CancerGeneCensus_07Jun2023.csv"
cancerGenes <- general.load.data(cancerGenes.file)
names(cancerGenes) <- gsub("[.]Symbol", "", names(cancerGenes))

.prepareAnnotation <- function(file){
  cohort.files <- general.load.data(file) ## Contains sample information and file storage locations
  cohort.files <- cohort.files[!is.na(cohort.files$Patient),] ## Remove empty rows that sometimes get attached at the bottom when saving from excel
  ## Remove normal (BC, NA) samples
  cohort.files <- cohort.files[!cohort.files$Sample.type %in% c("NA", "BC"),]
  cohort.files$Patient <- factor(cohort.files$Patient, levels = cohort$Patient[order(cohort$SurvivalSorting)])
  cohort.files <- cohort.files[order(cohort.files$Patient),]
  cohort.files$Tumor_Sample_Barcode <- factor(cohort.files$Sample_ID, levels = cohort.files$Sample_ID)
  return(cohort.files)
}

plotOnco_clonal <- function(){
  plottedMutationTypes <- c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del", "In_Frame_Del", "Splice_Site", "Nonstop_Mutation", "Translation_Start_Site")
  
  ## COSMIC gene set
  file.genes <- general.load.data(cancerGenes.file)
  names(file.genes) <- gsub("[.]Symbol", "", names(file.genes))
  set.name <- "COSMIC"
  
  ## Annotation set up
  cohort <- general.load.data(paste0(dir.data, "Annotation/Cohort.csv")) ## Contains sample IDs and clinical information for each patient
  cohort <- rbind(cohort, cohort)
  cohort$Tumor_Sample_Barcode <-  paste0(cohort$Patient, rep(c("_clonal", "_unique"), each = 29))
  cohort$Classification <- factor(ifelse(grepl("clonal", cohort$Tumor_Sample_Barcode), "clonal", "unique"), c("clonal", "unique"))
  cohort <- cohort[order(cohort$Number, cohort$Classification, decreasing = FALSE),]
  cohort$Tumor_Sample_Barcode <- factor(cohort$Tumor_Sample_Barcode, levels = cohort$Tumor_Sample_Barcode)
  cohort$Patient <- factor(cohort$Patient, levels = unique(cohort$Patient))

  ## Load data if available
  mutationCollection.file <- "DATA/allMutations_collection_preOnly.RData"
  res <- general.load.data(mutationCollection.file)
  
  allMutations <- data.frame()
  i <- 1
  total <- length(res)
  a <- sapply(names(res), function(patient){
    cat("\r", i, "/", total, "|", patient)
    i <<- i + 1
    tmp <- res[[patient]]$allMutations
    tmp$Tumor_Sample_Barcode <- ifelse(grepl("unique", tmp$Classification), paste0(patient, "_unique"), paste0(patient, "_clonal"))
    tmp$Patient <- patient
    tmp$Sample_ID <- NULL
    tmp <- unique(tmp)
    if(nrow(allMutations) == 0){
      allMutations <<- tmp
    }else{
      allMutations <<- rbind(allMutations, tmp)
    }
  })
  cat("\n")

  ## Fix Variant_Classification labels so that read.maf recognises them
  allMutations$Variant_Classification <- gsub("DEL", "Del", allMutations$Variant_Classification) 
  allMutations$Variant_Classification <- gsub("INS", "Ins", allMutations$Variant_Classification) 
  
  ## Remove genes with '.' in name -> pseudo genes
  allMutations <- allMutations[!grepl("[.]", allMutations$Hugo_Symbol),]
  
  ## Load gene set
  genes <- general.load.data(file.genes)

  ## Keep only genes that are in the gene list
  dataCOSMIC <- data <- allMutations[allMutations$Hugo_Symbol %in% genes$Gene,]
  
  ## Read as maf object
  maf <- read.maf(data, cohort)
  
  png(paste0(out.dir, "Oncoplot_", set.name, ".png"), width = 15*150, height = 20*150, res = 150, type = "cairo")
  oncoplot(maf = maf, top = nrow(maf@gene.summary[maf@gene.summary$MutatedSamples >= 2,]), showTumorSampleBarcodes = TRUE, removeNonMutated = FALSE, annotationDat = cohort, annotationOrder = levels(cohort$Patient), sortByAnnotation = TRUE, clinicalFeatures = c("Patient", "Classification"), barcodeSrt = 45, gene_mar = 15, sampleOrder = cohort$Tumor_Sample_Barcode, fontSize = 2, annotationFontSize = 2, legendFontSize = 2)
  dev.off()
  pdf(paste0(out.dir, "Oncoplot_", set.name, ".pdf"), width = 30, height = 25)
  oncoplot(maf = maf, top = nrow(maf@gene.summary[maf@gene.summary$MutatedSamples >= 2,]), showTumorSampleBarcodes = TRUE, removeNonMutated = FALSE, annotationDat = cohort, annotationOrder = levels(cohort$Patient), sortByAnnotation = TRUE, clinicalFeatures = c("Patient", "Classification"), barcodeSrt = 45, gene_mar = 15, sampleOrder = cohort$Tumor_Sample_Barcode, fontSize = 2, annotationFontSize = 2, legendFontSize = 2)
  dev.off()
  png(paste0(out.dir, "Oncoplot_min3samples_", set.name, ".png"), width = 15*150, height = 20*150, res = 150, type = "cairo")
  oncoplot(maf = maf, top = nrow(maf@gene.summary[maf@gene.summary$MutatedSamples >= 3,]), showTumorSampleBarcodes = TRUE, removeNonMutated = FALSE, annotationDat = cohort, annotationOrder = levels(cohort$Patient), sortByAnnotation = TRUE, clinicalFeatures = c("Patient", "Classification"), barcodeSrt = 45, gene_mar = 15, sampleOrder = cohort$Tumor_Sample_Barcode, fontSize = 2, annotationFontSize = 2, legendFontSize = 2)
  dev.off()
  pdf(paste0(out.dir, "Oncoplot_min3samples_", set.name, ".pdf"), width = 30, height = 25)
  oncoplot(maf = maf, top = nrow(maf@gene.summary[maf@gene.summary$MutatedSamples >= 3,]), showTumorSampleBarcodes = TRUE, removeNonMutated = FALSE, annotationDat = cohort, annotationOrder = levels(cohort$Patient), sortByAnnotation = TRUE, clinicalFeatures = c("Patient", "Classification"), barcodeSrt = 45, gene_mar = 15, sampleOrder = cohort$Tumor_Sample_Barcode, fontSize = 2, annotationFontSize = 2, legendFontSize = 2)
  dev.off()
  
  ## Stats clonal vs unique
  data <- allMutations[allMutations$Hugo_Symbol %in% maf@gene.summary$Hugo_Symbol[maf@gene.summary$MutatedSamples >= 3] & allMutations$Variant_Classification %in% plottedMutationTypes,]
  maf <- read.maf(data, cohort)
  png(paste0(out.dir, "Oncoplot_min3samples_", set.name, ".png"), width = 15*150, height = 20*150, res = 150, type = "cairo")
  oncoplot(maf = maf, top = nrow(maf@gene.summary), showTumorSampleBarcodes = TRUE, removeNonMutated = FALSE, annotationDat = cohort, annotationOrder = levels(cohort$Patient), sortByAnnotation = TRUE, clinicalFeatures = c("Patient", "Classification"), barcodeSrt = 45, gene_mar = 15, sampleOrder = cohort$Tumor_Sample_Barcode, fontSize = 2, annotationFontSize = 2, legendFontSize = 2)
  dev.off()
  pdf(paste0(out.dir, "Oncoplot_min3samples_", set.name, ".pdf"), width = 30, height = 25)
  oncoplot(maf = maf, top = nrow(maf@gene.summary), showTumorSampleBarcodes = TRUE, removeNonMutated = FALSE, annotationDat = cohort, annotationOrder = levels(cohort$Patient), sortByAnnotation = TRUE, clinicalFeatures = c("Patient", "Classification"), barcodeSrt = 45, gene_mar = 15, sampleOrder = cohort$Tumor_Sample_Barcode, fontSize = 2, annotationFontSize = 2, legendFontSize = 2)
  dev.off()
  
  dfCOSMIC <- df <- cohort[, c("Patient", "Tumor_Sample_Barcode", "Classification")]
  genes <- unique(data$Hugo_Symbol)
  df[, genes] <- 0
  sapply(genes, function(gen){
    df[df$Tumor_Sample_Barcode %in% data$Tumor_Sample_Barcode[data$Hugo_Symbol == gen], gen] <<- 1
    tab <- table(df[, gen], df$Classification)
    p <- NA
    try(
      p <- signif(fisher.test(tab)$p.value, 2)
    )
    df$Gen <- ifelse(df[, gen] == 0, "wt", "mt")
    ggbarstats(
      df, Gen, Classification,
      results.subtitle = FALSE,
      subtitle = paste0(
        "Fisher's exact test", ", p-value = ",
        ifelse(p < 0.001, "< 0.001", round(p, 3))
      ),
      legend.title = gen
    )
    ggsave(paste0(out.dir, "Barplot/Mutation_presentAbsent_", "_", gen, ".png"), width = 3, height = 3, type = "cairo")
    ggsave(paste0(out.dir, "Barplot/Mutation_presentAbsent_", "_", gen, ".pdf"), width = 3, height = 3)
  })
  df$COSMIC.number <- rowSums(df[, genes])
  df$COSMIC <- ifelse(df$COSMIC.number > 0, "mt", "wt")
  tab <- table(df$COSMIC, df$Classification)
  p <- NA
  try(
    p <- signif(fisher.test(tab)$p.value, 2)
  )
  ggbarstats(
    df, COSMIC, Classification,
    results.subtitle = FALSE,
    subtitle = paste0(
      "Fisher's exact test", ", p-value = ",
      ifelse(p < 0.001, "< 0.001", round(p, 3))
    ),
    legend.title = "COSMIC"
  )
  ggsave(paste0(out.dir, "Barplot/Mutation_presentAbsent_COSMIC.png"), width = 3, height = 3, type = "cairo")
  ggsave(paste0(out.dir, "Barplot/Mutation_presentAbsent_COSMIC.pdf"), width = 3, height = 3)
    
    
  ggplot(df, aes(x = Classification, y = COSMIC.number)) +
    geom_violin() +
    stat_summary(fun = "mean", geom = "crossbar", color = "black") +
    stat_compare_means(method = "t.test") +
    geom_beeswarm() +
    theme_minimal()
  ggsave(paste0(out.dir, "Barplot/Mutation_COSMICnumbers.png"), width = 3, height = 3, type = "cairo")
  ggsave(paste0(out.dir, "Barplot/Mutation_COSMICnumbers.pdf"), width = 3, height = 3)
  
  genes <- unique(dataCOSMIC$Hugo_Symbol)
  ## Remove synonymous variants (keep only mutations that would be plotted in an oncoplot)
  dataCOSMIC <- dataCOSMIC[dataCOSMIC$Variant_Classification %in% plottedMutationTypes,]
    
  dfCOSMIC[, genes] <- 0
  a <- sapply(genes, function(gen){
    dfCOSMIC[dfCOSMIC$Tumor_Sample_Barcode %in% dataCOSMIC$Tumor_Sample_Barcode[dataCOSMIC$Hugo_Symbol == gen], gen] <<- 1
  })
  muts.clonal <- colSums(dfCOSMIC[dfCOSMIC$Classification == "clonal", genes])
  muts.unique <- colSums(dfCOSMIC[dfCOSMIC$Classification == "unique", genes])
  ## Remove genes that don't have clonal or subclonal mutations
  muts.clonal <- muts.clonal[!names(muts.clonal) %in% exclude]
  muts.unique <- muts.unique[!names(muts.unique) %in% exclude]
    
  exclude <- intersect(names(muts.clonal)[muts.clonal == 0], names(muts.unique)[muts.unique == 0])
  df.plot.clonal <- data.frame(Classification = "clonal", Hugo_Symbol = names(muts.clonal), mutCount = as.numeric(muts.clonal))
  df.plot.unique <- data.frame(Classification = "unique", Hugo_Symbol = names(muts.unique), mutCount = as.numeric(muts.unique))
  df.plot <- rbind(df.plot.clonal, df.plot.unique)
  df.plot <- df.plot[!df.plot$Hugo_Symbol %in% exclude,]
  
  ## Violin plot number of clonal mutations in COSMIC gene vs number of unique mutations in COSMIC genes
  ggplot(df.plot, aes(x = Classification, y = mutCount)) +
    geom_violin() +
    stat_summary(fun = "mean", geom = "crossbar", color = "black") +
    stat_compare_means(method = "t.test", paired = TRUE) +
    theme_minimal()
    
  ggsave(paste0(out.dir, "Violin.png"), width = 3, height = 3, type = "cairo")
  ggsave(paste0(out.dir, "Violin.pdf"), width = 3, height = 3)
}

