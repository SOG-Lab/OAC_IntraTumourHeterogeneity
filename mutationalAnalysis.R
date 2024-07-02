cat("Project ID: OAC_Heterogeneity_MR_MT\n")
## R script to analyse the single nucleotide mutations ##
cat("Loading: R script to analyse the single nucleotide mutations\n")

## Loading external libraries ##
library(venn) ## Venn diagram plotting
library(ggbeeswarm) ## Nicer distribution of dots on ggplot

## Loading external scripts ##
## Load script with general project settings and background functionality
source("SCRIPTS/setup.R")


## Function to identify clonal and subclonal mutations per patient
.mutationClassification <- function(data){
  mutations <- list() ## List that will contain all mutation dataframes (1 dataframe per sample)
  allMutations <- data.frame() ## Dataframe that will contain all unique somatic mutations and will indicate if a mutation is clonal or subclonal
  mutationIDs <- list()
  patient <- as.character(data$Patient[1])
  ## Columns that are required for each mutation
  cols <- c("Hugo_Symbol", "mutationID", "Chromosome", "Start_Position", "End_Position", "Strand", "Variant_Classification", "Variant_Type", "Effect_Ontology", "Effect_Class", "Reference_Allele", "Tumor_Seq_Allele2")
  i <- 1
  total <- nrow(data)
  a <- apply(data, 1, function(row){
    sample <- row[["Sample_ID"]]
    cat("\r", i, "/", total, "|", sample, "- ")
    i <<- i + 1
    d <- general.load.data(paste0(row[["MAF.folder"]], row[["MAF.file"]]))
    
    ## Add mutation ID
    d$mutationID <- paste(d$Hugo_Symbol, d$Chromosome, d$Start_Position, d$End_Position, d$Reference_Allele, d$Tumor_Seq_Allele2, sep = "_")
    mutations[[sample]] <<- .filterVariants(d, dbSNP = NULL) ## Don't remove known SNPs
    mutationIDs[[sample]] <<- mutations[[sample]]$mutationID
    ## Add somatic mutations to allMutations dataframe
    if(!grepl("BC", sample) & !grepl("NA", sample)){
      tmp <- mutations[[sample]][, cols]
      tmp$Sample_ID <- sample
      if(nrow(allMutations) == 0){
        allMutations <<- tmp
      }else{
        allMutations <<- rbind(allMutations, tmp)
      }
    }    
  })
  
  ## Count number of occurrences for each unique mutation ID in the tumour samples
  cat("Count number of occurences for each unique mutation ID in the tumour samples \n")
  tmp <- aggregate(data.frame(sampleCount = allMutations$mutationID), list(mutationID = allMutations$mutationID), length)
  ## Merge the mutation counts with the mutation details
  tmp <- merge(tmp, unique(allMutations[, names(allMutations) != "Sample_ID"]))
  ## Count occurrences in the post treatment samples (if available)
  if(any(grepl("POST", allMutations$Sample_ID))){
    tmp.post <- aggregate(data.frame(sampleCount.post = allMutations$mutationID[grepl("POST", allMutations$Sample_ID)]), list(mutationID = allMutations$mutationID[grepl("POST", allMutations$Sample_ID)]), length)
    ## Count occurrences in the post treatment samples
    tmp.pre <- aggregate(data.frame(sampleCount.pre = allMutations$mutationID[!grepl("POST", allMutations$Sample_ID)]), list(mutationID = allMutations$mutationID[!grepl("POST", allMutations$Sample_ID)]), length)
    tmp <- merge(tmp, tmp.post[, c("mutationID", "sampleCount.post")], all.x = TRUE)
    tmp <- merge(tmp, tmp.pre[, c("mutationID", "sampleCount.pre")], all.x = TRUE)
  }
  
  ## Add mutation classification (clonal, subclonal, unique)
  ## Clonal: present in all samples from the patient
  ## Subclonal: present in at least 2 but not all samples from the patient
  ## Unique: present in only 1 sample from the patient
  tmp$Classification <- "subclonal"
  tmp$Classification[tmp$sampleCount == 1] <- "unique"
  tmp$Classification[tmp$sampleCount == length(unique(allMutations$Sample_ID))] <- "clonal"
  tmp$Classification.timepoint <- NA
  
  if(any(grepl("POST", allMutations$Sample_ID))){
    tmp$Classification.timepoint[tmp$sampleCount.post > 0 & is.na(tmp$sampleCount.pre)] <- "unique.post"
    tmp$Classification.timepoint[tmp$sampleCount.pre > 0 & is.na(tmp$sampleCount.post)] <- "unique.pre"
  }
  allMutations <- merge(allMutations, tmp[, c("mutationID", "Classification", "Classification.timepoint", "sampleCount")], all.x = TRUE)
  
  return(list(patient = as.character(patient), mutationDFs = mutations, allMutations = allMutations, allMutationsCounts = tmp, mutationIDs = mutationIDs))
}



prepareMutationCollections <- function(mode, insert = ""){
  
  ## Load annotation files and prepare them
  dir.out <- paste0("RESULTS/MutationVennDiagrams", insert, "/")
  cohort <- general.load.data(paste0(dir.data, "Annotation/Cohort.csv")) ## Contains sample IDs and clinical information for each patient
  cohort.files <- .prepareCohortFiles(paste0(dir.data, "Annotation/Cohort_sampleDetails.csv")) ## Contains sample information and file storage locations
  
  if(mode == "preOnly"){
    cohort.files <- cohort.files[!grepl("POST", cohort.files$Sample_ID) & !cohort.files$Sample_ID %in% exclude.pre,]
  }else if(mode == "prePost"){
    cohort <- cohort[cohort$Post_Samples > 0,]
    cohort.files <- cohort.files[cohort.files$Patient %in% cohort$Patient,]
  }else{
    stop("Please specify mode (preOnly/prePost). \n")
  }
  
  pdf(paste0(dir.out, "allVennDiagrams_", mode, insert, ".pdf"), width = 3, height = 3)
  
  ## Go through the patients and load the data for each of the associated samples
  i <- 1
  total <- nrow(cohort)
  res <- apply(cohort, 1, function(row){
    patient <- row[["Patient"]]
    cat("\n--- Patient", i, "/", total, "|", patient, "--- \n")
    i <<- i + 1
    cohort.files <- cohort.files[cohort.files$Patient == patient & !cohort.files$Sample.type %in% c("BC", "NA"),]
    cat("Samples:", cohort.files$Sample_ID, "\n")
    
    dataList <- .mutationClassification(cohort.files)
    dl <<- dataList
    
    ## Save list of mutation dataframes per patient
    general.save.data(dataList, paste0("DATA/MAFs/Collections/", patient, "_mutation_collection_", mode, insert), csv = FALSE)
    ## Plot Venn diagram
    cat("Plot Venn diagram ...\n")
    .plotVenn(dataList$mutationIDs, paste0(dir.out, patient, "_", mode, insert, ".pdf"))
    return(dataList)
  })
  dev.off()
  names(res) <- cohort$Patient
  
  ## Save data
  general.save.data(res, paste0("DATA/allMutations_collection_", mode, insert), csv = FALSE)
}

plotViolinClonalUnique <- function(mode, insert = ""){
  if(mode == "preOnly"){
    mutationCollection.file <- "DATA/allMutations_collection_preOnly.RData"
  }else if(mode == "prePost"){
    mutationCollection.file <- "DATA/allMutations_collection_prePost.RData"
  }
  
  df <- data.frame()
  
  i <- 1
  total <- length(res)
  ## Go through mutation dataframes to create pie charts for percent clonal/unique mutations per patient 
  a <- sapply(names(res), function(patient){
    cat("\r", i, "/", total, "|", patient, "        ")
    i <<- i + 1
    data <- res[[patient]]$allMutations
    tab <- table(data$Classification)
    
    tmp <- data.frame(Patient = patient, "clonal" = as.numeric(tab["clonal"])/2, "unique" = as.numeric(tab["unique"]))
    
    if(nrow(df) == 0){
      df <<- tmp
    }else{
      df <<- rbind(df, tmp)
    }   
  })
  cat("\n")
  df$total <- df$clonal + df$unique
  df$clonal.percent <- round(df$clonal/df$total*100, 2)
  df$unique.percent <- round(df$unique/df$total*100, 2)
  general.save.data(df, "RESULTS/MutationVennDiagrams/MutationNumbersClonalUnique")
  
  df.plot <- as.data.frame(pivot_longer(df, names_to = "Category", values_to = "Percent", cols = c("clonal.percent", "unique.percent")))
  df.plot$Category <- gsub("[.]percent", "", df.plot$Category)
  
  ggplot(df.plot, aes(x = Category, y = Percent)) +
    geom_violin() + 
    stat_summary(fun = "mean", geom = "crossbar", color = "black") + 
    geom_beeswarm() +
    stat_compare_means(method = "t.test", paired = TRUE)
  ggsave("RESULTS/MutationVennDiagrams/Violin_MutationNumbersClonalUnique.png", type = "cairo", width = 3, height = 4)
  ggsave("RESULTS/MutationVennDiagrams/Violin_MutationNumbersClonalUnique.pdf", width = 3, height = 4)
}


.plotVenn <- function(mutationIDsList, out){
  png(paste0("RESULTS/MutationVennDiagrams/", patient, ".png"), width = 5*300, height = 5*300, res = 300)
  venn(mutationIDsList, snames = names(mutationIDsList), zcolor = "style")
  dev.off()
  if(!grepl("[.]pdf", out)){
    out <- paste0(out, ".pdf")
  }
  pdf(out, width = 3, height = 3)
  venn(mutationIDsList, snames = names(mutationIDsList), zcolor = "style")
  dev.off()
  png(gsub("[.]pdf$", ".png", out), width = 3*300, height = 3*300, res = 300, type = "cairo")
  venn(mutationIDsList, snames = names(mutationIDsList), zcolor = "style")
  dev.off()
}
