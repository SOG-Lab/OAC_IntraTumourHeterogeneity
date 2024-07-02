cat("Project ID: OAC_Heterogeneity_MR_MT\n")
## R script to prepare PyClone VI files ##
cat("Loading: R script to prepare PyClone VI files\n")

## Loading external librarys ##
library(GenomicRanges) ## Mapping of genomic features
library(tidyr) ## Data frame manipulation (e.g. pivot_longer())

## Loading external scripts ##
## Load script with general project settings and background functionality
source("SCRIPTS/setup.R")

## Environment ##
pileupFiles <- list.files(paste0(dir.data, "Pileups/"), full.names = TRUE, pattern = "[.][tc]") ## only grab .txt or .csv files
bases <- c("A", "C", "G", "T")
dir.out <- "DATA/ForPyClone/"
.checkAndCreateDir(dir.out)

## Internally called functions
.preparePileupFile <- function(file, insert = ""){
  ## Load data
  df <- general.load.data(file, silent = TRUE)
  donor <- sapply(strsplit(basename(file), "[.]"), "[", 1)

  if(any(grepl("POST_", df$SnpId, ignore.case = TRUE))){
    cat(" | pre only")
    df_pre <- df[!grepl("POST_", df$SnpId, ignore.case = TRUE) & !grepl("POST", df$Donor, ignore.case = TRUE), ]
    out <- paste0(dirname(file), "/", donor, ".snp.pileup.preOnly")
    general.save.data(df_pre, out, silent = TRUE)
    .preparePileupFile(paste0(out, ".csv"), insert = "_preOnly")
    cat(" | all")
  }
  ## Create SNP ID without sample and gene ID
  df$AltBase <- .getLastNCharacters(df$SnpId, 1)
  df$Id <- paste(df$Chromosome, df$Start, df$End, df$RefBase, df$AltBase, sep = "_")
  tmp <- strsplit(df$SnpId, "_")
  a <- sapply(tmp, function(el){
    len <- length(el)
    alt <- el[[len]]
    hugo <- el[[len-5]]
    return(list(alt = alt, hugo = hugo))
  })
  a <- as.data.frame(t(a))
  df$AltBase <- as.character(a$alt)
  df$Hugo_Symbol <- as.character(a$hugo)
  
  df$SnpId <- paste(df$Hugo_Symbol, df$Chromosome, df$Start, df$End, df$RefBase, df$AltBase, sep = "_")
  samples <- unique(df$Donor)
  
  ## Remove duplicated entries
  df <- unique(df)
  
  ## Add variant information that will be helpful later
  ## Get count for exact mutation base
  a <- sapply(bases, function(base){
    df$alt.count[df$AltBase == base] <<- df[df$AltBase == base, paste0(base, "plus")]
  })
  ## Recalculate total non ref counts
  df$total.count <- df$alt.count + df$TotalRef
  ## Calculate VAF
  df$vaf <- df$alt.count/df$total.count
  
  ## Remove unnecessary columns
  df <- df[, c("Hugo_Symbol", "Chromosome", "Start", "End", "RefBase", "AltBase", "SnpId", "Donor", "TotalRef", "alt.count", "total.count", "vaf")]
  names(df) <- c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Reference", "AltBase", "id", "Donor", "ref.count", "alt.count", "total.count", "vaf")
  
  ## Change from long to wide format
  names(df) <- gsub("TotalRef", "ref.count", names(df))
  df <- as.data.frame(pivot_wider(df, values_from = c("ref.count", "alt.count", "total.count", "vaf"), names_from = "Donor", names_glue = "{Donor}.{.value}"))
  
  ## Change column order
  df <- df[, c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Reference", "AltBase", "id", paste(rep(samples, each = 4), c("ref.count", "alt.count", "total.count", "vaf"), sep = "."))]
  general.save.data(df, paste0(dirname(file), "/", donor, ".pileup.overlap", insert), silent = TRUE)
}

.formatPileup <- function(data, anno){
  ## 1 row per sample
  patient <- anno$Patient[1]
  names(data) <- gsub("PAH_", "PAH", names(data))
  dataLongRef <- as.data.frame(pivot_longer(data, names_to = "sample_id", values_to = "ref_counts", cols = grep("ref.count", grep(patient, names(data), value = TRUE), value = TRUE)))
  dataLongAlt <- as.data.frame(pivot_longer(data, names_to = "sample_id", values_to = "alt_counts", cols = grep("alt.count", grep(patient, names(data), value = TRUE), value = TRUE)))
  dataLongTotal <- as.data.frame(pivot_longer(data, names_to = "sample_id", values_to = "total_counts", cols = grep("total.count", grep(patient, names(data), value = TRUE), value = TRUE)))
  dataLongRef$sample_id <- gsub("[.]ref[.]count", "", dataLongRef$sample_id)
  dataLongAlt$sample_id <- gsub("[.]alt[.]count", "", dataLongAlt$sample_id)
  dataLongTotal$sample_id <- gsub("[.]total[.]count", "", dataLongTotal$sample_id)
  
  cols <- names(dataLongRef)[!grepl(patient, names(dataLongRef))]
  dataLong <- merge(dataLongRef[, cols], dataLongAlt[, c("id", "sample_id", "alt_counts")])
  dataLong <- merge(dataLong, dataLongTotal[, c("id", "sample_id", "total_counts")])

  a <- apply(anno, 1, function(row){
    dataLong$sample_id[dataLong$sample_id == row[["pileup_SampleID"]]] <<- row[["Sample_ID"]]    
  })
  
  return(dataLong)
}

## Add copy number data from ASCAT outputs
.addCNdata <- function(data, cnFile, cn_tool){
  cnData <- general.load.data(cnFile, silent = TRUE)

  if(nrow(cnData) == 0){
    return(data)
  }
  if(cn_tool %in% c("ascatNGS", "ASCATNGS", "ascatngs")){
    names(cnData) <- gsub("chromosome", "chr", names(cnData))
    names(cnData) <- gsub("start", "startpos", names(cnData))
    names(cnData) <- gsub("end", "endpos", names(cnData))
    names(cnData) <- gsub("tumour.minor.cn", "nMinor", names(cnData))
    names(cnData) <- gsub("tumour.major.cn", "nMajor", names(cnData))
    names(cnData) <- gsub("tumour.total.cn", "Copy", names(cnData))
  }
  if(!"Copy" %in% names(cnData)){
    cnData$Copy <- cnData$nMinor + cnData$nMajor
  }
  
  cnData <- cnData[, c("startpos", "endpos", "chr", "nMinor", "nMajor", "Copy")]

  cnData$chr[!grepl("chr", cnData$chr)] <- paste0("chr", cnData$chr[!grepl("chr", cnData$chr)])
  ## Keep autosomes only (chr 1-22)
  cnData <- cnData[cnData$chr %in% paste0("chr", 1:22),]
  ## Use GRanges to map nMinor and nMajor to mutations
  cnv.gr <- makeGRangesFromDataFrame(df = cnData, start.field = "startpos", end.field = "endpos", seqnames.field = "chr", keep.extra.columns = TRUE, ignore.strand = TRUE)
  data.gr <- makeGRangesFromDataFrame(df = data, start.field = "Start_Position", end.field = "End_Position", seqnames.field = "Chromosome", keep.extra.columns = TRUE)
  ## Find overlapping segments
  ov <- findOverlaps(query = cnv.gr, subject = data.gr)
  ## Add nMinor, nMajor and Copy
  data$minor_cn[subjectHits(ov)] <- cnv.gr[queryHits(ov)]$nMinor
  data$major_cn[subjectHits(ov)] <- cnv.gr[queryHits(ov)]$nMajor
  data$Copy[subjectHits(ov)] <- cnv.gr[queryHits(ov)]$Copy
  
  return(data)  
}
  
.addCN <- function(data, anno, cn_tool = "ascatNGS"){
  ## Add default values
  data$Copy <- data$normal_cn <- 2
  data$minor_cn <- data$major_cn <- 1

  ## Add actual CN for each sample
  dataCN <- data.frame()
  i <- 1
  total <- nrow(anno)
  a <- apply(anno, 1, function(row){
    sample_ID <- row[["Sample_ID"]]
    cat("\r", i, "/", total, "|", sample_ID)
    i <<- i + 1
    if(cn_tool == "GAP"){
      file <- paste(row[["CNA.GAP.folder"]], row[["CNA.GAP.file"]], sep = "/")
    }else if(cn_tool == "ASCAT"){
      file <- paste(row[["CNA.ASCAT.folder"]], row[["CNA.ASCAT.file"]], sep = "/")   
    }else if(cn_tool %in% c("ascatNGS", "ASCATNGS", "ascatngs")){
      file <- paste(row[["CNA.ascatNGS.folder"]], row[["CNA.ascatNGS.file"]], sep = "/")   
    }

    df <- .addCNdata(data[data$sample_id == sample_ID,], file, cn_tool)
    if(nrow(dataCN) == 0){
      dataCN <<- df
    }else{
      dataCN <<- rbind(dataCN, df)
    }
  })
#  cat("\n")
  return(dataCN)
}

.addTumourContent <- function(data, anno){
  ## Default tumour content
  data$tumour_content <- 1
  a <- sapply(anno$Sample_ID, function(sample){
    data$tumour_content[data$sample_id == sample] <<- anno$Cellularity[anno$Sample_ID == sample]/100
  })
  return(data)
}


correctPileupFormat <- function(){
  toFix <- pileupFiles[!grepl("overlap", pileupFiles) & !grepl("preOnly", pileupFiles)]#[c(1,8)]
  i <- 1
  total <- length(toFix)
  
  a <- sapply(toFix, function(file){
    donor <- sapply(strsplit(basename(file), "[.]"), "[", 1)
    if(grepl("preOnly", file)){
      donor <- paste0(donor, " | preOnly")
    }
    cat("\n", i, "/", total, "|", donor)
    i <<- i + 1
    .preparePileupFile(file)
  }) 
  cat("\n")
}

addPreOnly <- function(files){
  patient <- patients[1]
  annoSamples <- general.load.data("DATA/Annotation/Cohort_sampleDetails.csv")
  annoSamples <- annoSamples[annoSamples$Patient %in% patients & !is.na(annoSamples$Sample.type) & annoSamples$Sample.type != "BC" & annoSamples$Timepoint %in% c("Pre", "pre"),]
  sapply(patients, function(patient){
    anno <- annoSamples[annoSamples$Patient == patient, ]
    mutationIDs <- c()
    a <- sapply(paste(anno$MAF.folder, anno$MAF.file, sep = "/"), function(file){
      m <- general.load.data(file)
      mutationIDs <<- unique(c(mutationIDs, paste(m$Hugo_Symbol, paste0("chr", m$Chromosome), m$Start_Position, m$Start_Position, m$Reference_Allele, m$Tumor_Seq_Allele2, sep = "_")))
    })
    pile <- pile[pile$id %in% mutationIDs,]
    general.save.data(pile, paste0(dir.data, "Pileups/", patient, ".pileup.overlap_preOnly"))
  })
}


preparePileupFiles <- function(cn_tool = "ascatNGS"){
  pileupFiles <- pileupFiles[grepl("overlap", pileupFiles)]
 
  anno <- general.load.data("DATA/Annotation/Cohort.csv")
  annoSamples <- general.load.data("DATA/Annotation/Cohort_sampleDetails.csv")
  sampleMapped <- general.load.data("DATA/Annotation/SampleID_mapping.csv")
  
  patients <- unique(anno$Patient)
   i <- 1
  total <- length(patients)
  
  ## For each patient go through pileup file, filter, add copy number data and save
  a <- sapply(patients[i:total], function(patient){
    cat(i, "/", total, "|", patient, "\n")
    i <<- i + 1
    ## Get pileup file location
    tmp <- annoSamples[annoSamples$Patient == patient & !annoSamples$Sample.type %in% c("BC", "NA") & !is.na(annoSamples$Sample.type),]
    tmp <- merge(tmp, sampleMapped, by.x = c("Number", "Patient", "Timepoint", "Sample_ID"), by.y = c("Number", "Patient", "Timepoint", "clean_Sample_ID"))
    ## Load pileup data
    cat(" loading ...")
    pFile <- paste0(tmp$Pileup.folder[1], tmp$Pileup.file[1])
    pileup <- general.load.data(pFile, silent = TRUE)
    ## Filter
    cat(" filtering ...")
    pileup <- .filterVariants(pileup) ## from setup.R
    ## Format and amend sample IDs to clean sample IDs
    cat(" formatting ...")
    pileupLong <- .formatPileup(pileup, tmp)
    ## Add copy number data
    cat(" adding CN data ...")
    pileupLong <- .addCN(pileupLong, tmp, cn_tool)
    ## Add tumour content (0-1)
    pileupLong <- .addTumourContent(pileupLong, tmp)
    ## Save data
    cat(" saving ...")
    names(pileupLong)[names(pileupLong) == "id"] <- "mutation_id"
    cols <- c("mutation_id", "sample_id", "ref_counts", "alt_counts", "major_cn", "minor_cn", "normal_cn", "tumour_content")
    general.save.data(pileupLong, paste0(dir.out, patient, "_pileup_forPyClone"), silent = TRUE)
    write.table(pileupLong[, cols], paste0(dir.out, patient, "_pileup_forPyClone_", cn_tool, ".tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
    ## If pre-only file exists do the same with that file
    if((endsWith(pFile, "txt") & file.exists(gsub("[.]txt", "_preOnly.csv", pFile))) | (endsWith(pFile, "csv") & file.exists(gsub("[.]csv", "_preOnly.csv", pFile)))){
      cat(" PRE only ...")
      if(endsWith(pFile, "txt")){
        pFile <- gsub("[.]txt", "_preOnly.csv", pFile)
      }else{
        pFile <- gsub("[.]csv", "_preOnly.csv", pFile)
      }
      pOnly <- general.load.data(pFile, silent = TRUE)
      pileupLong <- pileupLong[pileupLong$mutation_id %in% pOnly$id & !grepl("POST", pileupLong$sample_id, ignore.case = TRUE), ]
      general.save.data(pileupLong, paste0(dir.out, patient, "_pileup_forPyClone_preOnly_", cn_tool), silent = TRUE)
      write.table(pileupLong[, cols], paste0(dir.out, patient, "_pileup_forPyClone_preOnly_", cn_tool, ".tsv"), quote = FALSE, sep = "\t", row.names = FALSE)
      pFile <- gsub("_preOnly.csv", ".csv", pFile)
    }
    cat(" done\n")
  })
}  
  
  
  
  
  
  
