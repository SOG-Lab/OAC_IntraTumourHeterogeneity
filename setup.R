## R script to set up the project environment for the multi-region, multi-timepoint heterogeneity and tumour evolution study ##
cat("Loading: R script to set up the project environment for the multi-region, multi-timepoint heterogeneity and tumour evolution study\n")

## Loading external libraries ##
library(ggplot2) ## Plots
library(ggpubr) ## Advanced plot settings e.g. ggarrange()
library(ggbeeswarm) ## replacement for geom_jitter()

## Loading external scripts ##
## Load script with general functions that are used across multiple projects
source("/lustre/tsbrosda/SCRIPTS/general.R")

## Environment ##
setwd("/lustre/tsbrosda/PROJECTS/OAC_Heterogeneity_MR_MT/")
dir.data <- "/lustre/tsbrosda/PROJECTS/OAC_Heterogeneity_MR_MT/DATA/"
## Load dbSnp mutations for filtering
dbSNP <- general.load.data(paste0(dir.data, "knownSNPs.RData"))
## Additional pre treatment samples to be excluded for more homogenous sample cohort
exclude.pre <- c("PAH001_TU", "PAH117_TL") ## Remove lowest cellularity sample
exclude.pre <- c("PAH001_TM", "PAH117_TM") ## Remove sample without matching RNAseq
prepost <- c("PAH001", "PAH003", "PAH008", "PAH009", "PAH025", "PAH040", "PAH053", "PAH096", "PAH113", "PAH125")
options(width = 190)

.prepareCohortFiles <- function(file){
  cohort.files <- general.load.data(file) ## Contains sample information and file storage locations
  cohort.files <- cohort.files[!is.na(cohort.files$Patient),] ## Remove empty rows that sometimes get attached at the bottom when saving from excel
  cohort.files$Sample.type[is.na(cohort.files$Sample.type)] <- "NA" ## Sample type was set to <NA> during import
  cohort.files <- cohort.files[!is.na(cohort.files$MAF.folder),]
  cohort.files <- cohort.files[order(cohort.files$Number), ]
  cohort.files$Patient <- factor(cohort.files$Patient, levels = unique(cohort.files$Patient))
  
  return(cohort.files)
}


## Filter genes by removing pseudo-genes, DBSNPs and unknown gene regions
.filterVariants <- function(data, dbSNP = NULL, geneFiltering = FALSE, snpOnly = FALSE){
  ## Remove known dbSNPs
  if(!is.null(dbSNP)){
    data <- data[!paste(data$Chromosome, data$Start_Position, sep = "_") %in% dbSNP$ID | data$Hugo_Symbol == "TP53",]
    if("DbSNP_RS" %in% names(data)){
      data <- data[data$DbSNP_RS == "novel" | data$Hugo_Symbol == "TP53",]
    }
    cat("removed known dbSNPs:", dim(data), "\n")
  }
  
  ## Keep only SNPs
  if(snpOnly){
    data <- data[data$Variant_Type == "SNP",]
  }
  
  ## Remove variants that are not on the autosomal chromosomes (1-22)
  data$Chromosome[!grepl("chr", data$Chromosome)] <- paste0("chr", data$Chromosome[!grepl("chr", data$Chromosome)])
  data <- data[data$Chromosome %in% paste0("chr", c(1:22)),]
  #cat("removed other chromosomes:", dim(data), "\n")
  
  ## Filtering for unknown and pseudo genes if requested
  ## This makes sense for RNAseq experiments but not necessarily if investigating mutational signatures etc.
  if(geneFiltering){
    ## Remove unknown genes
    data <- data[data$Hugo_Symbol != "unknown",]
    #cat("remove 'unknown':", dim(data), "\n")
    ## Remove further pseudo- and uncharacterised genes
    data <- data[!grepl("[.]", data$Hugo_Symbol),]
    exclude <- grep("[-]", data$Hugo_Symbol, value = TRUE)
    exclude <- exclude[!grepl("HLA-", exclude) & !grepl("BCL2L2", exclude) & !grepl("TRIM", exclude)]
    data <- data[!data$Hugo_Symbol %in% exclude,]
    #cat("remove pseudo-genes:", dim(data), "\n")
    ## Remove further genes
    data <- data[!grepl("mir", data$Hugo_Symbol),] 
    data <- data[!grepl("^LINC", data$Hugo_Symbol),]
  }
  return(data)
}

blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

