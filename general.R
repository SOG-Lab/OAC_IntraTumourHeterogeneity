## R script with general methods used in more than one package ##
cat("Loading: R script with general methods and variables.\n")

## Loading external libraries ##
library(tidyr) ## Data frame manipulation (e.g. pivot_longer())
library(data.table) ## for fread()
library(ggplot2) ## for general plots
library(RColorBrewer) ## for general colour use

## Loading external scripts ##

## Environment ##
dir.projects <- "PROJECTS/"
dir.data <- "DATA/"
prefix <- ""
suffix <- ""
nRowTesting <- 50
p.sig.value <- 0.05
p.fdr.sig.value <- 0.05
tumour.cols <- c(paste0(rep(c("Tumo", "tumo", "Primary tumo", "Primary Tumo", "primary tumo"), each = 2), c("r", "ur")), "TM", "TL", "TU")
chromosomeNames <- c(paste0("chr", c(1:22, "X", "Y")))
chromosomeArmNames <- paste(rep(chromosomeNames, each = 2), c("p", "q"), sep = ".")
cna.cellularity.cutoff <- 39
signatures.seqArtefacts <- paste0("SBS", c(27, 43, 45:60, 95))

## Functions ##

## Load the data from the input
## Required input 'input': Input containing the data
## Optional input 'silent': If TRUE, no output will be printed
## Return 'data': Loaded dataframe
general.load.data <- function(input, silent = FALSE, tbl = FALSE, stringsAsFactors = FALSE, header = NA, ...){
  ## Check whether the input is a dataframe or a file and load the data into a dataframe
  if(is.data.frame(input)){
    if(!silent){
      cat("Loading data from input dataframe ... ")
    }
    data <- input
  }else if(is.matrix(input)){
    if(!silent){
      cat("Loading data from input matrix ... ")
    }
    data <- as.data.frame(input)
  }else if(endsWith(input, ".gz") & !grepl("maf", input) & !grepl("vcf", input)){
    if(!silent){
      cat("Loading data from file", input,"... ")
    }
    data <- read.table(input, header = FALSE, fill = TRUE, sep = "\t", ...)
  }else if(endsWith(input, ".txt") & file.exists(input)){
    if(!silent){
      cat("Loading data from file", input,"... ")
    }
    if(is.na(header)){
      header = TRUE
    }
    ## SNP array data
    if(grepl("[.]array[.]", input) | grepl("_FinalReport", input)){
      data <- read.table(input, header = header, fill = TRUE, sep = "\t", stringsAsFactors = stringsAsFactors, skip = 10, ...)
    }else{
      data <- read.table(input, header = header, fill = TRUE, sep = "\t", stringsAsFactors = stringsAsFactors, ...)
    }
  }else if((endsWith(input, ".csv") | endsWith(input, ".CSV")) & file.exists(input)){
    if(!silent){
      cat("Loading data from file", input,"... ")
    }
    data <- read.csv(input, na.strings = c("NA", "", "n/a", "na", "#N/A", "N/A", ".", " ", "[Not Applicable]", "[Not Available]", "[Unknown]", "#NULL!", "[.]", "'--", "#VALUE!", "---", "#DIV/0!"), stringsAsFactors = stringsAsFactors, ...)
  }else if((grepl(".RData", input, ignore.case = TRUE) | grepl(".Rdat", input, ignore.case = TRUE)) & file.exists(input)){
    if(!silent){
      cat("Loading data from file", input,"... ")
    }
    data <- get(load(input))
  }else if((endsWith(input, ".maf") | endsWith(input, ".maf.gz")) & file.exists(input)){
    data <- as.data.frame(fread(input, skip = "Hugo_Symbol", stringsAsFactors = FALSE, showProgress = FALSE))
  }else if((endsWith(input, ".vcf") | endsWith(input, ".vcf.gz")) & file.exists(input)){
    if(!silent){
      cat("Loading data from file", input,"... ")
    }
    data <- as.data.frame(fread(input, skip = "#CHROM", fill = TRUE, header = TRUE, stringsAsFactors = FALSE))
  }else if(endsWith(input, ".tsv") & file.exists(input)){
    if(!silent){
      cat("Loading data from file", input,"... ")
    }
    data <- read.table(input, sep = "\t", header = TRUE)
  }else if((endsWith(input, ".dcc") | endsWith(input, ".dcc.gz")) & file.exists(input)){
    data <- as.data.frame(fread(input, skip = "analysis_id", stringsAsFactors = FALSE, showProgress = FALSE))
    if(!silent){
      cat("Loading data from file", input,"... ")
    }
  }else{
    stop(paste0("The input ", input, " does not have a supported format ('.csv', '.tsv', .txt', '.RData', '.dcc', '.maf', '.vcf', dataframe') or it doesn't exist.\n\n"))
  }
  if(!silent){
    if(is.data.frame(data) | is.matrix(data)){
      cat("Done:", nrow(data), "rows and", ncol(data), "columns\n")
    }else if(is.list(data) | is.vector(data)){
      cat("Done:", length(data), "entries\n")
    }else{
      cat("Done\n")
    }
  }
  if(tbl & is.data.frame(data)){
    data <- as_tibble(data)
  }
  invisible(data)
}

## Save the given data in a given file 
## Required input 'data': Data that schould be saved
## Required input 'out': Path to output file
## Optional input 'row.names': if FALSE, no row names will be exported in .csv
## Optional input 'csv': if FALSE, no additional .csv output will be generated
general.save.data <- function(data, out, row.names = FALSE, csv = TRUE, RData = TRUE, silent = FALSE){
  saved <- FALSE
  if(endsWith(out, ".csv")){
    out <- gsub("csv", "RData", out)
  }
  if(!endsWith(out, ".RData")){
    out <- paste0(out, ".RData")
  }
  if(RData){
    if(!silent){
      cat("Writing output to", out, "... ")
    }
    save(data, file = out)
    saved <- TRUE
    if(!silent){
      cat("Done.\n")
    }
  }
  if(csv && (is.data.frame(data) | is.matrix(data) | is.vector(data))){
    out <- gsub("RData", "csv", out)
    if(!silent){
      cat("Writing output to", out, "... ")
    }
    write.csv(data, file = out, row.names = row.names)
    saved <- TRUE
    if(!silent){
      cat("Done.\n")
    }
  }  
  if(!saved && !silent){
    warning("Sorry, data could not be saved.\n")
  }
}

## Internally called methods ##

.getLastNCharacters <- function(x, n){
  x <- as.character(x)
  return(substr(x = x, start = nchar(x)-n+1, stop = nchar(x)))
}

.setDateInsert <- function(date){
  if(date){
    insert <- gsub("-", ".", paste0("_", Sys.Date()))
  }else{
    insert <- ""
  }
  return(insert)
}

.checkAndCreateDir <- function(directory){
  if(!file.exists(directory)){
    dir.create(directory, recursive = TRUE)
  }
}

.extractPvalue <- function(s, mode){
  test <- 1
  if(mode == "survdiff"){
    ## from survival:::print.survdiff
    etmp <- s$exp
    tmp <- (sum(1 * (etmp > 0))) - 1
    p <- signif(1 - pchisq(s$chisq, tmp), digits = 2)
  }else if(mode == "survreg"){
    ## from survival:::print.survreg
    tmp <- sum(s$df) - s$idf
    chi <- 2 * diff(s$loglik)
    p <- signif(1 - pchisq(chi, tmp), digits = 2)
  }else if(mode == "coxph"){
    ## Wald test
    p <- signif(unlist(summary(s))$waldtest.pvalue, digits = 2)    
    ## Likelihood ratio test
    ## p <- signif(unlist(summary(s))$logtest.pvalue, digits = 2)
    ## Score (logrank) test
    ## p <- signif(unlist(summary(s))$sctest.pvalue, digits = 2)
    ## alternative from survival:::print.coxph
    ## se <- sqrt(diag(s$var))
    ## coef <- s$coefficients
    ## p <- signif(1 - pchisq((coef/se)^2, 1), digits = 2)
  }else{
    stop("Sorry, invalid mode.\n")
  }
  return(p)
}

.addString <- function(text, add){
  if(grepl(add, text)){
    return(text)
  }else{
    return(paste0(text, add))
  }
}
