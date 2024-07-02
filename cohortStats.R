cat("Project ID: OAC_Heterogeneity_MR_MT\n")
## R script to extract cohort statistics ##
cat("Loading: R script to extract cohort statistics\n")

## Loading external libraries ##

## Loading external scripts ##
## Load script with general project settings and background functionality
source("SCRIPTS/setup.R")

## Environment ##
anno <- general.load.data("DATA/Annotation/Cohort.csv")
anno <- anno[which(anno$Post_Samples != 0),]
total <- nrow(anno)

## Internally called functions ##
.percent <- function(i){
  return(round(i/total*100, 1))
}

## Functions ##
stats <- function(){

  cat("---\n")
  ## Age at diagnosis
  cat(paste0("Age:\nmedian: ", median(anno$Age), "; range: ", paste(range(anno$Age), collapse = "-"), "; mean: ", mean(anno$Age), "\n"))
  cat("---\n")
  
  ## Gender distribution
  tmp <- table(anno$Sex)
  cat(paste0("Gender:\n", names(tmp)[1], ": ", tmp[1], " (", .percent(tmp[1]), "%); ", names(tmp)[2], ": ", tmp[2], " (", .percent(tmp[2]), "%)\n"))
  cat("---\n")
  
  ## Clinical stage (cTNM)
  tmp <- table(gsub("[ABC]", "", anno$cTNM.8th.AJCC))
  cat("cTNM (8th AJCC edition):\n")
  a <- sapply(names(tmp), function(n){
    cat(paste0(n, ": ", tmp[n], " (", .percent(tmp[n]), "%); "))  
  })
  cat("\n---\n")
  
  ## Treatment distribution
  tmp <- table(anno$Treatment)
  cat("Treatment:\n")
  a <- sapply(names(tmp), function(n){
    cat(paste0(n, ": ", tmp[n], " (", .percent(tmp[n]), "%); "))  
  })
  cat("\n---\n")
 
  ## Pathological stage (ypTNM)
  tmp <- table(gsub("[ABC]", "", anno$ypTNM.8th.AJCC))
  cat("ypTNM (8th AJCC edition):\n")
  a <- sapply(names(tmp), function(n){
    cat(paste0(n, ": ", tmp[n], " (", .percent(tmp[n]), "%); "))  
  })
  cat("\n---\n")
  
  ## Follow-up
  alive <- anno[anno$OS == 0,]
  dead <- anno[anno$OS == 1,]
  cat(paste0("Follow-up (alive):\nn = ", nrow(alive), " (", .percent(nrow(alive)), "%); median: ", median(alive$OS.surgery), "; range: ", paste(range(alive$OS.surgery), collapse = "-"), "; mean: ", round(mean(alive$OS.surgery), 1), "\n"))
  cat(paste0("Follow-up (dead):\nn = ", nrow(dead), " (", .percent(nrow(dead)), "%); median: ", median(dead$OS.surgery), "; range: ", paste(range(dead$OS.surgery), collapse = "-"), "; mean: ", round(mean(dead$OS.surgery), 1), "\n"))
  cat("---\n")
  
  ## Recurrence
  cat("Recurrence:\n")
  rec <- anno$DOCTOR.site.of.recurrence
  rec[is.na(rec) | rec == "Dead without disease"] <- "No recurrence"
  tmp <- table(rec)
  cat("Recurrence location:\n")
  a <- sapply(names(tmp), function(n){
    cat(paste0(n, ": ", tmp[n], " (", .percent(tmp[n]), "%); "))  
  })
  cat("\n---\n")
  
  ## PET response
  cat("Recurrence:\n")
  pet <- anno$PET.response
  pet[is.na(pet)] <- "unknown"
  tmp <- table(pet)
  cat("PET response:\n")
  a <- sapply(names(tmp), function(n){
    cat(paste0(n, ": ", tmp[n], " (", .percent(tmp[n]), "%); "))  
  })
  cat("\n---\n")
  
  ## Pathological response
  cat(paste0("Pathological response (%):\nmedian: ", median(anno$Percentage.Response, na.rm = TRUE), "; range: ", paste(range(anno$Percentage.Response, na.rm = TRUE), collapse = "-"), "; mean: ", round(mean(anno$Percentage.Response, na.rm = TRUE), 1), "\n"))
  cat("---\n")
  
  ## Mandard score
  cat("Mandard score:\n")
  pResp <- anno$Mandard.simple
  pResp[is.na(pResp)] <- "unknown"
  tmp <- table(pResp)
  a <- sapply(names(tmp), function(n){
    cat(paste0(n, ": ", tmp[n], " (", .percent(tmp[n]), "%); "))  
  })
  cat("\n---\n")
}








