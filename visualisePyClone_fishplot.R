cat("Project ID: OAC_Heterogeneity_MR_MT\n")
## R script to visualise the PyClone results as Fishplots##
cat("Loading: R script to visualise the PyClone results as Fishplots\n")

## Loading external libraries ##
library(ggplot2)
library(fishplot) ## For fishplots (https://github.com/chrisamiller/fishplot)
library(RColorBrewer) ## For additional (larger) colour palettes

## Loading external scripts ##
source("SCRIPTS/setup.R")

## Environment ##
## Directory for all fishplots
dir.plots <- "RESULTS/PyCloneVI/Fishplots/"
.checkAndCreateDir(dir.plots)

## Internally called functions ##
.plotFish <- function(frac.table, parents, donor, mode, insert = "", swap = FALSE){
  if(mode == "5samples"){
    timepoints <- c(0, 10, 50, 60, 100, 110, 150, 160, 200, 210)
    lines <- c(5, 55, 105, 155, 205)
    frac.table <- frac.table[, c(1,1,2,2,3,3,4,4,5,5)]
    if(swap){
      frac.table <- frac.table[, 10:1] 
    }
  }else if(mode == "4samples"){
    timepoints <- c(0, 10, 50, 60, 100, 110, 150, 160)
    lines <- c(5, 55, 105, 155)
    frac.table <- frac.table[, c(1,1,2,2,3,3,4,4)]
    if(swap){
      frac.table <- frac.table[, 8:1] 
    }
  }else if(mode == "3samples"){
    timepoints <- c(0, 10, 50, 60, 100, 110)
    lines <- c(5, 55, 105)
    frac.table <- frac.table[, c(1,1,2,2,3,3)]
    if(swap){
      frac.table <- frac.table[, c(5, 5, 3, 3, 1, 1)] 
    }
  }else if(mode == "2samples"){
    timepoints <- c(0, 10, 50, 60)
    lines <- c(5, 55)
    frac.table <- frac.table[, c(1,1,2,2)]
    if(swap){
      frac.table <- frac.table[, c(3, 3, 1, 1)] 
    }
  }else if(mode == "2samples_extend"){
    timepoints <- c(0, 10, 50, 60, 100, 110, 150, 160)
    lines <- c(30, 130)
    frac.table <- frac.table[, c(1, 1, 1, 1, 2, 2, 2, 2)]
    if(swap){
      frac.table <- frac.table[, 8:1] 
    }
  }else if(mode == "1sample"){
    timepoints <- c(0, 10)
    lines <- 5
  }
  ## Create fish
  fish <- createFishObject(frac.table, parents, timepoints = timepoints, fix.missing.clones = TRUE)
  fish <- layoutClones(fish)
  if(nrow(frac.table) == 15){
    fish <- setCol(fish, c("darkgrey", brewer.pal(nrow(frac.table), "Set3"), "darkorange1", "yellowgreen"))
  }else if(nrow(frac.table) == 13){
    fish <- setCol(fish, c("darkgrey", brewer.pal(nrow(frac.table), "Set3")))
  }else if(nrow(frac.table) > 10){
    fish <- setCol(fish, brewer.pal(nrow(frac.table), "Set3"))
  }
  samples <- colnames(frac.table)
  insert <- paste(unique(colnames(frac.table)), collapse = ".")

  ## Plotting
  png(paste0(dir.plots, donor, "_fishplot_", insert, ".png"), width = 5*300, height = 3*300, res = 300, type = "cairo")
  fishPlot(fish,shape = "spline", title.btm = donor,
           cex.title = 0.5, vlines = lines, 
           vlab = unique(samples))
  dev.off()
  pdf(paste0(dir.plots, donor, "_fishplot_", insert, ".pdf"), width = 5, height = 3)
  fishPlot(fish,shape = "spline", title.btm = donor,
           cex.title = 0.5, vlines = lines, 
           vlab = unique(samples))
  dev.off()
}

## Same names for all pre post tumours
samples.prePost <- c("Pre", "Post")


## Plotted samples
cat("\n-------------------------------\n")
patient <- "OESO_0001"
cat(patient)
samples <- c("TL", "POST_TL", "TM", "POST_TM", "TU")
frac.table <- matrix(c(100, 100, 100, 100, 100,# C1
                       20, 25, 20, 20, 20, # C7
                       70, 65, 70, 0, 70, # C10
                       2, 55, 0, 0, 25, # C8
                       0, 30, 0, 0, 5, # C6
                       50, 6, 65, 0, 15, # C9
                       5, 0, 30, 0, 10, # C2
                       
                       40, 4, 25, 0, 0, # C4
                       
                       25, 2, 5, 0, 0, # C3
                       0, 0, 0, 70, 0 # C5
), ncol = length(samples), byrow = TRUE)
colnames(frac.table) <- samples
frac.table <- frac.table[, c(5,3,1,2,4)]
## clusterID 1  7 10  8  6  9  2  4  3  5
parents <- c(0, 1, 1, 3, 4, 3, 6, 6, 8, 1)
.plotFish(frac.table, parents, patient, mode = "5samples")
.plotFish(frac.table, parents, patient, mode = "5samples", swap = TRUE)
## Combine pres and posts
pre <- (frac.table[, "TU"] + frac.table[, "TM"] + frac.table[, "TL"]) / 3
post <- (frac.table[, "POST_TL"] + frac.table[, "POST_TM"]) / 2
frac.table <- matrix(c(pre, pre, post, post), ncol = 4, byrow = FALSE)
colnames(frac.table) <- rep(samples.prePost, each = 2)
.plotFish(frac.table, parents, patient, mode = "4samples")

cat("\n-------------------------------\n")
patient <- "OESO_0003"
cat(patient)
samples <- c("POST_TL", "TM", "POST_TM", "TU")
frac.table <- matrix(c(100, 100, 100, 100, # C1
                       30, 95, 15, 95, # C4
                       0, 50, 0, 35, # C2
                       65, 0, 80, 0, # C6
                       0, 0, 55, 0, # C5
                       50, 0, 0, 0 # C3
), ncol = length(samples), byrow = TRUE)
colnames(frac.table) <- samples
frac.table <- frac.table[, c(4,2,3,1)]
## clusterID 1  4  2  6  3  5 
parents <- c(0, 1, 2, 1, 4, 4)
.plotFish(frac.table, parents, patient, mode = "4samples")
.plotFish(frac.table, parents, patient, mode = "4samples", swap = TRUE)
## Combine pres and posts
pre <- (frac.table[, "TU"] + frac.table[, "TM"]) / 2
post <- (frac.table[, "POST_TL"] + frac.table[, "POST_TM"]) / 2
frac.table <- matrix(c(pre, pre, post, post), ncol = 4, byrow = FALSE)
colnames(frac.table) <- rep(samples.prePost, each = 2)
.plotFish(frac.table, parents, patient, mode = "4samples")

cat("\n-------------------------------\n")
patient <- "OESO_0008"
cat(patient)
samples <- c("TL", "POST_TM", "TU", "POST_TU")
frac.table <- matrix(c(100, 100, 100, 100, # C1
                       45, 50, 90, 90, # C9
                       0, 40, 0, 0, # C5
                       50, 2, 0, 0, # C8
                       18, 10, 20, 30, # C6    
                       5, 5, 5, 5, # C4
                       22, 20, 25, 2, # C7     
                       0, 0, 20, 0, # C10
                       2, 15, 40, 55, # C2    
                       0, 0, 35, 50, # C3
                       0, 0, 2, 45 # C11
), ncol = length(samples), byrow = TRUE)
colnames(frac.table) <- samples
frac.table <- frac.table[, c(1,3,4,2)]
## clusterID 1  9  5  8  6  4  7 10  2  3  11
parents <- c(0, 1, 1, 1, 2, 5, 2, 7, 2, 9, 10)
.plotFish(frac.table, parents, patient, mode = "4samples")
.plotFish(frac.table, parents, patient, mode = "4samples", swap = TRUE)
## Combine pres and posts
pre <- (frac.table[, "TU"] + frac.table[, "TL"]) / 2
post <- (frac.table[, "POST_TU"] + frac.table[, "POST_TM"]) / 2
frac.table <- matrix(c(pre, pre, post, post), ncol = 4, byrow = FALSE)
colnames(frac.table) <- rep(samples.prePost, each = 2)
.plotFish(frac.table, parents, patient, mode = "4samples")

cat("\n-------------------------------\n")
patient <- "OESO_0009"
cat(patient)
samples <- c("TL", "POST_TL", "TU", "POST_TU")
frac.table <- matrix(c(100, 100, 100, 100, # C1
                       25, 80, 0, 90, # C7
                       15, 25, 0, 10, # C3
                       10, 15, 0, 0, # C8
                       0, 0, 0, 65, # C6
                       70, 15, 95, 5, # C4
                       0, 0, 85, 0, # C9
                       0, 0, 50, 0, # C2
                       50, 0, 0, 2 # C5
), ncol = length(samples), byrow = TRUE)
colnames(frac.table) <- samples
frac.table <- frac.table[, c(3,1,2,4)]
## clusterID 1  7  3  8  6  4  9  2  5
parents <- c(0, 1, 2, 3, 2, 1, 6, 7, 6)
.plotFish(frac.table, parents, patient, mode = "4samples")
.plotFish(frac.table, parents, patient, mode = "4samples", swap = TRUE)
## Combine pres and posts
pre <- (frac.table[, "TU"] + frac.table[, "TL"]) / 2
post <- (frac.table[, "POST_TL"] + frac.table[, "POST_TU"]) / 2
frac.table <- matrix(c(pre, pre, post, post), ncol = 4, byrow = FALSE)
colnames(frac.table) <- rep(samples.prePost, each = 2)
.plotFish(frac.table, parents, patient, mode = "4samples")

cat("\n-------------------------------\n")
patient <- "OESO_0025"
cat(patient)
samples <- c("TL", "TM", "POST_TU")
frac.table <- matrix(c(100, 100, 100, # C1
                       85, 95, 95, # C6
                       75, 90, 90, # C2
                       55, 65, 85, # C8  
                       23, 30, 55, # C10 
                       5, 12, 40, # C9
                       15, 12, 10, # C4
                       27, 25, 30, # C5 
                       20, 20, 4, # C3
                       15, 2, 2 # C7
), ncol = length(samples), byrow = TRUE)
colnames(frac.table) <- samples
## clusterID 1  6  2  8 10  9  4  5  3  7 
parents <- c(0, 1, 2, 3, 4, 5, 5, 4, 8, 9)
.plotFish(frac.table, parents, patient, mode = "3samples")
.plotFish(frac.table, parents, patient, mode = "3samples", swap = TRUE)
## Combine pres and posts
pre <- (frac.table[, "TL"] + frac.table[, "TM"]) / 2
post <- frac.table[, "POST_TU"]
frac.table <- matrix(c(pre, pre, post, post), ncol = 4, byrow = FALSE)
colnames(frac.table) <- rep(samples.prePost, each = 2)
.plotFish(frac.table, parents, patient, mode = "4samples")

cat("\n-------------------------------\n")
patient <- "OESO_0040"
cat(patient)
samples <- c("TL", "POST_TM", "TU")
frac.table <- matrix(c(100, 100, 100, # C1
                       40, 35, 31, # C8
                       35, 10, 15, # C2
                       20, 0, 0,  # C10
                       10, 40, 35, # C7
                       5, 20, 30, # C11
                       0, 0, 25, # C9
                       2, 15, 0, # C6
                       45, 20, 28, # C5
                       40, 0, 0, # C3
                       20, 0, 0 # C4
                       
), ncol = length(samples), byrow = TRUE)
colnames(frac.table) <- samples
frac.table <- frac.table[, c(3,1,2)]
## clusterID 1  8  2 10  7 11  9  6  5  3  4
parents <- c(0, 1, 2, 3, 1, 5, 6, 6, 1, 9, 10)
.plotFish(frac.table, parents, patient, mode = "3samples")
.plotFish(frac.table, parents, patient, mode = "3samples", swap = TRUE)
## Combine pres and posts
pre <- (frac.table[, "TU"] + frac.table[, "TL"]) / 2
post <- frac.table[, "POST_TM"]
frac.table <- matrix(c(pre, pre, post, post), ncol = 4, byrow = FALSE)
colnames(frac.table) <- rep(samples.prePost, each = 2)
.plotFish(frac.table, parents, patient, mode = "4samples")

cat("\n-------------------------------\n")
patient <- "OESO_0053"
cat(patient)
samples <- c("TL", "POST_TL", "TU", "POST_TU") 
timepoints <- c(0, 50, 100, 150)
frac.table <- matrix(c(100, 100, 100, 100, # C1
                       25, 25, 30, 30, # C4
                       10, 10, 15, 15, # C9
                       35, 35, 55, 2, # C6
                       0, 0, 30, 0, # C3
                       35, 35, 10, 60, # C2
                       30, 30, 0, 55, # C8
                       20, 0, 0, 0, # C7
                       0, 25, 0, 45, # C10
                       0, 0, 0, 30, # C11
                       0, 10, 0, 0 # C5
), ncol = length(samples), byrow = TRUE)
colnames(frac.table) <- samples
frac.table <- frac.table[, c(1,3,4,2)]
## clusterID 1  4  9, 6  3  2  8  7 10 11  5
parents <- c(0, 1, 2, 1, 4, 1, 6, 7, 7, 9, 9)
.plotFish(frac.table, parents, patient, mode = "4samples")
.plotFish(frac.table, parents, patient, mode = "4samples", swap = TRUE)
frac.table <- frac.table[, c(2, 1, 4, 3)]
.plotFish(frac.table, parents, patient, mode = "4samples")
.plotFish(frac.table, parents, patient, mode = "4samples", swap = TRUE)
## Combine pres and posts
pre <- (frac.table[, "TU"] + frac.table[, "TL"]) / 2
post <- (frac.table[, "POST_TL"] + frac.table[, "POST_TU"]) / 2
frac.table <- matrix(c(pre, pre, post, post), ncol = 4, byrow = FALSE)
colnames(frac.table) <- rep(samples.prePost, each = 2)
.plotFish(frac.table, parents, patient, mode = "4samples")

cat("\n-------------------------------\n")
patient <- "OESO_0096"
cat(patient)
samples <- c("TM", "POST_TM", "TU", "POST_TU")
timepoints <- c(0, 50, 100, 150)
frac.table <- matrix(c(100, 100, 100, 100, # C1
                       35, 30, 35, 45, # C9 
                       25, 5, 25, 40, # C3
                       60, 30, 25, 2, # C4
                       35, 20, 15, 0, # C7
                       30, 0, 0, 0, # C5
                       2, 30, 35, 45, # C2
                       0, 20, 2, 0, # C8
                       0, 2, 15, 35 # C6
), ncol = length(samples), byrow = TRUE)
colnames(frac.table) <- samples
#frac.table <- frac.table[, c(3,1,2,4)]
frac.table <- frac.table[, c(1,3,4,2)]
## clusterID 1  9  3  4, 7  5  2  8  6 
parents <- c(0, 1, 2, 1, 4, 5, 1, 7, 7)
.plotFish(frac.table, parents, patient, mode = "4samples")
.plotFish(frac.table, parents, patient, mode = "4samples", swap = TRUE)
frac.table <- frac.table[, c(2,1,4,3)]
.plotFish(frac.table, parents, patient, mode = "4samples")
.plotFish(frac.table, parents, patient, mode = "4samples", swap = TRUE)
## Combine pres and posts
pre <- (frac.table[, "TU"] + frac.table[, "TM"]) / 2
post <- (frac.table[, "POST_TU"] + frac.table[, "POST_TM"]) / 2
frac.table <- matrix(c(pre, pre, post, post), ncol = 4, byrow = FALSE)
colnames(frac.table) <- rep(samples.prePost, each = 2)
.plotFish(frac.table, parents, patient, mode = "4samples")

cat("\n-------------------------------\n")
patient <- "OESO_0113"
cat(patient)
samples <- c("TM", "POST_TM", "TU", "POST_TU")
timepoints <- c(0, 50, 100, 150)
frac.table <- matrix(c(100, 100, 100, 100, # C1
                       40, 45, 0, 45, # C2
                       35, 40, 0, 0, # C8
                       0, 35, 0, 0, # C3
                       30, 0, 0, 0, # C9
                       0, 0, 0, 30, # C5
                       55, 50, 95, 50, # C4
                       0, 0, 55, 0, # C6
                       20, 30, 35, 25, # C10
                       10, 15, 15, 10 # C7
                       
), ncol = length(samples), byrow = TRUE)
colnames(frac.table) <- samples
frac.table <- frac.table[, c(3,1,2,4)]
## clusterID 1  2  5  8  3  9  4  6 10  7
parents <- c(0, 1, 2, 2, 4, 4, 1, 7, 7, 9)
## clusterID 1  2  8  3  9  5  4  6 10  7
parents <- c(0, 1, 2, 3, 3,2,  1, 7, 7, 9)

.plotFish(frac.table, parents, patient, mode = "4samples")
.plotFish(frac.table, parents, patient, mode = "4samples", swap = TRUE)
## Combine pres and posts
pre <- (frac.table[, "TU"] + frac.table[, "TM"]) / 2
post <- (frac.table[, "POST_TU"] + frac.table[, "POST_TM"]) / 2
frac.table <- matrix(c(pre, pre, post, post), ncol = 4, byrow = FALSE)
colnames(frac.table) <- rep(samples.prePost, each = 2)
.plotFish(frac.table, parents, patient, mode = "4samples")

cat("\n-------------------------------\n")
patient <- "OESO_0125"
cat(patient)
samples <- c("TL", "POST_TM", "TU", "POST_TU")
timepoints <- c(0, 50, 100, 150)
frac.table <- matrix(c(100, 100, 100, 100, # C1
                       95, 0, 95, 0, # C4
                       20, 0, 20, 0, # C5
                       0, 90, 0, 95, # C3
                       0, 2, 0, 50 # C2
), ncol = length(samples), byrow = TRUE)
colnames(frac.table) <- samples
frac.table <- frac.table[, c(1,3,4,2)]
## clusterID 1  4  5  3  2
parents <- c(0, 1, 2, 1, 4)
.plotFish(frac.table, parents, patient, mode = "4samples")
.plotFish(frac.table, parents, patient, mode = "4samples", swap = TRUE)
## Combine pres and posts
pre <- (frac.table[, "TU"] + frac.table[, "TL"]) / 2
post <- (frac.table[, "POST_TM"] + frac.table[, "POST_TU"]) / 2
frac.table <- matrix(c(pre, pre, post, post), ncol = 4, byrow = FALSE)
colnames(frac.table) <- rep(samples.prePost, each = 2)
.plotFish(frac.table, parents, patient, mode = "4samples")


## PRE ONLY ##
cat("\n-------------------------------\n")
patient <- "OESO_0014"
samples <- c("TL", "TM")
frac.table <- matrix(c(100, 100, # C1
                       90, 65, # C5
                       35, 25, # C3
                       0, 25, # C4
                       40, 0 # C2
), ncol = length(samples), byrow = TRUE)
colnames(frac.table) <- samples
## clusterID 1  5  3  4  2 
parents <- c(0, 1, 2, 2, 2)
.plotFish(frac.table, parents, patient, mode = "2samples_extend")
.plotFish(frac.table, parents, patient, mode = "2samples_extend", swap = TRUE)

cat("\n-------------------------------\n")
patient <- "OESO_0015"
samples <- c("TL", "TU")
frac.table <- matrix(c(100, 100, # C1
                       70, 95, # C2
                       65, 75, # C3
                       35, 30, # C5
                       0, 40, #C6
                       20, 0 # C4
                       
), ncol = length(samples), byrow = TRUE)  
colnames(frac.table) <- samples
## clusterID 1  2  3  5  6  4
parents <- c(0, 1, 2, 3, 3, 3)
.plotFish(frac.table, parents, patient, mode = "2samples_extend")
.plotFish(frac.table, parents, patient, mode = "2samples_extend", swap = TRUE)

cat("\n-------------------------------\n")
patient <- "OESO_0023"
samples <- c("TM", "TU")
frac.table <- matrix(c(100, 100, # C1
                       35, 30, # C2
                       0, 60, # C3
                       55, 5, # C4
                       45, 0 # C5
                       
), ncol = length(samples), byrow = TRUE)
colnames(frac.table) <- samples
## clusterID 1  2  3  4  5  
parents <- c(0, 1, 1, 1, 4)
.plotFish(frac.table, parents, patient, mode = "2samples_extend")
.plotFish(frac.table, parents, patient, mode = "2samples_extend", swap = TRUE)

cat("\n-------------------------------\n")
patient <- "OESO_0025_preOnly"
samples <- c("TL", "TM")
frac.table <- matrix(c(100, 100, # C1
                       85, 90, # C7
                       65, 85, # C8
                       60, 65, # C3
                       25, 15, # C4
                       30, 40, # C5
                       20, 2, # C2
                       5, 30 # C6
), ncol = length(samples), byrow = TRUE)
colnames(frac.table) <- samples
## clusterID 1  7  8  3  4  5  2  6
parents <- c(0, 1, 2, 3, 4, 4, 6, 6)
.plotFish(frac.table, parents, patient, mode = "2samples_extend")
.plotFish(frac.table, parents, patient, mode = "2samples_extend", swap = TRUE)

cat("\n-------------------------------\n")
patient <- "OESO_0001_preOnly"
samples <- c("TL", "TU")
frac.table <- matrix(c(100, 100, # C1
                       45, 15, # C5
                       20, 60, # C2
                       30, 20, # C3
                       5, 15, # C4
                       30, 0 # C6
), ncol = length(samples), byrow = TRUE)
colnames(frac.table) <- samples
## clusterID 1  5  2  3  4  6
parents <- c(0, 1, 1, 1, 3, 2)
.plotFish(frac.table, parents, patient, mode = "2samples_extend")
.plotFish(frac.table, parents, patient, mode = "2samples_extend", swap = TRUE)

cat("\n-------------------------------\n")
patient <- "OESO_0008_preOnly"
samples <- c("TL", "TU")
frac.table <- matrix(c(100, 100, # C1
                       50, 0, # C2
                       45, 90, # C5
                       0, 80, # C3
                       0, 55 # C4
), ncol = length(samples), byrow = TRUE)
colnames(frac.table) <- samples
## clusterID 1  2  5  3  4 
parents <- c(0, 1, 1, 3, 4)
.plotFish(frac.table, parents, patient, mode = "2samples_extend")
.plotFish(frac.table, parents, patient, mode = "2samples_extend", swap = TRUE)

cat("\n-------------------------------\n")
patient <- "OESO_0031"
samples <- c("TL", "TM")
frac.table <- matrix(c(100, 100, # C1
                       25, 30, # C5
                       10, 15, # C3
                       15, 0, # C2
                       2, 15 # C4
), ncol = length(samples), byrow = TRUE)
colnames(frac.table) <- samples
## clusterID 1  5  3  4  2
parents <- c(0, 1, 2, 1, 1)
.plotFish(frac.table, parents, patient, mode = "2samples_extend")
.plotFish(frac.table, parents, patient, mode = "2samples_extend", swap = TRUE)

cat("\n-------------------------------\n")
patient <- "OESO_0040_preOnly"
samples <- c("TL", "TU")
frac.table <- matrix(c(100, 100, # C1
                       40, 44, # C6
                       15, 20, # C9
                       40, 2, # C7
                       35, 0, # C3
                       20, 0, # C2
                       15, 49, # C4
                       5, 45, # C8
                       0, 30 # C5
), ncol = length(samples), byrow = TRUE)
colnames(frac.table) <- samples
## clusterID 1  6  9  7  3  2  4  8  5
parents <- c(0, 1, 2, 1, 4, 5, 1, 7, 8)
.plotFish(frac.table, parents, patient, mode = "2samples_extend")
.plotFish(frac.table, parents, patient, mode = "2samples_extend", swap = TRUE)

cat("\n-------------------------------\n")
patient <- "OESO_P0056"
samples <- c("TM", "TU")
frac.table <- matrix(c(100, 100, # C1
                       95, 75, # C3
                       70, 60, # C2
                       25, 25, # C4
                       5, 30, # C5
                       0, 25, # C6
                       35, 0 # C7
), ncol = length(samples), byrow = TRUE)
colnames(frac.table) <- samples
## clusterID 1  3  2  4  5  6  7
parents <- c(0, 1, 2, 3, 3, 5, 3)
.plotFish(frac.table, parents, patient, mode = "2samples_extend")
.plotFish(frac.table, parents, patient, mode = "2samples_extend", swap = TRUE)

cat("\n-------------------------------\n")
patient <- "OESO_0064"
samples <- c("TL", "TU")
frac.table <- matrix(c(100, 100, # C1
                       50, 7, # C2
                       30, 2, # C3
                       45, 70, # C6
                       5, 40, # C5
                       0, 25, # C7
                       30, 25 # C4
), ncol = length(samples), byrow = TRUE)
colnames(frac.table) <- samples
## clusterID 1  2  3  6  5  7  4
parents <- c(0, 1, 2, 1, 4, 5, 4)
.plotFish(frac.table, parents, patient, mode = "2samples_extend")
.plotFish(frac.table, parents, patient, mode = "2samples_extend", swap = TRUE)

cat("\n-------------------------------\n")
patient <- "OESO_0070"
samples <- c("TL", "TM")
frac.table <- matrix(c(100, 100, # C1
                       20, 15, # C4
                       0, 80, # C2
                       75, 0, # C3
                       40, 0 # C5
), ncol = length(samples), byrow = TRUE)
colnames(frac.table) <- samples
## clusterID 1  2  4  3  5  
parents <- c(0, 1, 1, 1, 4)
.plotFish(frac.table, parents, patient, mode = "2samples_extend")
.plotFish(frac.table, parents, patient, mode = "2samples_extend", swap = TRUE)

cat("\n-------------------------------\n")
patient <- "OESO_0073"
cat(patient)
samples <- c("TM", "TU")
frac.table <- matrix(c(100, 100, # C1
                       60, 95, # C3
                       50, 0, # C4
                       2, 80 # C2
), ncol = length(samples), byrow = TRUE)
colnames(frac.table) <- samples
## clusterID 1  3  4  2  
parents <- c(0, 1, 2, 2)
.plotFish(frac.table, parents, patient, mode = "2samples_extend")
.plotFish(frac.table, parents, patient, mode = "2samples_extend", swap = TRUE)

patient <- "OESO_0051"
cat(patient)
samples <- c("TL", "TU")
frac.table <- matrix(c(100, 100, # C1
                       95, 75, # C2
                       60, 30, # C4
                       50, 0, # C5
                       30, 40, # C6
                       0, 30 # C3
), ncol = length(samples), byrow = TRUE)
colnames(frac.table) <- samples
## clusterID 1  2  4  5  6  3  
parents <- c(0, 1, 2, 3, 2, 5)
.plotFish(frac.table, parents, patient, mode = "2samples_extend")
.plotFish(frac.table, parents, patient, mode = "2samples_extend", swap = TRUE)

cat("\n-------------------------------\n")
patient <- "OESO_0053_preOnly"
samples <- c("TL", "TU")
frac.table <- matrix(c(100, 100, # C1
                       45, 15, # C2
                       0, 60, # C4
                       50, 0 # C3
                       
), ncol = length(samples), byrow = TRUE)
colnames(frac.table) <- samples
## clusterID 1  2  4  3
parents <- c(0, 1, 1, 1)
.plotFish(frac.table, parents, patient, mode = "2samples_extend")
.plotFish(frac.table, parents, patient, mode = "2samples_extend", swap = TRUE)

cat("\n-------------------------------\n")
patient <- "OESO_0093"
cat(patient)
samples <- c("TM", "TU")
frac.table <- matrix(c(100, 100, # C1
                       0, 55, # C2
                       50, 40, # C4
                       35, 0 # C3
), ncol = length(samples), byrow = TRUE)
colnames(frac.table) <- samples
## clusterID 1  2  4  3 
parents <- c(0, 1, 1, 3)
.plotFish(frac.table, parents, patient, mode = "2samples_extend")
.plotFish(frac.table, parents, patient, mode = "2samples_extend", swap = TRUE)

cat("\n-------------------------------\n")
patient <- "OESO_0096_preOnly"
samples <- c("TM", "TU")
frac.table <- matrix(c(100, 100, # C1
                       15, 50, # C2
                       2, 45, # C5
                       0, 30, # C3
                       65, 45, # C6
                       55, 0 # C4
), ncol = length(samples), byrow = TRUE)
colnames(frac.table) <- samples
## clusterID 1  2  5  3  6  4
parents <- c(0, 1, 2, 3, 1, 5)
.plotFish(frac.table, parents, patient, mode = "2samples_extend")
.plotFish(frac.table, parents, patient, mode = "2samples_extend", swap = TRUE)

cat("\n-------------------------------\n")
patient <- "OESO_0109"
samples <- c("TL", "TM")
frac.table <- matrix(c(100, 100, # C1
                       65, 45, # C2
                       60, 0, # C3
                       0, 50 # C4
), ncol = length(samples), byrow = TRUE)
colnames(frac.table) <- samples
## clusterID 1  2  3  4
parents <- c(0, 1, 2, 1)
.plotFish(frac.table, parents, patient, mode = "2samples_extend")
.plotFish(frac.table, parents, patient, mode = "2samples_extend", swap = TRUE)

cat("\n-------------------------------\n")
patient <- "OESO_0113_preOnly"
samples <- c("TM", "TU")
frac.table <- matrix(c(100, 100, # C1
                       65, 95, # C4
                       25, 20, # C2
                       2, 70, # C3
                       0, 35, # C5
                       35, 0 # C6
), ncol = length(samples), byrow = TRUE)
colnames(frac.table) <- samples
## clusterID 1  4  2  3  5  6
parents <- c(0, 1, 2, 2, 4, 2)
.plotFish(frac.table, parents, patient, mode = "2samples_extend")
.plotFish(frac.table, parents, patient, mode = "2samples_extend", swap = TRUE)

cat("\n-------------------------------\n")
patient <- "OESO_0117"
samples <- c("TL", "TU")
frac.table <- matrix(c(100, 100, # C1
                       25, 30, # C2
                       0, 65, # C3
                       70, 0 # C4
), ncol = length(samples), byrow = TRUE)
colnames(frac.table) <- samples
## clusterID 1  2  3  4
parents <- c(0, 1, 1, 1)
.plotFish(frac.table, parents, patient, mode = "2samples_extend")
.plotFish(frac.table, parents, patient, mode = "2samples_extend", swap = TRUE)

cat("\n-------------------------------\n")
patient <- "OESO_0119"
samples <- c("TL", "TU")
frac.table <- matrix(c(100, 100, # C1
                       40, 95, # C2
                       35, 40, # C4
                       0, 50, # C5
                       55, 0 # C3
), ncol = length(samples), byrow = TRUE)
colnames(frac.table) <- samples
## clusterID 1  2  4  5  3
parents <- c(0, 1, 2, 2, 1)
.plotFish(frac.table, parents, patient, mode = "2samples_extend")
.plotFish(frac.table, parents, patient, mode = "2samples_extend", swap = TRUE)

cat("\n-------------------------------\n")
patient <- "OESO_0121"
samples <- c("TM", "TU")
frac.table <- matrix(c(100, 100, # C1
                       90, 80, # C7
                       45, 5, # C3
                       20, 0, # C5
                       40, 50, # C4
                       20, 35, # C2
                       0, 30 # C6
), ncol = length(samples), byrow = TRUE)
colnames(frac.table) <- samples
## clusterID 1  7  3  5  4  2  6
parents <- c(0, 1, 2, 3, 2, 5, 6)
.plotFish(frac.table, parents, patient, mode = "2samples_extend")
.plotFish(frac.table, parents, patient, mode = "2samples_extend", swap = TRUE)

cat("\n-------------------------------\n")
patient <- "OESO_6129"
samples <- c("TM", "TU")
frac.table <- matrix(c(100, 100, # C1
                       40, 95, # C5
                       15, 25, # C7
                       0, 65, # C2
                       0, 25, # C3
                       55, 0, # C4
                       25, 0 # C6
), ncol = length(samples), byrow = TRUE)
colnames(frac.table) <- samples
## clusterID 1  5  7  2  3  4  6
parents <- c(0, 1, 2, 2, 4, 1, 6)
.plotFish(frac.table, parents, patient, mode = "2samples_extend")
.plotFish(frac.table, parents, patient, mode = "2samples_extend", swap = TRUE)

cat("\n-------------------------------\n")
patient <- "OESO_0098"
samples <- c("TL", "TM")
frac.table <- matrix(c(100, 100, # C1
                       25, 20, # C4
                       70, 2, # C3
                       65, 0, # C2
                       30, 0, # C6
                       0, 65 # C5
), ncol = length(samples), byrow = TRUE)
colnames(frac.table) <- samples
## clusterID 1  4  3  2  6  5
parents <- c(0, 1, 1, 3, 4, 1)
.plotFish(frac.table, parents, patient, mode = "2samples_extend")
.plotFish(frac.table, parents, patient, mode = "2samples_extend", swap = TRUE)

cat("\n-------------------------------\n")
patient <- "OESO_0118"
samples <- c("TL", "TM")
frac.table <- matrix(c(100, 100, # C1
                       95, 45, # C5
                       90, 20, # C4
                       65, 10, # C2
                       0, 50 # C3
), ncol = length(samples), byrow = TRUE)
colnames(frac.table) <- samples
## clusterID 1  5  4  2  3
parents <- c(0, 1, 2, 3, 1)
.plotFish(frac.table, parents, patient, mode = "2samples_extend")
.plotFish(frac.table, parents, patient, mode = "2samples_extend", swap = TRUE)

cat("\n-------------------------------\n")
patient <- "OESO_0120"
samples <- c("TL", "TM")
frac.table <- matrix(c(100, 100, # C1
                       65, 95, # C4
                       20, 20, # C5
                       0, 65, # C3
                       0, 30, # C6
                       40, 5 # C2
), ncol = length(samples), byrow = TRUE)
colnames(frac.table) <- samples
## clusterID 1  4  5  3  6  2
parents <- c(0, 1, 2, 2, 4, 2)
.plotFish(frac.table, parents, patient, mode = "2samples_extend")
.plotFish(frac.table, parents, patient, mode = "2samples_extend", swap = TRUE)

cat("\n-------------------------------\n")
patient <- "OESO_0125_preOnly"
samples <- c("TL", "TU")
frac.table <- matrix(c(100, 100, # C1
                       2, 30, # C2
                       30, 0, # C3
                       20, 25 # C4
), ncol = length(samples), byrow = TRUE)
colnames(frac.table) <- samples
## clusterID 1  4  2  3
parents <- c(0, 1, 1, 1)
.plotFish(frac.table, parents, patient, mode = "2samples_extend")
.plotFish(frac.table, parents, patient, mode = "2samples_extend", swap = TRUE)

cat("\n-------------------------------\n")
patient <- "OESO_0227"
samples <- c("TM", "TU")
frac.table <- matrix(c(100, 100, # C1
                       30, 25, # C4
                       0, 70, # C2
                       0, 30, # C3
                       65, 0 # C5
), ncol = length(samples), byrow = TRUE)
colnames(frac.table) <- samples
## clusterID 1  4  2  3  5
parents <- c(0, 1, 1, 3, 1)
.plotFish(frac.table, parents, patient, mode = "2samples_extend")
.plotFish(frac.table, parents, patient, mode = "2samples_extend", swap = TRUE)

cat("\n-------------------------------\n")
patient <- "OESO_0009_preOnly"
samples <- c("TL", "TU")
frac.table <- matrix(c(100, 100, # C1
                       0, 90, # C3
                       0, 50, # C2
                       70, 0, # C4
                       25, 0 # C5
), ncol = length(samples), byrow = TRUE)
colnames(frac.table) <- samples
## clusterID 1  3  2  4  5
parents <- c(0, 1, 2, 1, 4)
.plotFish(frac.table, parents, patient, mode = "2samples_extend")
.plotFish(frac.table, parents, patient, mode = "2samples_extend", swap = TRUE)





