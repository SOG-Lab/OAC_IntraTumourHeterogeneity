cat("Project ID: OAC_Heterogeneity_MR_MT\n")
## R script to plot patient survival and recurrence summary (swimmer plot)##
cat("Loading: R script to plot patient survival and recurrence summary (swimmer plot)\n")

## Loading external libraries ##
library(ggplot2)
library(ggnewscale)
library(tidyr)

## Loading external scripts ##
source("SCRIPTS/general.R")

## Preparation ##
out.dir <- "RESULTS/Cohort/"

plot.swimmerPlot.survival <- function(mode = "endoscopy"){
  cohort <- general.load.data("DATA/Annotation/Cohort.csv")
  
  ## Remove first consultation date (endoscopy will be the start date)
  cohort$Date.First.Consultation <- NULL
  cohort$Baseline <- cohort$Date.Endoscopy
  
  ## Transform into long format
  cohort_long <- as.data.frame(pivot_longer(data = cohort, names_to = "Event", values_to = "EventDate", 
                                            cols = c("Date.Endoscopy", "Surgery.Date", "DOCTOR.Recurrence", "DOCTOR.Last.Followup")))
  ## Remove recurrence for Alive with NSR and Dead without disease patients
  cohort_long <- cohort_long[!(cohort_long$Event == "DOCTOR.Recurrence" & cohort_long$Current.Status == "Alive with NSR"), ]
  cohort_long <- cohort_long[!(cohort_long$Event == "DOCTOR.Recurrence" & cohort_long$Current.Status == "Dead without disease"), ]
  
  ## Rename events
  cohort_long$Event <- gsub("[.]", "", gsub("Date", "", gsub("Date.Endoscopy", "Baseline", gsub("DOCTOR", "", cohort_long$Event))))
  
  ## Add sample timepoints
  cohort_long$EventDate <- as.Date(cohort_long$EventDate, format = "%d/%m/%Y")
  cohort_long$Baseline <- as.Date(cohort_long$Baseline, format = "%d/%m/%Y")
  cohort_long$timeDays <- as.numeric(cohort_long$EventDate - cohort_long$Baseline)
  cohort_long$time <- cohort_long$timeDays/30.4167
  
  ## Amend survival information
  cohort_long$Status <- cohort_long$Current.Status
  cohort_long$Status[grepl("Dead", cohort_long$Status)] <- "Deceased" 
  
  ## Cut survival at 60 months after surgery for each patient
  if(mode == "surgery"){
    a <- sapply(cohort$Patient, function(patient){
      cohort_long$time[cohort_long$Event == "LastFollowup" & cohort_long$Patient == patient] <<- cohort_long$time[cohort_long$Event == "Surgery" & cohort_long$Patient == patient] + cohort_long$OS.surgery[cohort_long$Event == "Surgery" & cohort_long$Patient == patient]
    })
    insert <- "_surgery"
  }else if(mode == "endoscopy"){
    cohort_long$OS[cohort_long$OS == 1 & cohort_long$OS.endo > 60] <- 0
    cohort_long$Status[cohort_long$OS == 1 & cohort_long$OS.endo > 60] <- "Alive with NSR"
    cohort_long$OS.endo[cohort_long$OS.endo > 60] <- 60
    cohort_long$time[cohort_long$time > 60] <- 60
    insert <- "_endoscopy"
  }else{
    insert <- "_noSurvivalCut"
  }
  
  ## Sort donors by time
  cohort_long <- cohort_long[order(cohort_long$time, decreasing = TRUE),]
  cohort_long$Patient <- factor(cohort_long$Patient, levels = unique(cohort_long$Patient))
  
  ## Category labels
  cohort_long$Event[cohort_long$Event == "Baseline"] <- "Endoscopy"
  survivalEvents <- c("Alive with NSR", "Alive with disease", "Deceased")
  clinicalEvents <- c("Endoscopy", "Recurrence", "Surgery")
  events <- c(survivalEvents, clinicalEvents)
  survivalEvents <- factor(survivalEvents, levels = survivalEvents)
  clinicalEvents <- factor(clinicalEvents, levels = clinicalEvents)
  
  ## Colours
  eventColours <- c("lightgrey", "darkgrey", "black", "#4DAF4A", "#E41A1C", "blue")
  names(eventColours) <- events
  
  ## Plot
  ggplot(cohort_long, aes(x = time, y = Patient, group = Patient)) +
    geom_line(size = 2, aes(color = Status)) +
    geom_point(data = cohort_long[cohort_long$Event %in% clinicalEvents,], aes(x = time, y = Patient, color = Event, shape = Event), size = 4) +
    labs(x = "Time [months]", y = "Patient") +
    scale_x_continuous(breaks = seq(0, max(cohort_long$time, na.rm = TRUE), by = 12)) +
    scale_color_manual(name = "Status", 
                       labels = events,
                       values = eventColours) +
    scale_fill_manual(name = "Event",
                      labels = clinicalEvents,
                      values = eventColours) +    
    theme_test()
  ggsave(paste0(out.dir, "patientOverview_months", insert, ".png"), width = 17, height = 10, type = "cairo")
  ggsave(paste0(out.dir, "patientOverview_months", insert, ".pdf"), width = 17, height = 10)
  general.save.data(cohort_long, paste0(out.dir, "Plot_data_long", insert))
  cohort_long$EventDate <- cohort_long$timeDays <- NULL
  cohort_wide <- as.data.frame(pivot_wider(cohort_long, names_from = Event, values_from = time))
  general.save.data(cohort_wide, paste0(out.dir, "Plot_data_wide", insert))
}


