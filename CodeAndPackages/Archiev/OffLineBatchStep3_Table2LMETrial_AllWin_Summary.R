## ------------------ Libraries ------------------ ##
library(lmerTest)
library(dplyr)
library(ggplot2)

## ------------------ Base paths ------------------ ##
base_dir <- "//nimhlabstore1.nimh.nih.gov/UFNC/FNC2/Zhang/Projects/Project-LocalProcessing/Step3/awakeRefSpon/GroupSLM20-Nov-2025"

ValidSG <- read.csv(file.path(base_dir, "ValidFOVGroup.csv"))
valid_keys <- unique(ValidSG[c("Session", "Group")])

data_types <- c("deltaF", "spks")

## Detect all available WinX CSVs once
all_files <- list.files(base_dir, pattern = "Win[0-9]+\\.csv$", full.names = TRUE)
win_list  <- unique(gsub(".*Win([0-9]+)\\.csv$", "Win\\1",
                         basename(all_files)))
cat("Detected windows:", paste(win_list, collapse = ", "), "\n")

## Container to store all fixed-effect summaries
all_effects <- list()

## --------------- Helper: extract fixed effects --------------- ##
extract_fixed_effects <- function(model,
                                  data_type,
                                  window,
                                  condition,
                                  model_name) {
  if (is.null(model)) return(NULL)
  
  cf  <- coef(summary(model))
  df0 <- data.frame(
    Term      = rownames(cf),
    Estimate  = cf[, "Estimate"],
    StdError  = cf[, "Std. Error"],
    tValue    = cf[, "t value"],
    pValue    = cf[, "Pr(>|t|)"],
    stringsAsFactors = FALSE
  )
  
  df0$DataType  <- data_type
  df0$Window    <- window
  df0$Condition <- condition   # "NonSensory" or "Sensory"
  df0$Model     <- model_name  # "Base" or "Extended"
  
  # Numeric window index (e.g., Win1 -> 1)
  df0$WindowNum <- as.numeric(gsub("Win", "", window))
  
  df0
}

## ------------------ Main loop ------------------ ##
for (data_type in data_types) {
  for (win in win_list) {
    
    csv_path <- file.path(base_dir,
                          paste0(data_type, "SLMGroupResTrialDyn", win, ".csv"))
    
    if (!file.exists(csv_path)) {
      message("Skipping missing file: ", csv_path)
      next
    }
    
    cat("Processing:", data_type, win, "\n")
    
    dataTable <- read.csv(csv_path)
    dataTable <- merge(dataTable, valid_keys,
                       by = c("Session", "Group"),
                       all = FALSE)
    
    dataTable$Group   <- factor(dataTable$Group)
    dataTable$Cell    <- factor(dataTable$Cell)
    dataTable$Sensory <- factor(dataTable$Sensory)
    dataTable$Session <- factor(dataTable$Session)
    
    ## ---------- Non-sensory trials (Sensory == 0) ---------- ##
    dat0 <- subset(dataTable, Sensory == 0 & NonTargetCell == 1)
    l1 <- l2 <- NULL
    
    if (nrow(dat0) > 10) {
      l1 <- lmer(Response ~ (SpeedR + SensoryR + TargetCellN) *
                   (TargetSpeedR + TargetSensoryR) +
                   SpeedR * Speed + (1|Session),
                 data = dat0, REML = FALSE)
      
      l2 <- lmer(Response ~ (SpeedR + SensoryR + TargetCellN) *
                   (TargetSpeedR + TargetSensoryR) +
                   SpeedR * Speed + AveDist + MinDist + (1|Session),
                 data = dat0, REML = FALSE)
      
      sink(file.path(base_dir,
                     paste0(data_type, "_", win, "_TrialBased_NonSensoryLME.txt")))
      print(summary(l1))
      print(summary(l2))
      print(anova(l1, l2))
      sink()
      
      all_effects[[length(all_effects) + 1]] <-
        extract_fixed_effects(l1, data_type, win, "NonSensory", "Base")
      all_effects[[length(all_effects) + 1]] <-
        extract_fixed_effects(l2, data_type, win, "NonSensory", "Extended")
    }
    
    ## ---------- Sensory trials (Sensory == 1) ---------- ##
    dat1 <- subset(dataTable, Sensory == 1 & NonTargetCell == 1)
    l3 <- l4 <- NULL
    
    if (nrow(dat1) > 10) {
      l3 <- lmer(Response ~ (SpeedR + SensoryR + TargetCellN) *
                   (TargetSpeedR + TargetSensoryR) +
                   SpeedR * Speed + (1|Session),
                 data = dat1, REML = FALSE)
      
      l4 <- lmer(Response ~ (SpeedR + SensoryR + TargetCellN) *
                   (TargetSpeedR + TargetSensoryR) +
                   SpeedR * Speed + AveDist + MinDist + (1|Session),
                 data = dat1, REML = FALSE)
      
      sink(file.path(base_dir,
                     paste0(data_type, "_", win, "_TrialBased_SensoryLME.txt")))
      print(summary(l3))
      print(summary(l4))
      print(anova(l3, l4))
      sink()
      
      all_effects[[length(all_effects) + 1]] <-
        extract_fixed_effects(l3, data_type, win, "Sensory", "Base")
      all_effects[[length(all_effects) + 1]] <-
        extract_fixed_effects(l4, data_type, win, "Sensory", "Extended")
    }
  }
}

## ------------------ Combine & Save Summary ------------------ ##
effects_df <- do.call(rbind, all_effects)

# Reorder columns nicely
effects_df <- effects_df %>%
  select(DataType, Window, WindowNum, Condition, Model,
         Term, Estimate, StdError, tValue, pValue)

write.csv(effects_df,
          file.path(base_dir, "LME_FixedEffectsSummary_AllWindows.csv"),
          row.names = FALSE)

cat("Saved fixed-effect summary to LME_FixedEffectsSummary_AllWindows.csv\n")

## ------------------ Plotting across windows ------------------ ##

}
