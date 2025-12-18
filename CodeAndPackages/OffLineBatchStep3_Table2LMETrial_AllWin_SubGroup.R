library(lmerTest)
library(dplyr)

## ------------------ Base paths ------------------ ##
base_dir <- "//nimhlabstore1.nimh.nih.gov/UFNC/FNC2/Zhang/Projects/Project-LocalProcessing/Step3/awakeRefSpon/GroupSLM20-Nov-2025"

ValidSG <- read.csv(file.path(base_dir, "ValidFOVGroup.csv"))
valid_keys <- unique(ValidSG[c("Session", "Group")])

data_types <- c("deltaF", "spks")

## Detect all available WinX CSVs
all_files <- list.files(base_dir, pattern = "Win[0-9]+\\.csv$", full.names = TRUE)
win_list  <- unique(gsub(".*Win([0-9]+)\\.csv$", "Win\\1",
                         basename(all_files)))
cat("Detected windows:", paste(win_list, collapse = ", "), "\n")

## Container to store fixed-effect summaries
all_effects <- list()

## ---------- Helper: extract fixed effects into a tidy row set ---------- ##
extract_fixed_effects <- function(model,
                                  data_type,
                                  window,
                                  condition,    # "NonSensory"/"Sensory"
                                  subset_label, # e.g. "NonSensory_G1"
                                  model_name) { # "Base"/"Extended"
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
  df0$Condition <- condition
  df0$Subset    <- subset_label
  df0$Model     <- model_name
  
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
    
    ## ------------------ Factor conversions ------------------ ##
    dataTable$Group   <- factor(dataTable$Group)
    dataTable$Cell    <- factor(dataTable$Cell)
    dataTable$Sensory <- factor(dataTable$Sensory)
    dataTable$Session <- factor(dataTable$Session)
    
    ## ------------------ Z-score continuous predictors ------------------ ##
    ## Define which columns to scale (modify this vector if needed)
    scale_vars <- c("SpeedR",
                    "SensoryR",
                    "TargetCellN",
                    "TargetSpeedR",
                    "TargetSensoryR",
                    "Speed",
                    "AveDist",
                    "MinDist")
    
    for (v in scale_vars) {
      if (v %in% names(dataTable)) {
        new_name <- paste0(v, "_z")
        dataTable[[new_name]] <- as.numeric(scale(dataTable[[v]]))
      } else {
        warning("Variable ", v, " not found in dataTable for ",
                data_type, " ", win)
      }
    }
    
    ## -------------------------------------------------------
    ## Loop over Sensory (0/1) and Group (ALL, 1, 2, 3)
    ## -------------------------------------------------------
    for (sens_val in c(0, 1)) {
      cond_label <- ifelse(sens_val == 0, "NonSensory", "Sensory")
      
      for (grp_val in c(NA, 1, 2, 3)) {
        
        if (is.na(grp_val)) {
          subset_label <- paste0(cond_label, "_All")  # e.g. "NonSensory_All"
          file_tag     <- subset_label
          dat_sub <- subset(dataTable,
                            Sensory == sens_val & NonTargetCell == 1)
        } else {
          subset_label <- paste0(cond_label, "_G", grp_val)  # e.g. "NonSensory_G1"
          file_tag     <- subset_label
          dat_sub <- subset(dataTable,
                            Sensory == sens_val &
                            NonTargetCell == 1 &
                            Group == grp_val)
        }
        
        if (nrow(dat_sub) <= 10) {
          message("Not enough rows for ",
                  data_type, " ", win, " subset ", file_tag,
                  " (n = ", nrow(dat_sub), ") – skipping.")
          next
        }
        
        ## ---------------- Fit models with *scaled* predictors ---------------- ##
        ## NOTE: We now use the *_z versions of continuous covariates
        
        l_base <- lmer(Response ~ (SpeedR_z + SensoryR_z + TargetCellN_z) *
                         (TargetSpeedR_z + TargetSensoryR_z) +
                         (SpeedR_z+TargetSpeedR_z)* Speed_z + (1 | Session),
                       data = dat_sub,
                       REML = FALSE)
        
        l_ext  <- lmer(Response ~ (SpeedR_z + SensoryR_z + TargetCellN_z) *
                         (TargetSpeedR_z + TargetSensoryR_z) +
                         (SpeedR_z+TargetSpeedR_z) * Speed_z + AveDist_z + MinDist_z +
                         (1 | Session),
                       data = dat_sub,
                       REML = FALSE)
        
        ## ---- Save text summary for this (window, condition, group) ---- ##
        out_file <- file.path(
          base_dir,
          paste0(data_type, "_", win, "_TrialBased_", file_tag, "_LME_scaled.txt")
        )
        sink(out_file)
        cat("### SCALED MODEL ###\n")
        cat("Data type:", data_type,
            "| Window:", win,
            "| Subset:", file_tag, "\n\n")
        
        cat("---- Base model (scaled predictors) ----\n")
        print(summary(l_base))
        cat("\n\n---- Extended model (scaled predictors) ----\n")
        print(summary(l_ext))
        cat("\n\n---- Likelihood ratio test (Extended vs Base) ----\n")
        print(anova(l_base, l_ext))
        sink()
        
        ## ---- Store fixed effects for later comparison ---- ##
        all_effects[[length(all_effects) + 1]] <-
          extract_fixed_effects(l_base, data_type, win,
                                cond_label, subset_label, "Base_scaled")
        all_effects[[length(all_effects) + 1]] <-
          extract_fixed_effects(l_ext, data_type, win,
                                cond_label, subset_label, "Extended_scaled")
      }
    }
  }
}

## ------------------ Combine & Save Summary ------------------ ##
effects_df <- do.call(rbind, all_effects)

effects_df <- effects_df %>%
  select(DataType, Window, Condition, Subset, Model,
         Term, Estimate, StdError, tValue, pValue)

write.csv(effects_df,
          file.path(base_dir,
                    "LME_FixedEffectsSummary_AllWindows_Groups_SCALED.csv"),
          row.names = FALSE)

cat("Saved fixed-effect summary to LME_FixedEffectsSummary_AllWindows_Groups_SCALED1.csv\n")
