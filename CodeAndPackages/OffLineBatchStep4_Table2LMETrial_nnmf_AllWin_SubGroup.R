library(lmerTest)
library(dplyr)

## ------------------ Base path ------------------ ##
base_dir <- "//nimhlabstore1.nimh.nih.gov/UFNC/FNC2/Zhang/Projects/Project-LocalProcessing/Step4/awakeRefSpon/08-Dec-2025/NonNegMatFac"

## only deltaF, but two score versions encoded in file names
## e.g. deltaFAllSLMScoreTrialDynWin1.csv
##      deltaFNonTargetSLMScoreTrialDynWin1.csv
all_files <- list.files(
  base_dir,
  pattern = "^deltaF.*SLMScoreTrialDynWin[0-9]+\\.csv$",
  full.names = TRUE
)

if (length(all_files) == 0) {
  stop("No deltaF*SLMScoreTrialDynWinX.csv files found in base_dir.")
}

cat("Found files:\n", paste(basename(all_files), collapse = "\n"), "\n\n")

## response variables to model
response_vars <- c("SpeedScore",
                   "SensoryScore",
                   "SpeedScore_SpeedReg",
                   "SensoryScore_SpeedReg")

## container for fixed-effect summaries
all_effects <- list()

## helper: extract fixed effects into tidy rows
extract_fixed_effects <- function(model,
                                  window,
                                  condition,    # "NonSensory"/"Sensory"
                                  subset_label, # e.g. "NonSensory_G1"
                                  model_name,   # "Base"/"Extended"
                                  response_var,
                                  cell_subset)  # "All" or "NonTarget"
{
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
  
  df0$DataType    <- "deltaF"
  df0$Window      <- window
  df0$Condition   <- condition
  df0$Subset      <- subset_label
  df0$Model       <- model_name
  df0$ResponseVar <- response_var
  df0$CellSubset  <- cell_subset   # All vs NonTarget
  
  df0
}

## ------------------ Main loop over files ------------------ ##
for (csv_path in all_files) {
  
  fname <- basename(csv_path)
  
  ## window: Win1, Win2, ...
  window <- gsub(".*(Win[0-9]+)\\.csv$", "\\1", fname)
  
  ## identify cell subset from file name
  ##   "NonTarget" -> CellSubset = "NonTarget"
  ##   "All"       -> CellSubset = "All"
  cell_subset <- if (grepl("NonTarget", fname, ignore.case = TRUE)) {
    "NonTarget"
  } else if (grepl("All", fname, ignore.case = TRUE)) {
    "All"
  } else {
    "Unknown"
  }
  
  cat("Processing file:", fname,
      "| Window:", window,
      "| CellSubset:", cell_subset, "\n")
  
  dataTable <- read.csv(csv_path)
  
  ## factors
  dataTable$Group   <- factor(dataTable$Group)
  dataTable$Cell    <- factor(dataTable$Cell)
  dataTable$Sensory <- factor(dataTable$Sensory)
  dataTable$Session <- factor(dataTable$Session)
  
  ## ---------- z-score selected predictors ----------
  scale_vars <- c("Speed",
                  "TargetCellN",
                  "TargetSpeedR",
                  "TargetSensoryR",
                  "AveDist",
                  "MinDist")
  
  for (v in scale_vars) {
    if (v %in% names(dataTable)) {
      dataTable[[paste0(v, "_z")]] <- as.numeric(scale(dataTable[[v]]))
    }
  }
  
  ## sanity check: response columns
  missing_resp <- setdiff(response_vars, names(dataTable))
  if (length(missing_resp) > 0) {
    warning("File ", fname,
            " is missing response columns: ",
            paste(missing_resp, collapse = ", "),
            " – those will be skipped.")
  }
  
  ## -------------------------------------------------------
  ## Loop over Sensory (0/1) and Group (ALL, 1, 2, 3)
  ## -------------------------------------------------------
  for (sens_val in c(0, 1)) {
    cond_label <- ifelse(sens_val == 0, "NonSensory", "Sensory")
    
    for (grp_val in c(NA, 1, 2, 3)) {
      
      if (is.na(grp_val)) {
        subset_label <- paste0(cond_label, "_All")
        file_tag     <- subset_label
        dat_sub <- subset(dataTable, Sensory == sens_val)
      } else {
        subset_label <- paste0(cond_label, "_G", grp_val)
        file_tag     <- subset_label
        dat_sub <- subset(dataTable,
                          Sensory == sens_val &
                            Group == grp_val)
      }
      
      if (nrow(dat_sub) <= 10) {
        message("  Not enough rows for ",
                fname, " subset ", file_tag,
                " (n = ", nrow(dat_sub), ") – skipping this subset.")
        next
      }
      
      ## ---------- loop over response variables ---------- ##
      for (resp in response_vars) {
        
        if (!resp %in% names(dat_sub)) {
          message("  Response ", resp, " not found in file ",
                  fname, " subset ", file_tag, " – skipping.")
          next
        }
        
        cat("    Response:", resp, "| subset:", file_tag, "\n")
        
        ## BASE model:
        ##   Resp ~ (Speed_z + TargetCellN_z) * (TargetSpeedR_z + TargetSensoryR_z) + (1|Session)
        base_formula <- as.formula(
          paste0(resp,
                 " ~ (Speed_z + TargetCellN_z) * (TargetSpeedR_z + TargetSensoryR_z) + (1|Session)")
        )
        
        ## EXTENDED model: + AveDist_z + MinDist_z
        ext_formula <- as.formula(
          paste0(resp,
                 " ~ (Speed_z + TargetCellN_z) * (TargetSpeedR_z + TargetSensoryR_z) +
                    AveDist_z + MinDist_z + (1|Session)")
        )
        
        l_base <- try(lmer(base_formula, data = dat_sub, REML = FALSE),
                      silent = TRUE)
        if (inherits(l_base, "try-error")) {
          message("      Base model failed for ", resp,
                  " subset ", file_tag, " – skipping.")
          next
        }
        
        l_ext <- try(lmer(ext_formula, data = dat_sub, REML = FALSE),
                     silent = TRUE)
        if (inherits(l_ext, "try-error")) {
          message("      Extended model failed for ", resp,
                  " subset ", file_tag, " – skipping extended.")
          l_ext <- NULL
        }
        
        ## ---- write text summary ---- ##
        out_file <- file.path(
          base_dir,
          paste0("deltaF_", window, "_", cell_subset, "_",
                 resp, "_TrialBased_", file_tag, "_LME.txt")
        )
        sink(out_file)
        cat("### File:", fname,
            "| Window:", window,
            "| CellSubset:", cell_subset,
            "| Response:", resp,
            "| Subset:", file_tag, "\n\n")
        
        cat("---- Base model ----\n")
        print(summary(l_base))
        
        if (!is.null(l_ext)) {
          cat("\n\n---- Extended model ----\n")
          print(summary(l_ext))
          cat("\n\n---- Likelihood ratio test (Extended vs Base) ----\n")
          print(anova(l_base, l_ext))
        } else {
          cat("\n\nExtended model did not converge / failed.\n")
        }
        sink()
        
        ## ---- store fixed effects ---- ##
        all_effects[[length(all_effects) + 1]] <-
          extract_fixed_effects(l_base,
                                window,
                                cond_label,
                                subset_label,
                                "Base",
                                resp,
                                cell_subset)
        
        if (!is.null(l_ext)) {
          all_effects[[length(all_effects) + 1]] <-
            extract_fixed_effects(l_ext,
                                  window,
                                  cond_label,
                                  subset_label,
                                  "Extended",
                                  resp,
                                  cell_subset)
        }
      } # end resp loop
    }   # end group loop
  }     # end Sensory loop
}       # end file loop

## ------------------ Combine & save summary ------------------ ##
effects_df <- do.call(rbind, all_effects)

effects_df <- effects_df %>%
  select(DataType, Window, CellSubset, Condition, Subset,
         Model, ResponseVar, Term, Estimate, StdError, tValue, pValue)

out_csv <- file.path(base_dir,
                     "LME_FixedEffectsSummary_SpeedSensory_ZSCALED_AllVsNonTarget_AllWindows.csv")
write.csv(effects_df, out_csv, row.names = FALSE)

cat("Saved fixed-effect summary to:\n", out_csv, "\n")
