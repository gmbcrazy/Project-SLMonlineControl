library(lmerTest)
library(dplyr)

# Base directory
base_dir <- "//nimhlabstore1.nimh.nih.gov/UFNC/FNC2/Zhang/Projects/Project-LocalProcessing/Step3/awakeRefSpon/GroupSLM20-Nov-2025"

# Read in valid session/group list
ValidSG <- read.csv(file.path(base_dir, "ValidFOVGroup.csv"))
valid_keys <- unique(ValidSG[c("Session", "Group")])

# Data types to process
data_types <- c("deltaF", "spks")

# Detect available WinX CSVs automatically
all_files <- list.files(base_dir, pattern = "Win[0-9]+\\.csv$", full.names = TRUE)

# Extract unique window identifiers (e.g., Win1, Win2, Win3)
win_list <- unique(gsub(".*Win([0-9]+)\\.csv$", "Win\\1",
                        basename(all_files)))

cat("Detected windows:", paste(win_list, collapse=", "), "\n")

# ---------------------------------------------------------
# Loop through data types and windows
# ---------------------------------------------------------
for (data_type in data_types) {
  
  for (win in win_list) {
    
    # Build filename
    csv_path <- file.path(base_dir,
                          paste0(data_type, "SLMGroupResTrialDyn", win, ".csv"))
    
    if (!file.exists(csv_path)) {
      message("Skipping missing file: ", csv_path)
      next
    }
    
    cat("Processing:", csv_path, "\n")
    
    # Read trial data
    dataTable <- read.csv(csv_path)
    
    # Keep only valid session–group pairs
    dataTable <- merge(dataTable, valid_keys,
                       by = c("Session", "Group"),
                       all = FALSE)
    
    # Convert to factors
    dataTable$Group    <- factor(dataTable$Group)
    dataTable$Cell     <- factor(dataTable$Cell)
    dataTable$Sensory  <- factor(dataTable$Sensory)
    dataTable$Session  <- factor(dataTable$Session)
    
    # ------------------------------
    # Build output filenames
    # ------------------------------
    out_nonSens <- file.path(base_dir,
                             paste0(data_type, "_", win, "_NonSensoryLME.txt"))
    out_sens <- file.path(base_dir,
                          paste0(data_type, "_", win, "_SensoryLME.txt"))
    
    # ---------------------------------------------------------
    # Model set 1: Sensory == 0
    # ---------------------------------------------------------
    dat0 <- subset(dataTable, Sensory == 0 & NonTargetCell == 1)
    
    if (nrow(dat0) > 10) {
      l1 <- lmer(Response ~ (SpeedR + SensoryR + TargetCellN) *
                   (TargetSpeedR + TargetSensoryR) + 
                   SpeedR * Speed + (1|Session),
                 data = dat0, REML = FALSE)
      
      l2 <- lmer(Response ~ (SpeedR + SensoryR + TargetCellN) *
                   (TargetSpeedR + TargetSensoryR) +
                   SpeedR * Speed + AveDist + MinDist + (1|Session),
                 data = dat0, REML = FALSE)
      
      sink(out_nonSens)
      print(summary(l1))
      print(summary(l2))
      print(anova(l1, l2))
      sink()
    }
    
    # ---------------------------------------------------------
    # Model set 2: Sensory == 1
    # ---------------------------------------------------------
    dat1 <- subset(dataTable, Sensory == 1 & NonTargetCell == 1)
    
    if (nrow(dat1) > 10) {
      l3 <- lmer(Response ~ (SpeedR + SensoryR + TargetCellN) *
                   (TargetSpeedR + TargetSensoryR) + 
                   SpeedR * Speed + (1|Session),
                 data = dat1, REML = FALSE)
      
      l4 <- lmer(Response ~ (SpeedR + SensoryR + TargetCellN) *
                   (TargetSpeedR + TargetSensoryR) +
                   SpeedR * Speed + AveDist + MinDist + (1|Session),
                 data = dat1, REML = FALSE)
      
      sink(out_sens)
      print(summary(l3))
      print(summary(l4))
      print(anova(l3, l4))
      sink()
    }
  }
}
