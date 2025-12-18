#library(R.matlab)
library(lmerTest)
library(dplyr)


#library(lme4)

# Define base directory
base_dir <- "//nimhlabstore1.nimh.nih.gov/UFNC/FNC2/Zhang/Projects/Project-LocalProcessing/Step3/awakeRefSpon/GroupSLM20-Nov-2025"

ValidSG <- read.csv(file.path(base_dir, "ValidFOVGroup.csv"))
valid_keys <- unique(ValidSG[c("Session", "Group")])

# List subdirectories (assuming numeric names only)
subdirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)

data_types <- c("deltaF", "spks") 

for (data_type in data_types) {
  csv_path <- file.path(base_dir, paste0(data_type,"SLMGroupResTrialDynWin1.csv"))
  

#csv_path <- file.path(base_dir, "deltaFSLMGroupResTrialDynWin1.csv")
  
# Check if file exists before proceeding
if (!file.exists(csv_path)) {
   message(paste("File not found:", csv_path))
   next
}
  
# Read data
dataTable <- read.csv(csv_path)
dataTable <- merge(dataTable,valid_keys, by = c("Session", "Group"),all = FALSE)          # inner join

  
# Convert to factor
dataTable$Group <- factor(dataTable$Group)
dataTable$Cell <- factor(dataTable$Cell)
dataTable$Sensory <- factor(dataTable$Sensory)
dataTable$Session <- factor(dataTable$Session)

# Mixed linear model for Sensory == 0
l1 <- lmer(Response ~ (SpeedR+SensoryR+TargetCellN) * (TargetSpeedR+TargetSensoryR) + SpeedR*Speed + (1|Session),data = subset(dataTable, Sensory == 0 & NonTargetCell == 1),REML = FALSE)
l2 <- lmer(Response ~ (SpeedR+SensoryR+TargetCellN) * (TargetSpeedR+TargetSensoryR) + SpeedR*Speed + AveDist + MinDist+(1|Session),data = subset(dataTable, Sensory == 0 & NonTargetCell == 1),REML = FALSE)


sink(file.path(base_dir, paste0(data_type,"_TrialBased_NonSensoryLME.txt")), append = FALSE)
print(summary(l1))
print(summary(l2))
print(anova(l1,l2))
sink()

  # Mixed linear model for Sensory == 1
l3 <- lmer(Response ~ (SpeedR+SensoryR+TargetCellN) * (TargetSpeedR+TargetSensoryR) + SpeedR*Speed + (1|Session),data = subset(dataTable, Sensory == 1 & NonTargetCell == 1),REML = FALSE)
l4 <- lmer(Response ~ (SpeedR+SensoryR+TargetCellN) * (TargetSpeedR+TargetSensoryR) + SpeedR*Speed + AveDist + MinDist+(1|Session),data = subset(dataTable, Sensory == 1 & NonTargetCell == 1),REML = FALSE)


sink(file.path(base_dir, paste0(data_type,"_TrialBased_SensoryLME.txt")), append = FALSE)
print(summary(l3))
print(summary(l4))
print(anova(l4,l3))
sink()
}
