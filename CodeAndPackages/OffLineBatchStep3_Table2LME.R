#library(R.matlab)
library(lmerTest)

#library(lme4)

# Define base directory
base_dir <- "//nimhlabstore1.nimh.nih.gov/UFNC/FNC2/Zhang/Projects/Project-LocalProcessing/Step3/awakeRefSpon/GroupSLM23-Oct-2025"

# List subdirectories (assuming numeric names only)
subdirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)

for (subdir in subdirs) {
  # Define path to CSV
  csv_path <- file.path(subdir, "deltaFSLMGroupRes_Win1.csv")
  
  # Check if file exists before proceeding
  if (!file.exists(csv_path)) {
    message(paste("File not found:", csv_path))
    next
  }
  
  # Read data
  dataTable <- read.csv(csv_path)
  
  # Convert to factor
  dataTable$Group <- factor(dataTable$Group)
  dataTable$Cell <- factor(dataTable$Cell)
  dataTable$Sensory <- factor(dataTable$Sensory)
  dataTable$Session <- factor(dataTable$Session)

  # Mixed linear model for Sensory == 0
  l1 <- lmer(Response ~ (SpeedR+SensoryR) * (TargetSpeedR+TargetSensoryR) + SpeedR*Speed + MinDist + AveDist+(1|Session),data = subset(dataTable, Sensory == 0 & NonTargetCell == 1),REML = FALSE)
  l2 <- lmer(Response ~ (SpeedR+SensoryR) * (TargetSpeedR+TargetSensoryR) + SpeedR*Speed +(1|Session),data = subset(dataTable, Sensory == 0 & NonTargetCell == 1),REML = FALSE)


  sink(file.path(subdir, "NonSensoryLME.txt"), append = FALSE)
  print(summary(l1))
  print(summary(l2))
  print(anova(l1,l2))
  sink()

  # Mixed linear model for Sensory == 1
  l3 <- lmer(Response ~ (SpeedR+SensoryR) * (TargetSpeedR+TargetSensoryR) + SpeedR*Speed + MinDist + AveDist+(1|Session),data = subset(dataTable, Sensory == 1 & NonTargetCell == 1),REML = FALSE)
  l4 <- lmer(Response ~ (SpeedR+SensoryR) * (TargetSpeedR+TargetSensoryR) + SpeedR*Speed +(1|Session),data = subset(dataTable, Sensory == 1 & NonTargetCell == 1),REML = FALSE)


  sink(file.path(subdir, "SensoryLME.txt"), append = FALSE)
  print(summary(l3))
  print(summary(l4))
  print(anova(l4,l3))
  sink()
}
