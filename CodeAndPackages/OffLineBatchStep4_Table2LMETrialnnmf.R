#library(R.matlab)
library(lmerTest)
library(emmeans)
dataTable=read.csv("//nimhlabstore1.nimh.nih.gov/UFNC/FNC2/Zhang/Projects/Project-LocalProcessing/Step4/awakeRefSpon/NonNegMatFac/NonTargetSLMScoreTrialDynWin1.csv");

dataTable$Group=factor(dataTable$Group)
dataTable$Cell=factor(dataTable$Cell)
#dataTable$Sensory=factor(dataTable$Sensory)
dataTable$Session=factor(dataTable$Session)
#dataTable$PowerZero=factor(dataTable$PowerZero)

#dataTable$Session=factor(dataTable$Session)

l11 <- lmer(SpeedScore ~ Group*Speed + (1|Session),data = subset(dataTable, Sensory == 0 & PowerZero == 0),REML = FALSE)
summary(l11)

emtr <-emtrends(l11, ~ Group, var = "Speed")
summary(emtr, infer = TRUE)

dataT1<-dataTable
dataT2<-dataTable
l11 <- lmer(SpeedScore ~ Group * Speed + (1 | Session),
            data = subset(dataTable, Sensory == 0 & PowerZero == 0),
            REML = FALSE)
summary(l11)

dataT2$Group <- relevel(dataT2$Group, ref = "2")

l11 <- lmer(SpeedScore_SpeedReg ~ Group  + (1 | Session),
            data = subset(dataT1, Sensory == 0 & PowerZero == 0 & abs(Speed) <0.5),
            REML = TRUE)
summary(l11)



l11 <- lm(WhiskerScore ~ Group+Speed + (1|Session),data = subset(dataT1, Sensory == 1 & PowerZero == 0),REML = FALSE)
summary(l11)






l11 <- lm(WhiskerScore ~ Group+Speed + (1|Session),data = subset(dataTable, Sensory == 1 & PowerZero == 0),REML = FALSE)
summary(l11)


aov_model <- aov(WhiskerScore ~ Group,
                 data = subset(dataTable, Sensory == 1 & PowerZero == 0))
summary(aov_model)






l11 <- lmer(SpeedScore_SpeedReg ~ TargetSpeedR+Speed + (1|Session),data = subset(dataTable, Sensory == 0 & PowerZero == 0),REML = FALSE)
summary(l11)

l11 <- lmer(Response ~ (SpeedR+SensoryR) * (TargetSpeedR+TargetSensoryR+Speed) + (1|Session),data = subset(dataTable, Sensory == 0 & NonTargetCell == 1),REML = FALSE)
l1 <- lmer(Response ~ (SpeedR+SensoryR+Speed) * (TargetSpeedR+TargetSensoryR) + (1|Session),data = subset(dataTable, Sensory == 0 & NonTargetCell == 1),REML = FALSE)

l1 <- lmer(Response ~ (SpeedR+SensoryR) * (TargetSpeedR+TargetSensoryR) + Speed*SpeedR+(1|Session),data = subset(dataTable, Sensory == 0 & NonTargetCell == 1),REML = FALSE)
l2 <- lmer(Response ~ (SpeedR+SensoryR) * (TargetSpeedR+TargetSensoryR) + Speed*SpeedR+(1|Session),data = subset(dataTable, Sensory == 1 & NonTargetCell == 1),REML = FALSE)


l11 <- lmer(Response ~ (SpeedR) * (TargetSpeedR+TargetSensoryR) + Speed*SpeedR + (1|Session),data = subset(dataTable, Sensory == 0 & NonTargetCell == 1),REML = FALSE)
l11 <- lmer(Response ~ (SpeedR) * (TargetSpeedR+TargetSensoryR) + Speed*SpeedR + (1|Session),data = subset(dataTable, Sensory == 0 & NonTargetCell == 1),REML = FALSE)

#l1 <- lmer(Response ~ SpeedR * (TargetSpeedR+TargetSensoryR) + Speed + (1|Session),data = subset(dataTable, Sensory == 0 & NonTargetCell == 1),REML = FALSE)
#l2 <- lmer(Response ~ SpeedR * TargetSpeedR+SpeedR *TargetSensoryR + Speed + (1|Session),data = subset(dataTable, Sensory == 0 & NonTargetCell == 1),REML = FALSE)
#l3 <- lmer(Response ~ SpeedR * TargetSpeedR+SpeedR *TargetSensoryR + Speed + (1|Session),data = subset(dataTable, Sensory == 0 & NonTargetCell == 1),REML = FALSE)


#l2 <- lmer(Response ~ (SpeedR+SensoryR) * (TargetSpeedR+TargetSensoryR) + Speed + (1|(Session),data = subset(dataTable, Sensory == 1 & NonTargetCell == 1),REML = FALSE)

l2 <- lmer(Response ~ (SpeedR+SensoryR) * (TargetSpeedR+TargetSensoryR) + Speed*SpeedR + (1|Session),data = subset(dataTable, Sensory == 1 & NonTargetCell == 1),REML = FALSE)

l2 <- lmer(Response ~ (SpeedR+SensoryR+Speed) * (TargetSpeedR+TargetSensoryR) + (1|Session),data = subset(dataTable, Sensory == 1 & NonTargetCell == 1),REML = FALSE)

l2 <- lmer(Response ~ (SpeedR+SensoryR+Spee) * (TargetSpeedR+TargetSensoryR) + Speed*SpeedR+(1|Session),data = subset(dataTable, Sensory == 1 & NonTargetCell == 1),REML = FALSE)
l21 <- lmer(Response ~ (SpeedR+SensoryR) * (TargetSpeedR+TargetSensoryR) + Speed+(1|Session),data = subset(dataTable, Sensory == 1 & NonTargetCell == 1),REML = FALSE)

sink("//nimhlabstore1.nimh.nih.gov/UFNC/FNC2/Zhang/Projects/Project-LocalProcessing/Step3/awakeRefSpon/GroupSLM17-Oct-2025/NonSensoryLMESpeedTh1.txt",append=TRUE);
print(summary(l1));
sink()

sink("//nimhlabstore1.nimh.nih.gov/UFNC/FNC2/Zhang/Projects/Project-LocalProcessing/Step3/awakeRefSpon/GroupSLM17-Oct-2025/SensoryLMESpeedTh1.txt",append=TRUE);
print(summary(l2));
sink()


sink()
sink()
summary(l1)
summary(l2)

l3 <- lmer(Response ~ (SpeedR+SensoryR)*Group + Speed + (1|(Session)),data = subset(dataTable, Sensory == 0 & NonTargetCell == 1),REML = FALSE)
l4 <- lmer(Response ~ (SpeedR+SensoryR)*Group + Speed + (1|(Session)),data = subset(dataTable, Sensory == 1 & NonTargetCell == 1),REML = FALSE)


summary(l3)
summary(l4)


l3 <- lmer(Response ~ (SpeedR+SensoryR) * (TargetSpeedR+TargetSensoryR) + Speed + (1|Session),data = subset(dataTable, Sensory == 0 & NonTargetCell == 1 & Group != 3),REML = FALSE)

summary(l1)
summary(l2)

NoW1HF=lmer(Speed ~ TargetSpeedR*Time+(1|Session), data=subset(dataTable,Whisk==0&TargetCell==0),REML=FALSE)
summary(NoW1HF)


NoW1HF=lmer(Response ~ SpeedR*TargetSpeedR*Speed*Time+(1|Session), data=subset(dataTable,Whisk==0&TargetCell==0),REML=FALSE)
summary(NoW1HF)



NoW1HF=lmer(Response ~ SpeedR*TargetSpeedR*Speed*Time+(1|Session), data=subset(dataTable,Whisk==0&TargetCell==0&Time>0.5),REML=FALSE)
summary(NoW1HF)

NoW1HF=lmer(Response ~ SpeedR*TargetSpeedR*Speed+(1|Session), data=subset(dataTable,Whisk==0&TargetCell==0&Time<0.2&abs(Speed)<0.1),REML=FALSE)
summary(NoW1HF)

NoW1HF=lmer(Response ~ SpeedR*TargetSpeedR*Speed+(1|Session), data=subset(dataTable,Whisk==0&TargetCell==0&Time<0.4&Time>0.0&abs(Speed)<0.2),REML=FALSE)
summary(NoW1HF)

NoW2HF=lmer(Response ~ SpeedR*TargetSpeedR+SpeedR*Speed+(1|Session), data=subset(dataTable,Whisk==0&TargetCell==0&Time<0.4&Time>0.0&abs(Speed)<0.2),REML=FALSE)
summary(NoW2HF)

anova(NoW1HF,NoW2HF)


NoW1HF=lmer(Response ~ SpeedR*TargetSpeedR*Speed+(1|Session), data=subset(dataTable,Whisk==0&TargetCell==0&abs(Speed)<0.1&Time<0.4),REML=FALSE)
summary(NoW1HF)

NoW2HF=lmer(Response ~ SpeedR*TargetSpeedR*Time+Speed+(1|Session), data=subset(dataTable,Whisk==0&TargetCell==0&Time>0.50),REML=FALSE)
summary(NoW2HF)

NoW1=lmer(Response ~ SpeedR*TargetSpeedR*Time*Speed+(1|Session), data=subset(dataTable,Whisk==0&TargetCell==0&Speed<0.1),REML=TRUE)
summary(NoW1)

NoW1=lmer(Response ~ SpeedR*TargetSpeedR*Time+(1|Session), data=subset(dataTable,Whisk==0&TargetCell==0),REML=TRUE)
summary(NoW1)


NoW1=lmer(Response ~ SpeedR*TargetSpeedR*Time+SpeedR*TargetSpeedR*Speed+SpeedR*Speed*Time+(1|Session), data=subset(dataTable,Whisk==0&TargetCell==0&Speed<0.1),REML=FALSE)
summary(NoW2)

NoW2=lmer(Response ~ SpeedR*TargetSpeedR*Time+SpeedR*TargetSpeedR*Speed+SpeedR*Speed*Time+(1|Session), data=subset(dataTable,Whisk==0&TargetCell==0&Speed<0.1),REML=TRUE)
summary(NoW2)

anova(NoW2,NoW1)

NoW3=lmer(Response ~ SpeedR*TargetSpeedR*Time+SpeedR*Speed+(1|Session), data=subset(dataTable,Whisk==0&TargetCell==0),REML=FALSE)
summary(NoW3)


NoW4=lmer(Response ~ SpeedR*TargetSpeedR*Time+SpeedR*Speed*Time+(1|Session), data=subset(dataTable,Whisk==0&TargetCell==0&abs(Speed)<0.2),REML=FALSE)
summary(NoW4)

anova(NoW1,NoW4)



NoWFull=lmer(Response ~ StimR*TargetStimR*Time+SpeedR*TargetSpeedR*Time+Speed+(1|Session), data=subset(dataTable,Whisk==0&TargetCell==0),REML=FALSE)
summary(NoWFull)
anova(NoW,NoWFull)

l11=lmer(Response ~ SpeedR*TargetSpeedR*Time+StimR*TargetStimR*Time+Speed+(1|Session), data=subset(dataTable,Whisk==0&TargetCell==0),REML=FALSE)

summary(l1)

anova(l1,l11)

l1=lmer(Response ~ StimR*TargetStimR*Time+Speed+(1|Session), data=subset(dataTable,Whisk==1&TargetCell==0),REML=FALSE)
summary(l1)
l11=lmer(Response ~ SpeedR*TargetSpeedR*Time+StimR*TargetStimR*Time+Speed+(1|Session), data=subset(dataTable,Whisk==1&TargetCell==0),REML=FALSE)



W1HF=lmer(Response ~ StimR*TargetStimR*Time+Speed+(1|Session), data=subset(dataTable,Whisk==1&TargetCell==0&Time<0.50),REML=FALSE)
summary(W1HF)
W2HF=lmer(Response ~ StimR*TargetStimR*Time+Speed+(1|Session), data=subset(dataTable,Whisk==1&TargetCell==0&Time>0.50),REML=FALSE)
summary(W2HF)

W=lmer(Response ~ StimR*TargetStimR*Time+Speed+(1|Session), data=subset(dataTable,Whisk==1&TargetCell==0),REML=FALSE)
summary(W)

WFull=lmer(Response ~ StimR*TargetStimR*Time+SpeedR*TargetSpeedR*Time+Speed+(1|Session), data=subset(dataTable,Whisk==1&TargetCell==0),REML=FALSE)
summary(WFull)

W1=lmer(Response ~ StimR*TargetStimR*Time+SpeedR*Time+Speed+(1|Session), data=subset(dataTable,Whisk==1&TargetCell==0&abs(Speed)<2),REML=FALSE)

anova(W,WFull)




sink("//nimhlabstore1.nimh.nih.gov/UFNC/FNC2/Zhang/Projects/Project-LocalProcessing/Step3/awakeRefSpon/GroupSLM6Sessions/NonWhiskLME.txt",append=FALSE);
#print(summary(l1));
l11=lmer(Response ~ SpeedR*TargetSpeedR+Speed+(1|Session), data=subset(dataTable,Whisk==1),REML=FALSE)
#sink("//nimhlabstore1.nimh.nih.gov/UFNC/FNC2/Zhang/Projects/Project-LocalProcessing/Step3/awakeRefSpon/GroupSLM6Sessions/NonWhiskLME.txt",append=TRUE);
#print(summary(l1));

sink()

anova(l1,l11)


l2=lmer(Response ~ SpeedR*TargetSpeedR+StimR*TargetStimR+Speed+(1|Session), data=subset(dataTable,Whisk==2&TargetCell==0),REML=FALSE)
l2=lmer(Response ~ StimR*TargetStimR*Time+SpeedR*TargetSpeedR*Speed*Time+(1|Session), data=subset(dataTable,Whisk==1&TargetCell==0),REML=FALSE)


sink("//nimhlabstore1.nimh.nih.gov/UFNC/FNC2/Zhang/Projects/Project-LocalProcessing/Step3/awakeRefSpon/GroupSLM6Sessions/WhiskLME.txt",append=FALSE);
print(summary(l2));
sink()


#l2=lmer(Response ~ SpeedR*StimR*StimGroup+Dist+(1|Cell), data=dataTable,REML=FALSE)
#l3=lmer(Firing ~ Speed*Spa+(1|CellSubj), data=dataTable,REML=FALSE)
#A=(anova(l3,l2));
#sink("C:/Users/lzhang481/SingerLab/Projects/Project7-Aging/SF180WRfiringRate/statisText/tempLME.txt",append=TRUE);
#print(A);
#sink();