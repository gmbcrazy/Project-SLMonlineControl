
clear all


TBLDataPath='\\nimhlabstore1.nimh.nih.gov\UFNC\FNC2\Zhang\Projects\Project-LocalProcessing\Step3\awakeRefSpon\GroupSLM23-Oct-2025\';

SavePath=['\\nimhlabstore1.nimh.nih.gov\UFNC\FNC2\Zhang\Projects\Project-LocalProcessing\Step4\' Suite2pDataKeywords '\'];
mkdir(SavePath);


load([TBLDataPath 'FOVoutputs.mat'])
load([TBLDataPath 'ValidSessionGroup.mat'])

animalpaths={FOVUpdate.DataFolder};

for i = 1:length(animalpaths)
    % Extract animal ID (e.g., SL0886)
    animalMatch = regexp(animalpaths{i}, '(SL\d{4}|L\d{5})', 'match');
    if ~isempty(animalMatch)
        animalIDs{i} = animalMatch{1};
    else
        animalIDs{i} = '';
    end

    % Extract date (8-digit number)
    dateMatch = regexp(animalpaths{i}, '\d{8}', 'match');
    if ~isempty(dateMatch)
        sessionDates{i} = dateMatch{1};
    else
        sessionDates{i} = '';
    end

    % suite2pFOVPath{i}=[Suite2pSaveFolderAll '\' animalIDs{i} sessionDates{i} '\'];
    % mkdir(suite2pFOVPath{i})
    % suite2pFOVPathLocal{i}=[Suite2pSaveFolderAllLocal '\' animalIDs{i} sessionDates{i} '\'];
    % mkdir(suite2pFOVPathLocal{i})

end

ProcessPar.PlotFigure=0;
ProcessPar.GroupLabel={'L','S','N'};
ProcessPar.GroupList=[1 2 3];
ProcessPar.GroupColor=[255 51 153;91 20 212;121 247 111]/255;
ProcessPar.PowerZeroColor=[0.5 0.5 0.5];
ProcessPar.PowerZero=[0 1];
ProcessPar.PowerZeroLabel={'SLM','FakeSLM'};
ProcessPar.VolOut=[0 1];
ProcessPar.VolOutLabel={'NoWhisk','Whisk'}
ProcessPar.AwakeState=[1];
ProcessPar.AwakeStateLabel={'Awake'};

SLMGroupTableTrial=readtable([TBLDataPath 'deltaFSLMGroupResTrialDynWin1.csv']);



DimentionReductionMethod='NonNegMatFac';
SaveFunFOVDim=[SavePath DimentionReductionMethod '\'];
mkdir(SaveFunFOVDim)


SponInd=[1:8200];
ExcludeStimInd=[1:2000 2800:4100];  
ExcludeStimInd=union(ExcludeStimInd,ExcludeStimInd+4100);

% AddStimNoise=random('norm',0,0.1,length(SponInd),1);

SaveFunSub=[SaveFunFOVDim 'FOV\'];
mkdir(SaveFunSub)
% clear FOVres

SaveFunFOVDim
iData=1;
lowDimN=40;
ScoreLabel={'SpeedScore','WhiskScore'}
NDataName={'deltaF','spks'};


% for iFOV=3:3

IndexFOVNeed=unique(ValidSG.Session)
clear FOVres
for iFOV=1:length(IndexFOVNeed)
    FOVres(iFOV) = FOVnnmf_SpeedWStim(IndexFOVNeed(iFOV), SLMGroupTableTrial, NeuroTrace, BehTrace, iData, lowDimN, SponInd, ExcludeStimInd);
    FOVsavePath=[SaveFunSub 'FOV' num2str(IndexFOVNeed(iFOV)) animalIDs{IndexFOVNeed(iFOV)} sessionDates{IndexFOVNeed(iFOV)}];
    plotFOVnnmf(FOVres(iFOV), SponInd, BehTrace(IndexFOVNeed(iFOV)), ProcessPar, FOVsavePath);
    close all
end

TBLall=[];
TBLallNontarget=[];
TBLall_subNontarget=[];
clear FOVresTest
for iFOV=1:length(IndexFOVNeed)
    [FOVresTest(iFOV), ScoreGroup, SpeedGroup, maxScore, TBL_rows] = SLMGroup_NNMFScoresFOV(FOVres(iFOV).FOVid, SLMGroupTableTrial,FOVres(iFOV), ProcessPar);
    TBLall=[TBLall;TBL_rows.all];
    TBLallNontarget=[TBLallNontarget;TBL_rows.nontarget];
    TBLall_subNontarget=[TBLall_subNontarget;TBL_rows.sub_NonTarget];
end
AllTbl=innerjoin(TBLall , unique(ValidSG(:,{'Session','Group'})), ...
                     'Keys', {'Session','Group'});
NTTbl=innerjoin(TBLallNontarget , unique(ValidSG(:,{'Session','Group'})), ...
                     'Keys', {'Session','Group'});
NTTbl_sub=innerjoin(TBLall_subNontarget , unique(ValidSG(:,{'Session','Group'})), ...
                     'Keys', {'Session','Group'});
DataForCompare={AllTbl NTTbl_sub};



save([SaveFunFOVDim date NDataName{iData} 'nnmfResult.mat'],'-v7.3');


iWin=1;
writetable(TBLall_subNontarget,[SaveFunFOVDim 'NonTargetSLMScoreTrialDynWin' num2str(iWin) '.csv']);
writetable(TBLall,[SaveFunFOVDim 'AllSLMScoreTrialDynWin' num2str(iWin) '.csv']);


clear tempD statsSpeed

close all
figure;
   GroupPair.CorrName='fdr';
   GroupPair.Q=0.1;
   GroupPair.Pair=[];
   GroupPair.SignY=0.006;
   GroupPair.Plot=1;
   GroupPair.Std=0;      %%%%%%%%%using standard deviation as errorbar
   GroupPair.SamplePlot=1; %%%%%%%%%Plot Individual Sample Point
   GroupPair.SamplePairedPlot=1; %%%%%%%%%Dash line for paired comparison sample
   GroupPair.Marker={'o'};   % Datatype=varargin{1};  %%%%%%Datatype=0, non-paired data; 1 for paired data

GroupPair.SignY=0.03;
GroupPair.LimY=[-0.0001 GroupPair.SignY*1.2];
GroupPair.ViolinLR=[0 0 0 0];
GroupPair.ANOVA=0;
GroupPair.Test='RankTest';
% GroupPair.Crit_p
clear statsSpeed;
close all
SpeedTh=0.5;
showYLim=[-1;1]*GroupPair.LimY(end)
for iWhisk=1:length(ProcessPar.VolOut)
    for jtype=1:length(DataForCompare)
        clear tempDGroup
        
        tempD0=DataForCompare{jtype};
        tempD=tempD0;

        % tempD=groupsummary(tempD0,{'Session','Group','Sensory','PowerZero'},'mean');
        % tempD= groupsummaryBack2OldNames(tempD0, tempD, 'mean');
        % 

        % tempDNonZ=tempD(tempD.Sensory==ProcessPar.VolOut(iWhisk)&tempD.PowerZero==0,:);
        % tempDZero=tempD(tempD.Sensory==ProcessPar.VolOut(iWhisk)&tempD.PowerZero==1,:);
        % 

        % tempD=tempD(tempD.Sensory==ProcessPar.VolOut(iWhisk)&tempD.PowerZero==0&abs(tempD.Speed)<SpeedTh,:);
        % tempDNonZ=tempD(tempD.Sensory==ProcessPar.VolOut(iWhisk)&tempD.PowerZero==0&tempD.Speed<SpeedTh,:);
        % tempDZero=tempD(tempD.Sensory==ProcessPar.VolOut(iWhisk)&tempD.PowerZero==1&tempD.Speed<SpeedTh,:);

        tempDNonZ=tempD(tempD.Sensory==ProcessPar.VolOut(iWhisk)&tempD.PowerZero==0&abs(tempD.Speed)<SpeedTh,:);
        tempDZero=tempD(tempD.Sensory==ProcessPar.VolOut(iWhisk)&tempD.PowerZero==1&abs(tempD.Speed)<SpeedTh,:);
        % 
        for iGroup=1:3
            tempDGroup{iGroup}=tempDNonZ.SpeedScore_SpeedReg(tempDNonZ.Group==iGroup);
        end
            tempDGroup{iGroup+1}=tempDZero.SpeedScore_SpeedReg;
        % % 
        for iGroup=1:3
            tempDGroup{iGroup}=tempDNonZ.SpeedScore(tempDNonZ.Group==iGroup);
        end
            tempDGroup{iGroup+1}=tempDZero.SpeedScore;

      
        subplotLU(2,length(DataForCompare),iWhisk,jtype);
        % statsSpeed(iWhisk,jtype)=ErrorBarPlotLU(1:4,tempDGroup,[],[ProcessPar.GroupColor;0.7 0.7 0.7],2,0,[],GroupPair,1:4);
        statsSpeed(iWhisk,jtype)=ErrorViolinHalf(1:4,tempDGroup,[ProcessPar.GroupColor;0.7 0.7 0.7],0,[],GroupPair,1:4);

        % statsSpeed(iWhisk,jtype)=ErrorBoxPlotLU(1:4,tempDGroup,[ProcessPar.GroupColor;0.7 0.7 0.7],[],GroupPair,1:4)
        ylim(showYLim);
        yticks(union(0,showYLim))
        if jtype==1
        yticklabels(union(0,showYLim))
        end
         if iWhisk==2
             xticklabels({'L','S','N','0'})
        end
    end

end
statsSpeed(1,2)


close all
clear statsSensory
GroupPair.SignY=0.05;
GroupPair.LimY=[-0.0001 GroupPair.SignY*1.2];

showYLim=[-1;1]*GroupPair.LimY(end)

SpeedTh=1;
showYLim=[-1;1]*GroupPair.LimY(end)
for iWhisk=1:length(ProcessPar.VolOut)
    for jtype=1:length(DataForCompare)
        clear tempDGroup
        
        tempD0=DataForCompare{jtype};
        tempD=tempD0;

        % tempD=groupsummary(tempD0,{'Session','Group','Sensory','PowerZero'},'mean');
        % tempD= groupsummaryBack2OldNames(tempD0, tempD, 'mean');
        % 

        % tempDNonZ=tempD(tempD.Sensory==ProcessPar.VolOut(iWhisk)&tempD.PowerZero==0,:);
        % tempDZero=tempD(tempD.Sensory==ProcessPar.VolOut(iWhisk)&tempD.PowerZero==1,:);
        % 

        % tempD=tempD(tempD.Sensory==ProcessPar.VolOut(iWhisk)&tempD.PowerZero==0&abs(tempD.Speed)<SpeedTh,:);
        tempDNonZ=tempD(tempD.Sensory==ProcessPar.VolOut(iWhisk)&tempD.PowerZero==0&tempD.Speed<SpeedTh,:);
        tempDZero=tempD(tempD.Sensory==ProcessPar.VolOut(iWhisk)&tempD.PowerZero==1&tempD.Speed<SpeedTh,:);

        tempDNonZ=tempD(tempD.Sensory==ProcessPar.VolOut(iWhisk)&tempD.PowerZero==0&abs(tempD.Speed)<SpeedTh,:);
        tempDZero=tempD(tempD.Sensory==ProcessPar.VolOut(iWhisk)&tempD.PowerZero==1&abs(tempD.Speed)<SpeedTh,:);

        for iGroup=1:3
            tempDGroup{iGroup}=tempDNonZ.WhiskerScore_SpeedReg(tempDNonZ.Group==iGroup);
        end
            tempDGroup{iGroup+1}=tempDZero.WhiskerScore_SpeedReg;
        % % 
        % for iGroup=1:3
        %     tempDGroup{iGroup}=tempDNonZ.SpeedScore(tempDNonZ.Group==iGroup);
        % end
        %     tempDGroup{iGroup+1}=tempDZero.SpeedScore;
        % 
      
        subplotLU(2,length(DataForCompare),iWhisk,jtype);
        % statsSpeed(iWhisk,jtype)=ErrorBarPlotLU(1:4,tempDGroup,[],[ProcessPar.GroupColor;0.7 0.7 0.7],2,0,[],GroupPair,1:4);
        statsSensory(iWhisk,jtype)=ErrorViolinHalf(1:4,tempDGroup,[ProcessPar.GroupColor;0.7 0.7 0.7],0,[],GroupPair,1:4);

        % statsSpeed(iWhisk,jtype)=ErrorBoxPlotLU(1:4,tempDGroup,[ProcessPar.GroupColor;0.7 0.7 0.7],[],GroupPair,1:4)
        ylim(showYLim);
        yticks(union(0,showYLim))
        if jtype==1
        yticklabels(union(0,showYLim))
        end
         if iWhisk==2
             xticklabels({'L','S','N','0'})
        end
    end

end
statsSensory(1,2)


figure;

SpeedTh=1;
clear statsSensory
for iWhisk=1:length(ProcessPar.VolOut)
    for jtype=1:length(DataForCompare)
        clear tempDGroup
        tempD=DataForCompare{jtype};
        % tempD=tempD(tempD.Sensory==ProcessPar.VolOut(iWhisk)&tempD.PowerZero==0&abs(tempD.Speed)<SpeedTh,:);
        % tempDNonZ=tempD(tempD.Sensory==ProcessPar.VolOut(iWhisk)&tempD.PowerZero==0&tempD.Speed<SpeedTh,:);
        % tempDZero=tempD(tempD.Sensory==ProcessPar.VolOut(iWhisk)&tempD.PowerZero==1&tempD.Speed<SpeedTh,:);

        tempDNonZ=tempD(tempD.Sensory==ProcessPar.VolOut(iWhisk)&tempD.PowerZero==0&abs(tempD.Speed)<SpeedTh,:);
        tempDZero=tempD(tempD.Sensory==ProcessPar.VolOut(iWhisk)&tempD.PowerZero==1&abs(tempD.Speed)<SpeedTh,:);

        for iGroup=1:3
            tempDGroup{iGroup}=tempDNonZ.WhiskerScore_SpeedReg(tempDNonZ.Group==iGroup);
        end
            tempDGroup{iGroup+1}=tempDZero.WhiskerScore_SpeedReg;


      
        subplotLU(2,length(DataForCompare),iWhisk,jtype);
        % statsSpeed(iWhisk,jtype)=ErrorBarPlotLU(1:3,tempDGroup,[],ProcessPar.GroupColor,2,0,[],GroupPair);
        % statsSpeed(iWhisk,jtype)=ErrorViolinHalf(1:3,tempDGroup,ProcessPar.GroupColor,0,[],GroupPair,1:3);

        statsSensory(iWhisk,jtype)=ErrorBoxPlotLU(1:4,tempDGroup,[ProcessPar.GroupColor;0.7 0.7 0.7],[],GroupPair,1:4)
        ylim(showYLim);
        yticks(union(0,showYLim))
        if jtype==1
        yticklabels(union(0,showYLim))
        end
         if iWhisk==2
             xticklabels({'L','S','N','0'})
        end
    end

end




stats=ErrorViolinHalf(x,y_data,PlotColor,varargin)



figure;
GroupPair.SignY=0.00005;
GroupPair.LimY=[-0.00005 0.00005];
SpeedTh=1;
for iWhisk=1:length(ProcessPar.VolOut)

    for jtype=1:length(DataForCompare)
        tempD=DataForCompare{jtype};
        tempD=tempD(tempD.Sensory==ProcessPar.VolOut(iWhisk)&tempD.PowerZero==0&abs(tempD.Speed)<SpeedTh,:);
        for iGroup=1:3
            tempDGroup{iGroup}=tempD.WhiskerScore(tempD.Group==iGroup);
        end
          
        subplotLU(2,length(DataForCompare),iWhisk,jtype);
        statsSensory(iWhisk,jtype)=ErrorBarPlotLU(1:3,tempDGroup,[],ProcessPar.GroupColor,2,0,[],GroupPair);
    end


end




NTTbl=NTTbl(NTTbl.Sensory==0&NTTbl.Group==2,:);
[r,p]=corr(NTTbl.TargetSpeedR,NTTbl.SpeedScore)
figure;
plot(NTTbl.TargetSpeedR,NTTbl.SpeedScore,'.')


AllTbl=AllTbl(AllTbl.Sensory==0&AllTbl.Group==1,:);
[r,p]=partialcorr(AllTbl.TargetSpeedR,AllTbl.SpeedScore,[AllTbl.Speed AllTbl.SpeedR])
figure;
plot(AllTbl.TargetSpeedR,AllTbl.SpeedScore,'.')




   GroupPair.CorrName='fdr';
   GroupPair.Q=0.1;
   GroupPair.Pair=[];
   GroupPair.SignY=0.006;
   GroupPair.Plot=1;
   GroupPair.Std=0;      %%%%%%%%%using standard deviation as errorbar
   GroupPair.SamplePlot=1; %%%%%%%%%Plot Individual Sample Point
   GroupPair.SamplePairedPlot=1; %%%%%%%%%Dash line for paired comparison sample
   GroupPair.LimY=[-0.006 GroupPair.SignY*1.2];
   GroupPair.Marker={'o'};   % Datatype=varargin{1};  %%%%%%Datatype=0, non-paired data; 1 for paired data


GroupMethod='mean';
SpeedTh=1
TBLallavg = groupsummary(TBLall(AllTbl.Sensory==0&AllTbl.PowerZero==0&abs(AllTbl.Speed)<SpeedTh,:), {'Session', 'Group'}, GroupMethod);
TBL_subavg= groupsummary(NTTbl(NTTbl.Sensory==0&NTTbl.PowerZero==0&abs(NTTbl.Speed)<SpeedTh,:),{'Session', 'Group'}, GroupMethod);

TBLallavg = groupsummary(TBLall(AllTbl.Sensory==0&AllTbl.PowerZero==0&AllTbl.Speed<SpeedTh,:), {'Session', 'Group'}, GroupMethod);
TBLallavg=groupsummaryBack2OldNames(AllTbl, TBLallavg,GroupMethod);

TBL_subavg= groupsummary(NTTbl(NTTbl.Sensory==0&NTTbl.PowerZero==0&NTTbl.Speed<SpeedTh,:),{'Session', 'Group'}, GroupMethod);
TBL_subavg=groupsummaryBack2OldNames(NTTbl, TBL_subavg,GroupMethod);


% Define grouping
[G, sessionVals, TargetcellVals] = findgroups(newtbl2.Session, newtbl2.PointTargetCell);
% [G, sessionVals, TargetcellVals] = findgroups(newtbl2.Session, newtbl2.PointTargetCell);

% Define custom function
myfun = @(x,y) sum(x.*y) / sum(x);
myfun2 = @(x,y,z) sum((x.*z)-sum(y.*z)) / sum(x+y);

% Apply for each pair
AIabs_SpeedR  = splitapply(myfun, newtbl2.ActivateI,    newtbl2.SpeedR, G);
AIP_SpeedR = splitapply(myfun, newtbl2.ActivateIPos, newtbl2.SpeedR, G);
AIN_SpeedR = splitapply(myfun, newtbl2.ActivateINeg, newtbl2.SpeedR, G);
AI_SpeedR =  splitapply(myfun2, newtbl2.ActivateIPos,newtbl2.ActivateINeg, newtbl2.SpeedR, G);


ProcessPar.GroupColor

   Param.Color=ProcessPar.GroupColor;
   Param.Marker='o';
   Param.MarkerSize=10;
   Param.Rtype='pearson';
   Param.xLim=[min(TBLallavg.TargetSpeedR) max(TBLallavg.TargetSpeedR)];
   Param.yLim=[min(TBLallavg.SpeedScore) max(TBLallavg.SpeedScore)];
   Param.xLabel=[];
   Param.yLabel=[];

figure;
[OutPut,r,p,h]=LuPairRegressPlot_Group_Cov(TBLallavg.TargetSpeedR,TBLallavg.SpeedScore,TBLallavg.Speed,TBLallavg.Group,Param)
  
figure;
[OutPut,r,p,h]=LuPairRegressPlot_Group_Cov(TBL_subavg.TargetSpeedR,TBL_subavg.SpeedScore,TBL_subavg.Speed,TBL_subavg.Group,Param)

SpeedTh=1;

close all
Param2=Param;
Param2.Color=Param.Color(1:2,:);
figure;
[OutPut,r,p,h]=LuPairRegressPlot_Group_Cov(TBLallavg.TargetSpeedScore,TBLallavg.SpeedScore,TBLallavg.Speed,TBLallavg.Group,Param2)
  
figure;
[OutPut,r,p,h]=LuPairRegressPlot_Group_Cov(TBLall.TargetSpeedR,TBLall.SpeedScore,TBLall.Speed,TBLall.Group,Param)
figure;
[OutPut,r,p,h]=LuPairRegressPlot_Group_Cov(TBLall.TargetSpeedR,TBLall.SpeedScore,TBLall.Speed,TBLall.Group,Param)



SpeedTh=2;
TBLallTemp=TBLall(TBLall.Whisk==0&TBLall.PowerZero==0&abs(TBLall.Speed)<SpeedTh,:);
TBLsubTemp=TBLall_sub(TBLall_sub.Whisk==0&TBLall_sub.PowerZero==0&abs(TBLall_sub.Speed)<SpeedTh,:);
figure;
[OutPut,r,p,h]=LuPairRegressPlot_Group_Cov(TBLallTemp.TargetSpeedR,TBLallTemp.SpeedScore,TBLallTemp.Speed,TBLallTemp.Group,Param)
figure;
[OutPut,r,p,h]=LuPairRegressPlot_Group_Cov(TBLsubTemp.TargetSpeedR,TBLsubTemp.SpeedScore,TBLsubTemp.Speed,TBLsubTemp.Group,Param)


SpeedTh=2;
TBLallTemp=TBLall(TBLall.Whisk==0&TBLall.PowerZero==0&TBLall.Speed<SpeedTh,:);
TBLsubTemp=TBLall_sub(TBLall_sub.Whisk==0&TBLall_sub.PowerZero==0&TBLall_sub.Speed<SpeedTh,:);
figure;
[OutPut,r,p,h]=LuPairRegressPlot_Group_Cov(TBLallTemp.TargetSpeedR,TBLallTemp.SpeedScore,TBLallTemp.Speed,TBLallTemp.Group,Param)
figure;
[OutPut,r,p,h]=LuPairRegressPlot_Group_Cov(TBLsubTemp.TargetSpeedR,TBLsubTemp.SpeedScore,TBLsubTemp.Speed,TBLsubTemp.Group,Param)


TBLallTemp=TBLall(TBLall.Whisk==1&TBLall.PowerZero==0&TBLall.Speed<SpeedTh,:);
TBLsubTemp=TBLall_sub(TBLall_sub.Whisk==1&TBLall_sub.PowerZero==0&TBLall_sub.Speed<SpeedTh,:);

figure;
[OutPut,r,p,h]=LuPairRegressPlot_Group_Cov(TBLallTemp.TargetStimR,TBLallTemp.WhiskerScore,TBLallTemp.Speed,TBLallTemp.Group,Param)
figure;
[OutPut,r,p,h]=LuPairRegressPlot_Group_Cov(TBLsubTemp.TargetStimR,TBLsubTemp.WhiskerScore,TBLsubTemp.Speed,TBLsubTemp.Group,Param)

figure;
[OutPut,r,p,h]=LuPairRegressPlot_Group_Cov(TBL_subavg.TargetSpeedScore,TBL_subavg.SpeedScore,TBL_subavg.Speed,TBL_subavg.Group,Param)


GroupMethod='median';

TBLallavg = groupsummary(TBLall(TBLall.Whisk==1&TBLall.Group<=3&TBLall.PowerZero==0&TBLall.Speed<SpeedTh,:), {'Session', 'Group'}, GroupMethod);
TBL_subavg= groupsummary(TBLall_sub(TBLall.Whisk==1&TBLall.Group<=3&TBLall.PowerZero==0&TBLall.Speed<SpeedTh,:),{'Session', 'Group'}, GroupMethod);
oldNames = TBLallavg.Properties.VariableNames;
% Find which names start with 'mean_'
meanVars = startsWith(oldNames, [GroupMethod '_']);
% Remove 'mean_' prefix
newNames = oldNames;
newNames(meanVars) = erase(oldNames(meanVars), [GroupMethod '_']);
TBLallavg.Properties.VariableNames=newNames;
TBL_subavg.Properties.VariableNames=newNames;



Param.Color=ProcessPar.GroupColor;
Param.Marker='o';
Param.MarkerSize=10;
Param.Rtype='pearson';
Param.xLim=[min(TBLallavg.TargetStimR) max(TBLallavg.TargetStimR)];
Param.yLim=[min(TBLallavg.WhiskerScore) max(TBLallavg.WhiskerScore)];
figure;
[OutPut,r,p,h]=LuPairRegressPlot_Group_Cov(TBLallavg.TargetStimR,TBLallavg.WhiskerScore,TBLallavg.Speed,TBLallavg.Group,Param)
figure;
[OutPut,r,p,h]=LuPairRegressPlot_Group_Cov(TBL_subavg.TargetStimR,TBL_subavg.WhiskerScore,TBL_subavg.Speed,TBL_subavg.Group,Param)




   GroupPair.CorrName='fdr';
   GroupPair.Q=0.1;
   GroupPair.Pair=[];
   GroupPair.SignY=0.01;
   GroupPair.Plot=1;
   GroupPair.Std=0;      %%%%%%%%%using standard deviation as errorbar
   GroupPair.SamplePlot=0; %%%%%%%%%Plot Individual Sample Point
   GroupPair.SamplePairedPlot=0; %%%%%%%%%Dash line for paired comparison sample
   GroupPair.LimY=[-0.005 GroupPair.SignY*1.6];
   GroupPair.Marker={'o'};   % Datatype=varargin{1};  %%%%%%Datatype=0, non-paired data; 1 for paired data
   % SavePath=varargin{2};  %%%%%%Save the stats to txt file
   % GroupPair=varargin{3};
   % RGroupID=varargin{4};
    GroupPair.ANOVA=0;

   GroupPairSpeed=GroupPair;
   GroupPairSpeed.SignY=0.01;
   GroupPairSpeed.LimY=[-0.2 GroupPairSpeed.SignY*1.2];



SpeedTh=1;
% A=TBLallNontarget(TBLallNontarget.Whisk==0&TBLallNontarget.PowerZero==0&TBLallNontarget.Speed<SpeedTh,:);
% B=TBLall(TBLall.Whisk==0&TBLall.PowerZero==0&TBLall.Speed<SpeedTh,:);
figure;
clear dataAll dataNonT dataAll_sub statsSpeed
for iWhisk=1:length(ProcessPar.VolOut)

tempNonT=TBLallNontarget(TBLallNontarget.Sensory==ProcessPar.VolOut(iWhisk)&TBLallNontarget.PowerZero==0&TBLallNontarget.Speed<SpeedTh,:);
tempAll=TBLall(TBLall.Sensory==ProcessPar.VolOut(iWhisk)&TBLall.PowerZero==0&TBLall.Speed<SpeedTh,:);
tempAll_sub=TBLall_subNontarget(TBLall_subNontarget.Sensory==ProcessPar.VolOut(iWhisk)&TBLall_subNontarget.PowerZero==0&TBLall_subNontarget.Speed<SpeedTh,:);
% 
% tempNonT=TBLallNontarget(TBLallNontarget.Whisk==ProcessPar.VolOut(iWhisk)&TBLallNontarget.PowerZero==0&abs(TBLallNontarget.Speed)<SpeedTh,:);
% tempAll=TBLall(TBLall.Whisk==ProcessPar.VolOut(iWhisk)&TBLall.PowerZero==0&abs(TBLall.Speed)<SpeedTh,:);
% tempAll_sub=TBLall_subNontarget(TBLall_subNontarget.Sensory==ProcessPar.VolOut(iWhisk)&TBLall_subNontarget.PowerZero==0&abs(TBLall_subNontarget.Speed<SpeedTh),:);


for iGroup=1:3
    dataAll{iGroup,iWhisk}=tempAll.SpeedScore(tempAll.Group==iGroup);
    dataNonT{iGroup,iWhisk}=tempNonT.SpeedScore(tempNonT.Group==iGroup);
    dataAll_sub{iGroup,iWhisk}=tempAll_sub.SpeedScore(tempAll_sub.Group==iGroup);

end

subplotLU(2,3,iWhisk,1)
% statsSpeed(iWhisk,1)=ErrorBoxPlotLU(1:3,dataAll(:,iWhisk),ProcessPar.GroupColor,[],GroupPair)
statsSpeed(iWhisk,1)=ErrorBarPlotLU(1:3,dataAll(:,iWhisk),[],ProcessPar.GroupColor,2,0,[],GroupPair)

figure;
a=ErrorBoxPlotLU(1:3,dataAll(:,iWhisk),ProcessPar.GroupColor,[],GroupPair)




set(gca,'ylim',GroupPair.LimY,'ytick',GroupPair.LimY);
subplotLU(2,3,iWhisk,2)
% statsSpeed(iWhisk,2)=ErrorBoxPlotLU(1:3,dataNonT(:,iWhisk),ProcessPar.GroupColor,[],GroupPair)
statsSpeed(iWhisk,2)=ErrorBarPlotLU(1:3,dataNonT(:,iWhisk),[],ProcessPar.GroupColor,2,0,[],GroupPair)

set(gca,'ylim',GroupPair.LimY,'ytick',GroupPair.LimY);
subplotLU(2,3,iWhisk,3)
% statsSpeed(iWhisk,3)=ErrorBoxPlotLU(1:3,dataAll_sub(:,iWhisk),ProcessPar.GroupColor,[],GroupPair)
statsSpeed(iWhisk,3)=ErrorBarPlotLU(1:3,dataAll_sub(:,iWhisk),[],ProcessPar.GroupColor,2,0,[],GroupPair)

set(gca,'ylim',GroupPair.LimY,'ytick',GroupPair.LimY);



end
statsSpeed(1,1)
statsSpeed(1,3)

figure;
clear dataAll dataNonT dataAll_sub statsSpeed
for iWhisk=1:length(ProcessPar.VolOut)

tempNonT=TBLallNontarget(TBLallNontarget.Sensory==ProcessPar.VolOut(iWhisk)&TBLallNontarget.PowerZero==0&TBLallNontarget.Speed<SpeedTh,:);
tempAll=TBLall(TBLall.Sensory==ProcessPar.VolOut(iWhisk)&TBLall.PowerZero==0&TBLall.Speed<SpeedTh,:);
tempAll_sub=TBLall_subNontarget(TBLall_subNontarget.Sensory==ProcessPar.VolOut(iWhisk)&TBLall_subNontarget.PowerZero==0&TBLall_subNontarget.Speed<SpeedTh,:);
% 
tempNonT=TBLallNontarget(TBLallNontarget.Sensory==ProcessPar.VolOut(iWhisk)&TBLallNontarget.PowerZero==0&abs(TBLallNontarget.Speed)<SpeedTh,:);
tempAll=TBLall(TBLall.Sensory==ProcessPar.VolOut(iWhisk)&TBLall.PowerZero==0&abs(TBLall.Speed)<SpeedTh,:);
tempAll_sub=TBLall_subNontarget(TBLall_subNontarget.Sensory==ProcessPar.VolOut(iWhisk)&TBLall_subNontarget.PowerZero==0&abs(TBLall_subNontarget.Speed<SpeedTh),:);


for iGroup=1:3
    dataAll{iGroup,iWhisk}=tempAll.SpeedScore(tempAll.Group==iGroup);
    dataNonT{iGroup,iWhisk}=tempNonT.SpeedScore(tempNonT.Group==iGroup);
    dataAll_sub{iGroup,iWhisk}=tempAll_sub.SpeedScore(tempAll_sub.Group==iGroup);

end

subplotLU(2,3,iWhisk,1)
statsSpeed(iWhisk,1)=ErrorBoxPlotLU(1:3,dataAll(:,iWhisk),ProcessPar.GroupColor,[],GroupPair)
set(gca,'ylim',GroupPair.LimY,'ytick',GroupPair.LimY);
subplotLU(2,3,iWhisk,2)
statsSpeed(iWhisk,2)=ErrorBoxPlotLU(1:3,dataNonT(:,iWhisk),ProcessPar.GroupColor,[],GroupPair)
set(gca,'ylim',GroupPair.LimY,'ytick',GroupPair.LimY);
subplotLU(2,3,iWhisk,3)
statsSpeed(iWhisk,3)=ErrorBoxPlotLU(1:3,dataAll_sub(:,iWhisk),ProcessPar.GroupColor,[],GroupPair)
set(gca,'ylim',GroupPair.LimY,'ytick',GroupPair.LimY);



end
statsSpeed(1,3)


   GroupPair.SignY=0.003;
   GroupPair.Plot=1;
   GroupPair.Std=0;      %%%%%%%%%using standard deviation as errorbar
   GroupPair.SamplePlot=1; %%%%%%%%%Plot Individual Sample Point
   GroupPair.SamplePairedPlot=0; %%%%%%%%%Dash line for paired comparison sample
   GroupPair.LimY=[-0.001 GroupPair.SignY*1.2];

figure;
clear dataAll dataNonT dataAll_sub
for iWhisk=1:length(ProcessPar.VolOut)

tempNonT=TBLallNontarget(TBLallNontarget.Sensory==ProcessPar.VolOut(iWhisk)&TBLallNontarget.PowerZero==0&TBLallNontarget.Speed<SpeedTh,:);
tempAll=TBLall(TBLall.Sensory==ProcessPar.VolOut(iWhisk)&TBLall.PowerZero==0&TBLall.Speed<SpeedTh,:);
tempAll_sub=TBLall_subNontarget(TBLall_subNontarget.Sensory==ProcessPar.VolOut(iWhisk)&TBLall_subNontarget.PowerZero==0&TBLall_subNontarget.Speed<SpeedTh,:);




for iGroup=1:3
    dataAll{iGroup,iWhisk}=tempAll.SensoryerScore(tempAll.Group==iGroup);
    dataNonT{iGroup,iWhisk}=tempNonT.SensoryerScore(tempNonT.Group==iGroup);
    dataAll_sub{iGroup,iWhisk}=tempAll_sub.SensoryerScore(tempAll_sub.Group==iGroup);

end

subplotLU(2,3,iWhisk,1)
statsWhisker(iWhisk,1)=ErrorBoxPlotLU(1:3,dataAll(:,iWhisk),ProcessPar.GroupColor,[],GroupPair)
set(gca,'ylim',[-0.001 0.005],'ytick',[-0.001:0.001:0.004]);

subplotLU(2,3,iWhisk,2)
statsWhisker(iWhisk,2)=ErrorBoxPlotLU(1:3,dataNonT(:,iWhisk),ProcessPar.GroupColor,[],GroupPair)
set(gca,'ylim',[-0.001 0.005],'ytick',[-0.001:0.001:0.004]);

subplotLU(2,3,iWhisk,3)
statsWhisker(iWhisk,3)=ErrorBoxPlotLU(1:3,dataAll_sub(:,iWhisk),ProcessPar.GroupColor,[],GroupPair)
set(gca,'ylim',[-0.001 0.005],'ytick',[-0.001:0.001:0.004]);


end



SpeedTh=1;
A=TBLallNontarget(TBLallNontarget.Sensory==1&TBLallNontarget.PowerZero==0&abs(TBLallNontarget.Speed)<SpeedTh,:);
B=TBLall(TBLall.Sensory==1&TBLall.PowerZero==0&abs(TBLall.Speed)<SpeedTh,:);

for iGroup=1:3
    tempData{iGroup}=A.WhiskerScore(A.Group==iGroup);
    tempData2{iGroup}=B.WhiskerScore(B.Group==iGroup);
end

   GroupPair.CorrName='fdr';
   GroupPair.Q=0.1;
   GroupPair.Pair=[];
   GroupPair.SignY=2;
   GroupPair.Plot=1;
   GroupPair.Std=0;      %%%%%%%%%using standard deviation as errorbar
   GroupPair.SamplePlot=1; %%%%%%%%%Plot Individual Sample Point
   GroupPair.SamplePairedPlot=0; %%%%%%%%%Dash line for paired comparison sample
   GroupPair.LimY=[-1 GroupPair.SignY*1.2];
   GroupPair.Marker={'o'};   % Datatype=varargin{1};  %%%%%%Datatype=0, non-paired data; 1 for paired data
   % SavePath=varargin{2};  %%%%%%Save the stats to txt file
   % GroupPair=varargin{3};
   % RGroupID=varargin{4};
   

% figure;
% stats=ErrorBoxPlotLU(1:3,tempData,ProcessPar.GroupColor,[])

figure;
subplot(1,2,1)
stats=ErrorBoxPlotLU(1:3,tempData2,ProcessPar.GroupColor,[])
subplot(1,2,2)
stats=ErrorBoxPlotLU(1:3,tempData,ProcessPar.GroupColor,[])


figure;
ErrorBarPlotLU(1:3,dataAll_sub(:,1),[],ProcessPar.GroupColor,2,0,[])

figure;
subplot(1,2,1)
stats=ErrorBarPlotLU(1:3,tempData2,[],ProcessPar.GroupColor,2,0,[],GroupPair)
subplot(1,2,2)
stats=ErrorBarPlotLU(1:3,tempData,[],ProcessPar.GroupColor,2,0,[],GroupPair)

stats=ErrorViolinHalf(1:3,tempData,ProcessPar.GroupColor,GroupPair)

x=1:3
GroupPair.CorrName='fdr';
GroupPair.Test='Ttest';
GroupPair.Q=0.1;
GroupPair.Pair=[];
GroupPair.SignY=2;
GroupPair.Plot=1;
GroupPair.Std=0;      %%%%%%%%%using standard deviation as errorbar
GroupPair.SamplePlot=1; %%%%%%%%%Plot Individual Sample Point
GroupPair.SamplePairedPlot=0; %%%%%%%%%Dash line for paired comparison sample
GroupPair.LimY=[0 GroupPair.SignY*1.2];
GroupPair.Marker={'o'};
GroupPair.ViolinLR=zeros(1,length(x));
GroupPairSpeed=GroupPair;
GroupPairSpeed.SignY=3;
GroupPairSpeed.LimY=[-1 GroupPairSpeed.SignY*1.2];

figure;
stats=ErrorViolinHalf(1:3,tempData2,ProcessPar.GroupColor,0,[],GroupPair,[1 2 3])

figure;
stats=ErrorViolinHalf(1:3,tempData,ProcessPar.GroupColor,0,[],GroupPair,[1 2 3])
set(gca,'ylim',[-3 3])


iData=1;
lowDimN=20;
ScoreLabel={'SpeedScore','WhiskScore'}
DataTable=[];
for iFOV=1:length(IndexFOVNeed)
tblFOV=SLMGroupTableTrial(SLMGroupTableTrial.Session==iFOV,:);
TrialID=unique(tblFOV.TrialID);
tbltemp=SLMGroupTableTrial(SLMGroupTableTrial.Session==iFOV&SLMGroupTableTrial.TrialID==TrialID(1),:);
% [W,H] = nnmf(NeuroTrace{iFOV}{iData}(tbltemp.TargetCell==0,SponInd),lowDimN);

I1=find(tblFOV.TrialID==TrialID(1));
NonTarget=find(tblFOV.TargetCell(I1)==0);

StructToVars(FOVres(iFOV));
pcaIspeed=ComI(1);
pcaIstim=ComI(2);



tblFOV=SLMGroupTableTrial(SLMGroupTableTrial.Session==iFOV,:);
TrialID=unique(tblFOV.TrialID);
WSub=W(:,[pcaIspeed pcaIstim]);
WNontargetSub=WNontarget(:,[pcaIspeed pcaIstim]);




Score=DataVec'*WSub;
ScoreNontarget=DataVecNonTarget'*WNontargetSub;


InfoTable(:,2:end)




end






% close all
% for PCI=1:10
%     figure;
%     [~,rankI]=sort((W(:,PCI)),'descend');
%     imagesc(NeuroTrace{iFOV}{iData}(rankI,SponInd));
%     % subplot(1,2,1)
%     % plot(AmpNormalize(coeff(SponInd,PCI),[1 99]),'k');
%     % hold on;plot(AmpNormalize(nanmean(BehTrace(iFOV).Speed(SponInd,:),2)),'g')
%     % % plot(AmpNormalize(coeff(SponInd,1),[1 99]),'r');
%     % hold on;plot(AmpNormalize(nanmean(BehTrace(iFOV).Stim(SponInd,:),2)),'r')
% end




close all
% figure;
for PCI=6:11
    figure;
    % subplot(1,2,1)
    plot(AmpNormalize(H(SponInd,PCI),[1 99]),'k');
    hold on;plot(AmpNormalize(nanmean(BehTrace(iFOV).Speed(SponInd,:),2)),'g')
    % plot(AmpNormalize(coeff(SponInd,1),[1 99]),'r');
    hold on;plot(AmpNormalize(nanmean(BehTrace(iFOV).Stim(SponInd,:),2)),'r')

end


figure;
plot(lags,c(:,10))
figure;
plot(lags,c(:,2))






close all
for PCI=1:5
    figure;
    [~,rankI]=sort((score(:,PCI)),'descend');
    imagesc(NeuroTrace{iFOV}{iData}(rankI,SponInd));
    % subplot(1,2,1)
    % plot(AmpNormalize(coeff(SponInd,PCI),[1 99]),'k');
    % hold on;plot(AmpNormalize(nanmean(BehTrace(iFOV).Speed(SponInd,:),2)),'g')
    % % plot(AmpNormalize(coeff(SponInd,1),[1 99]),'r');
    % hold on;plot(AmpNormalize(nanmean(BehTrace(iFOV).Stim(SponInd,:),2)),'r')
end
[coeff,score,latent,tsquared,explained,mu] = pca(NeuroTrace{iFOV}{iData}(:,SponInd));
explained(1:10);


lowDimN=10;
rSpeed=corr(coeff(ExcludeStimInd,1:lowDimN),nanmean(BehTrace(iFOV).Speed(ExcludeStimInd,:),2),'rows','complete')

maxLag=10;
clear rStim c
for PCI = 1:lowDimN
    [c(:,PCI), lags] = xcorr(coeff(SponInd,PCI), AmpNormalize(nanmean(BehTrace(iFOV).Stim(SponInd,:),2),[0 100]), maxLag, 'coeff');
    PostI = find(lags >= 0);
    PreI=find(lags<0);
    [~, i1] = max(abs(c(PostI,PCI))-mean(c(PreI,PCI)));
    rStim(PCI, 1) = c(PostI(i1),PCI);
end

[~,pcaIstim]=max(rStim)

[rankSpeedI,pcaIstim]=sort(rSpeed,'descend')

close all;
for i=1:12
figure;
plot(lags,c(:,i))
end
figure;
% plot(lags,c(:,1))
plot(lags,c(:,:))

close all
% figure;
for PCI=1:11
    figure;
    % subplot(1,2,1)
    plot(AmpNormalize(H(SponInd,PCI),[1 99]),'k');
    hold on;plot(AmpNormalize(nanmean(BehTrace(iFOV).Speed(SponInd,:),2)),'g')
    % plot(AmpNormalize(coeff(SponInd,1),[1 99]),'r');
    hold on;plot(AmpNormalize(nanmean(BehTrace(iFOV).Stim(SponInd,:),2)),'r')

end