
clear all


TBLDataPath='\\nimhlabstore1.nimh.nih.gov\UFNC\FNC2\Zhang\Projects\Project-LocalProcessing\Step3\awakeRefSpon\GroupSLM20-Nov-2025\';


load([TBLDataPath 'FOVoutputs.mat'])
load([TBLDataPath 'ValidSessionGroup.mat'])


SavePath=['\\nimhlabstore1.nimh.nih.gov\UFNC\FNC2\Zhang\Projects\Project-LocalProcessing\Step4\' Suite2pDataKeywords '\'];
mkdir(SavePath);
SavePath=[SavePath date '\'];
mkdir(SavePath);


% animalpaths={FOVUpdate.DataFolder};
% 
% for i = 1:length(animalpaths)
%     % Extract animal ID (e.g., SL0886)
%     animalMatch = regexp(animalpaths{i}, '(SL\d{4}|L\d{5})', 'match');
%     if ~isempty(animalMatch)
%         animalIDs{i} = animalMatch{1};
%     else
%         animalIDs{i} = '';
%     end
% 
%     % Extract date (8-digit number)
%     dateMatch = regexp(animalpaths{i}, '\d{8}', 'match');
%     if ~isempty(dateMatch)
%         sessionDates{i} = dateMatch{1};
%     else
%         sessionDates{i} = '';
%     end
% 
%     % suite2pFOVPath{i}=[Suite2pSaveFolderAll '\' animalIDs{i} sessionDates{i} '\'];
%     % mkdir(suite2pFOVPath{i})
%     % suite2pFOVPathLocal{i}=[Suite2pSaveFolderAllLocal '\' animalIDs{i} sessionDates{i} '\'];
%     % mkdir(suite2pFOVPathLocal{i})
% 
% end

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
SpeedPercentile=10;
% iData=1;
lowDimN=40;
ScoreLabel={'SpeedScore','WhiskScore'}
NDataName={'deltaF','spks'};

parpool('local', 4);

% for iFOV=3:3
for iData=1:length(NDataName)

    SaveFunFOVDimData=[SaveFunFOVDim NDataName{iData} '\'];
    mkdir(SaveFunFOVDimData);
    SaveFunSub=[SaveFunFOVDimData 'FOV\'];
    mkdir(SaveFunSub)
    SLMGroupTableTrial=readtable([TBLDataPath NDataName{iData} 'SLMGroupResTrialDynWin1.csv']);   %%Win# has no impact on the results. As training only use spontouse data
    IndexFOVNeed=unique(ValidSG.Session)
    clear FOVres
        parfor iFOV=1:length(IndexFOVNeed)
        % for iFOV=7:8

            % FOVres(iFOV) = FOVnnmf_SpeedWStim(IndexFOVNeed(iFOV), SLMGroupTableTrial, NeuroTrace, BehTrace, iData, lowDimN, SponInd, ExcludeStimInd);
            FOVres(iFOV) = FOVnnmf_SpeedWStim_SpeedSubset(IndexFOVNeed(iFOV), SLMGroupTableTrial, NeuroTrace, BehTrace, iData, lowDimN, SponInd, ExcludeStimInd, SpeedPercentile)
            FOVsavePath=[SaveFunSub 'FOV' num2str(IndexFOVNeed(iFOV)) animalIDs{IndexFOVNeed(iFOV)} sessionDates{IndexFOVNeed(iFOV)} ];
            plotFOVnnmf(FOVres(iFOV), SponInd, BehTrace(IndexFOVNeed(iFOV)), ProcessPar, FOVsavePath);
            close all
        end
    save([SaveFunFOVDimData 'nnmfResult.mat'],'-v7.3');
end









iData=1;
SaveFunFOVDimData=[SaveFunFOVDim NDataName{iData} '\'];
load([SaveFunFOVDimData 'nnmfResult.mat']);
iWin=1;
SLMGroupTableTrial=readtable([TBLDataPath NDataName{iData} 'SLMGroupResTrialDynWin1.csv']);   %%Win# has no impact on the results. As training only use spontouse data


DualScore=[];
DualScore_NonTarget=[];

DualGroup=[];
DualSession=[];
for iFOV=1:length(IndexFOVNeed)
    DualScore = [DualScore;FOVres(iFOV).SpeedScore FOVres(iFOV).SensoryScore];
    DualScore_NonTarget = [DualScore_NonTarget;FOVres(iFOV).SpeedScore_NonTarget FOVres(iFOV).SensoryScore_NonTarget];

    DualGroup = [DualGroup;1;2];
    DualSession=[DualSession;IndexFOVNeed(iFOV);IndexFOVNeed(iFOV);];
end



nnmfSG=array2table([DualSession DualGroup DualScore DualScore_NonTarget],"VariableNames",{'Session','Group','SpeedScore','SensoryScore','SpeedScore_Nontarget','SensoryScore_Nontarget'});



% ValidSG_DualScore(isnan(ValidSG_DualScore.SpeedScore),:)=[];
% 

figure;
subplot(1,2,1)
gscatterLU(DualScore(:,1), DualScore(:,2), DualGroup, ...
           ProcessPar.GroupColor(1:3,:), 'o', 12, 'filled');

subplot(1,2,2)
gscatterLU(DualScore_NonTarget(:,1), DualScore_NonTarget(:,2), DualGroup, ...
           ProcessPar.GroupColor(1:3,:), 'o', 12, 'filled');

HighSensorySpeedCut=0.3;
LowSensoryCut=0.05;
LowSpeedCut=0.4;


ValidSG_DualScore = outerjoin(ValidSG,nnmfSG, 'Keys', {'Session','Group'},'MergeKeys', true);


InvalidSensoryI=find((ValidSG_DualScore.SpeedScore_Nontarget>HighSensorySpeedCut|ValidSG_DualScore.SensoryScore_Nontarget<LowSensoryCut)&ValidSG_DualScore.Group==2);
InvalidSpeedI=find(ValidSG_DualScore.SpeedScore_Nontarget<LowSpeedCut&ValidSG_DualScore.Group==1);
InvalidI=union(InvalidSensoryI,InvalidSpeedI);
ValidI=setdiff(setdiff([1:size(ValidSG_DualScore,1)]',InvalidI),find(isnan(ValidSG_DualScore.SpeedScore)==1));

InValidnnmfSG=ValidSG_DualScore(InvalidI,:);
ValidnnmfSG=ValidSG_DualScore(ValidI,:);


TempValidSG_DualScore1=ValidSG_DualScore(ValidI,:)
TempValidSG_DualScore2=ValidSG_DualScore;
TempValidSG_DualScore2(InvalidI,:)=[];


MergedValidSG=TempValidSG_DualScore2;
MergedValidSG=MergedValidSG(MergedValidSG.GroupTargetCellN>=5,:);


SessionUni=unique(MergedValidSG.Session)
illG3=[];
for iS=1:length(SessionUni)
    I1=find(MergedValidSG.Session==SessionUni(iS));
    if length(I1)==1
       if MergedValidSG.Group(I1)==3
          illG3=[illG3;I1];
       end
    end
end
MergedValidSG (illG3,:)=[];
MergedValidSG = MergedValidSG (:,{'Session','Group','GroupTargetCellN'});



% % ValidSG_DualScore = MergedValidSG;
% % ValidSG_DualScore(isnan(ValidSG_DualScore.SpeedScore),:)=[];


% MergedValidSG = MergedValidSG(:,{'Session','Group','GroupTargetCellN'})

save([SaveFunFOVDimData 'ValidnnmfSG.mat'],'ValidSG_DualScore','MergedValidSG','ValidSG')


load([SaveFunFOVDimData 'ValidnnmfSG.mat'],'ValidSG_DualScore','MergedValidSG','ValidSG')



%%
figure;
gscatterLU(ValidSG_DualScore.SpeedScore_Nontarget(InvalidSensoryI), ValidSG_DualScore.SensoryScore_Nontarget(InvalidSensoryI), ValidSG_DualScore.Group(InvalidSensoryI), ProcessPar.GroupColor(2,:), 'x', 15);
hold on;
gscatterLU(ValidSG_DualScore.SpeedScore_Nontarget(InvalidSpeedI), ValidSG_DualScore.SensoryScore_Nontarget(InvalidSpeedI), ValidSG_DualScore.Group(InvalidSpeedI), ProcessPar.GroupColor(1,:), 'x', 15);
hold on;
gscatterLU(ValidSG_DualScore.SpeedScore_Nontarget(ValidI), ValidSG_DualScore.SensoryScore_Nontarget(ValidI), ValidSG_DualScore.Group(ValidI), ProcessPar.GroupColor, 'o', 15,'filled');
hold on;

% gscatterLU(DualScore_NonTarget(ValidI,1), DualScore_NonTarget(ValidI,2), DualGroup(ValidI), ProcessPar.GroupColor, 'o', 15,'filled');

% hold on;
% gscatterLU(DualScore_NonTarget(:,1), DualScore_NonTarget(:,2), DualGroup, ...
%            ProcessPar.GroupColor(1:2,:), 'o', 12, 'filled');

hold on;
plot([-0.2 0.8],[LowSensoryCut LowSensoryCut],':','Color',ProcessPar.GroupColor(2,:));
plot([HighSensorySpeedCut HighSensorySpeedCut],[-0.2 0.4],':','Color',ProcessPar.GroupColor(1,:));

legend({'Excluded sensory component','Excluded speed component','Speed component','Sensory component'})
xlabel('SpeedScore');
ylabel('SensoryScore');
set(gca,'xlim',[-0.2 0.8],'ylim',[-0.2 0.3])
papersizePX = [0 0 15 15];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', papersizePX, 'PaperSize', papersizePX(3:4));
print(gcf, [SaveFunFOVDimData 'NNMFscores.svg'], '-dsvg', '-painters');
print(gcf, [SaveFunFOVDimData 'NNMFscores.tif'], '-dtiffn', '-painters');



