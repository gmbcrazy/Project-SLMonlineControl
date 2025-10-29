clear all


BatchSavePath='D:\Project1-LocalProcessing\Step1\';
load([BatchSavePath '07-Oct-2025FOV.mat'])
Suite2pDataKeywords='awakeRefSpon';
DataSavePath='\\nimhlabstore1.nimh.nih.gov\UFNC\FNC2\Zhang\Projects\Project-LocalProcessing\Step3\';
mkdir(DataSavePath);
DataSavePath=[DataSavePath Suite2pDataKeywords '\'];
mkdir(DataSavePath);
SaveFunDate=[DataSavePath 'GroupSLM' date '\'];
mkdir(SaveFunDate)

VolEventLabel={}

PSTHparam.PreSLMCal = 15; 
PSTHparam.PostSLMCal = 15;
PSTHparam.pTh = 0.05; 
PSTHparam.TestMethod = 'ranksum';
PSTHparam.MPFrameJump = 2;
PSTHparam.TestStepFrame = 4;    %%post-slm frames for Test whether SLM works
PSTHparam.iData = 1;    %%post-slm frames for Test whether SLM works
% PSTHGroupSLM.TestStepFrame = 3;    %%post-slm frames for Test whether SLM works
% PSTHGroupSLM.PostSLMCal = 15;
PSTHparam.PreTestFrame = 13;    %%post-slm frames for Test whether SLM works
umPerlPixel=700/512;
offTargetum=15;
% PSTHparam.OffTargetPixel=ceil(offTargetum/umPerlPixel);
PSTHparam.TestWinNum=3;
PSTHparam.OffTargetUm=offTargetum;


ResponseMap=slanCM('seismic',64);
TimBinFrame = -PSTHparam.PreSLMCal:PSTHparam.PostSLMCal-1;

% IndexFOVNeed=[1 2 3 4 7 8 9 10 11 12 13 14 15 16 ];
% IndexFOVNeed=[1 2 3 4 7 8 9 10 11 12];
IndexFOVNeed=1:length(FOVUpdate)
% IndexFOVNeed=1:2



% Output=OfflineSLM_ExtractFOVs(FOVUpdate(IndexFOVNeed), Suite2pDataKeywords,suite2pFOVPathLocal(IndexFOVNeed),PSTHparam);
% [Output,NeuroTrace,BehTrace]=OfflineSLM_ExtractFOVs(FOVUpdate(IndexFOVNeed), Suite2pDataKeywords,suite2pFOVPath(IndexFOVNeed),PSTHparam);
% [Output,NeuroTrace,BehTrace]=OfflineSLM_ExtractFOVs(FOVUpdate(IndexFOVNeed), Suite2pDataKeywords,suite2pFOVPath(IndexFOVNeed),PSTHparam);

for iData=1:2
    PSTHparam.iData = iData;    %%post-slm frames for Test whether SLM works
    [Output(iData),NeuroTrace,BehTrace]=OfflineSLM_ExtractFOVs(FOVUpdate(IndexFOVNeed), Suite2pDataKeywords,suite2pFOVPath(IndexFOVNeed),PSTHparam);
end


save([SaveFunDate 'FOVoutputs.mat'],'-v7.3');
load([SaveFunDate 'FOVoutputs.mat'])


clear SLMSucc;
for iFOV=1:length(Output(1).GroupTargetCellMeta)
    FinalTargetCell=[];
    for igroup=1:length(Output(1).GroupTargetCellMeta{iFOV})
        FinalTargetCell=[FinalTargetCell(:);Output(1).GroupTargetCellMeta{iFOV}{igroup}(:)];
        GroupNumFOV(iFOV,igroup)=sum(Output(1).ActCellFunGroupFOV{iFOV}==igroup);
    end
    SLMSucc(iFOV,1)=length(Output(1).ActCellFunGroupFOV{iFOV})/length(FinalTargetCell);
end
NumCellSucc=sum(GroupNumFOV,2);
MaxGroupNum=max(GroupNumFOV,[],2)
NPerGrouTh=5;
NPerSessTh=12;
SessTh=0.7;

ValidInd=find(NumCellSucc>=NPerSessTh&SLMSucc>=SessTh&MaxGroupNum>=NPerGrouTh);
ValidInd=setdiff(ValidInd,[5 6 20]);  %%Session 5, 6 are from SL0886, WhiskerTrigger Runnning; Session 20, PV crash and FOV shift between Power test and SLM Group experiment
ValidSG=[];
for jFOV=1:length(ValidInd)
    for iGroup=1:size(GroupNumFOV,2)
        if GroupNumFOV((ValidInd(jFOV)),iGroup)>=NPerGrouTh
           ValidSG=[ValidSG;ValidInd(jFOV) iGroup];
        end
    end
end

ValidSG=array2table(ValidSG,'VariableNames',{'Session','Group'});

save([SaveFunDate 'ValidSessionGroup.mat'],'ValidSG','SLMSucc','NumCellSucc','MaxGroupNum','NPerGrouTh','NPerSessTh','SessTh');

load([SaveFunDate 'ValidSessionGroup.mat']);



PostPreDiffSpeedTh=[1 2 10000];
iSpeedTh=2;
ProcessPar.PostPreDiffSpeedTh=PostPreDiffSpeedTh(iSpeedTh);    %%noticed that this speed threshold is not applied when use OfflineSLM_FOVmeta2NeuroDeltaTrial.
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
%%Making table for mixed linear models with cell and trial as sample
NDataName={'deltaF','spks'};



tblDist= OfflineSLM_FOVmeta2CellSLMGroupDist(Output); 
writetable(tblDist,[SaveFunDate 'SLMCellGroupTarDistWin' num2str(iWin) '.csv']);


for iData=1:2
    PSTHparam.iData=iData;
    tblResponse=[];
    for iWin=1:PSTHparam.TestWinNum
        PSTHparam.PostWinN=iWin;
        [tbl,InfoTable]=OfflineSLM_FOVmeta2NeuroDeltaTrial(Output(iData),ProcessPar,PSTHparam);
    
        % SLMGroupTableTrial=[tbl InfoTable];
        NeedInfo={'Group','VolOut','PowerZero','AwakeState'};
    
        SLMGroupTableTrial = [tbl InfoTable(:,NeedInfo)];
        SLMGroupTableTrial(isnan(SLMGroupTableTrial.Group),:)=[];
        SLMGroupTableTrial=removevars(SLMGroupTableTrial,{'PointSpeedR','PointSensoryR'});
        
        SLMGroupTableTrial.Properties.VariableNames{'VolOut'} = 'Sensory';
        
        SLMGroupTrialID=SLMGroupTableTrial.TrialID;
        EndSessInd=[find(diff(SLMGroupTrialID)<0);length(SLMGroupTrialID)];
        StartSessInd=[1;find(diff(SLMGroupTrialID)<0)+1];
        Time=[];
        for iT=1:length(EndSessInd)
            TempSeq=StartSessInd(iT):EndSessInd(iT);
            TempSeq=(TempSeq-min(TempSeq))/(max(TempSeq)-min(TempSeq));
            Time=[Time;TempSeq(:)];
        end
        SLMGroupTableTrial.Time=Time;
        SLMGroupTableTrial=innerjoin(SLMGroupTableTrial,tblDist,"Keys",{'Cell','Group'});


        writetable(SLMGroupTableTrial,[SaveFunDate NDataName{iData} 'SLMGroupResTrialDynWin' num2str(iWin) '.csv']);


        for iSpeedTh=1:length(PostPreDiffSpeedTh)
            SaveP1=[SaveFunDate num2str(PostPreDiffSpeedTh(iSpeedTh)) '\'] ;
            mkdir(SaveP1);
            SLMtemp=SLMGroupTableTrial(SLMGroupTableTrial.Speed<=PostPreDiffSpeedTh(iSpeedTh),:);             %%Apply speed threshold to all trials
            AveTable=TblTrial2AveTbl(SLMtemp);
            writetable(AveTable,[SaveP1 NDataName{iData} 'SLMGroupRes_Win' num2str(iWin) '.csv']);

            s1=SLMGroupTableTrial.Speed<=PostPreDiffSpeedTh(iSpeedTh)&SLMGroupTableTrial.PowerZero==0;
            s2=SLMGroupTableTrial.PowerZero==1;
            SLMtemp=SLMGroupTableTrial(s1|s2,:); %%Apply speed threshold to SLMpower>0 trials but not SLMpower=0 trial
            % SLMtemp=SLMGroupTableTrial(SLMGroupTableTrial.Speed<=PostPreDiffSpeedTh(iSpeedTh),:);             %%Apply speed threshold to all trials
            AveTable=TblTrial2AveTbl(SLMtemp);
            AveTable= groupsummaryBack2OldNames(SLMtemp, AveTable, 'mean');

            writetable(AveTable,[SaveP1 NDataName{iData} 'SLMGroupRes_PowerZeroNoSpeedThWin' num2str(iWin) '.csv']);
         

        end

    end

end


for iData=1:2
    PSTHparam.iData=iData;
    tblResponse=[];
    for iWin=1:PSTHparam.TestWinNum

        SLMGroupTableTrial=readtable([SaveFunDate NDataName{iData} 'SLMGroupResTrialDynWin' num2str(iWin) '.csv']);


        for iSpeedTh=1:length(PostPreDiffSpeedTh)
            SaveP1=[SaveFunDate num2str(PostPreDiffSpeedTh(iSpeedTh)) '\'] ;
            mkdir(SaveP1);
            SLMtemp=SLMGroupTableTrial(SLMGroupTableTrial.Speed<=PostPreDiffSpeedTh(iSpeedTh),:);             %%Apply speed threshold to all trials
            AveTable=TblTrial2AveTbl(SLMtemp);
            AveTable= groupsummaryBack2OldNames(SLMtemp, AveTable, 'mean');
            writetable(AveTable,[SaveP1 NDataName{iData} 'SLMGroupRes_Win' num2str(iWin) '.csv']);

            s1=SLMGroupTableTrial.Speed<=PostPreDiffSpeedTh(iSpeedTh)&SLMGroupTableTrial.PowerZero==0;
            s2=SLMGroupTableTrial.PowerZero==1;
            SLMtemp=SLMGroupTableTrial(s1|s2,:); %%Apply speed threshold to SLMpower>0 trials but not SLMpower=0 trial
            % SLMtemp=SLMGroupTableTrial(SLMGroupTableTrial.Speed<=PostPreDiffSpeedTh(iSpeedTh),:);             %%Apply speed threshold to all trials
            AveTable=TblTrial2AveTbl(SLMtemp);
            AveTable= groupsummaryBack2OldNames(SLMtemp, AveTable, 'mean');

            writetable(AveTable,[SaveP1 NDataName{iData} 'SLMGroupRes_PowerZeroNoSpeedThWin' num2str(iWin) '.csv']);
         

        end

    end

end


PSTHparam.PostWinN=1;



% CellN=size(Output.NeuroPos3DMeta,1);
% CellNNonTarget=size(Output.NeuroPos3DMeta,1)-size(Output.GroupTargetCellAll,1);


% SaveFunDate=[DataSavePath 'GroupSLM' num2str(length(IndexFOVNeed)) 'Sessions\'];

PostPreDiffSpeedTh=[1 2 10000];

% PostPreDiffSpeedTh=2;


% tbl=readtable(['\\nimhlabstore1.nimh.nih.gov\UFNC\FNC2\Zhang\Projects\Project-LocalProcessing\Step3\awakeRefSpon\GroupSLM17-Oct-2025\1\deltaFSLMGroupRes_Win1.csv']);
% Keep only rows whose (Session, Group) are listed in ValidSG
% tbl_keep = innerjoin(tbl, unique(ValidSG(:,{'Session','Group'})), ...
%                      'Keys', {'Session','Group'});

PSTHparam.PostWinN=1;
%%Plot distributions
for iSpeedTh=1:length(PostPreDiffSpeedTh)
SaveP1=[SaveFunDate num2str(PostPreDiffSpeedTh(iSpeedTh)) '\'] ;
tbl=readtable([SaveP1 'deltaFSLMGroupRes_Win1.csv']);


mkdir(SaveP1);

ProcessPar.PostPreDiffSpeedTh=PostPreDiffSpeedTh(iSpeedTh);
ProcessPar.PlotFigure=0;
ProcessPar.GroupLabel={'L','S','N'};
ProcessPar.GroupList=[1 2 3];
ProcessPar.GroupColor=[255 51 153;91 20 212;121 247 111]/255;
ProcessPar.PowerZeroColor=[0.5 0.5 0.5];
ProcessPar.PowerZero=[0 1];
ProcessPar.PowerZeroLabel={'SLM','FakeSLM'};
ProcessPar.VolOut=[0 1];
ProcessPar.VolOutLabel={'NoSensory','Sensory'}
ProcessPar.AwakeState=[1];
ProcessPar.AwakeStateLabel={'Awake'};
% % GroupMetaName=[ProcessPar.GroupLabel {'FakeSLM'}];
% % GroupMetaColor=[ProcessPar.GroupColor;PowerZeroColor];
% % 

[OutResult,CellFOV]=OfflineSLM_FOVmeta2NeuroDelta(Output(iData),ProcessPar,PSTHparam);
% NonTargetCell=setdiff(1:size(Output.rSpeed,1),Output.GroupTargetCellAll(:,1));
% TargetCell=Output.GroupTargetCellAll(:,1);


% data1=OutResult(1).delta(:,4);
% data2=Output.rSpeed(:,1,1);
   Param.Color=ProcessPar.GroupColor;
   Param.Marker='o';
   Param.MarkerSize=12;
   Param.Rtype='pearson';
   % Param.xLim=[min(data1) max(data1)];
   % Param.yLim=[min(data2) max(data2)];
   Param.xLabel=[];
   Param.yLabel=[];
   Param.ExcludeColor=[0.5 0.5 0.5];
   Param.PlotExclude=0;
% [~,r,p]=LuPairRegressPlotGroup_ExcludeDots(data1,data2,Output.NeuroPos3DMeta(:,4),Output.GroupTargetCellAll(:,1), Param)    







% figure;
% for iWisk=1:2
%     clear TempData;
%     if iWisk==1
%        [~,rankI]=sort(Output(1).rSpeed(:,1,1),'descend');
%     else
%        [~,rankI]=sort(Output(1).rStim(:,1,1),'descend');
%     end
%     for iGroup=1:4
%         subplotLU(2,4,iWisk,iGroup)
%         imagesc(OutResult(iWisk).GroupResponse(rankI,:,iGroup));colormap(ResponseMap);set(gca,'clim',[-0.2 0.2])
%     end
% 
% end
% close all
   P1=[1 1 3;3 5 5];
   P1=[P1 P1+1];
   P2=[1 3 5;2 4 6];
   GroupPair.Pair=[];
   GroupPair.SignY=0.12;
   GroupPair.Plot=1;
   GroupPair.Std=0;      %%%%%%%%%using standard deviation as errorbar
   GroupPair.SamplePlot=1; %%%%%%%%%Plot Individual Sample Point
   GroupPair.SamplePairedPlot=1; %%%%%%%%%Dash line for paired comparison sample
   GroupPair.LimY=[0 GroupPair.SignY*1.2];
   GroupPair.Marker={'o'};
   GroupPair.ViolinLR=[0 0 0 0];
   GroupPair.GroupName=ProcessPar.GroupLabel;
   GroupPair.Q=0.1;

   GroupPair.CorrName='fdr';
   GroupPair.SignY=0.12;
   GroupPair.Plot=1;
   GroupPair.Std=0;      %%%%%%%%%using standard deviation as errorbar
   GroupPair.SamplePlot=1; %%%%%%%%%Plot Individual Sample Point
   GroupPair.SamplePairedPlot=1; %%%%%%%%%Dash line for paired comparison sample
   GroupPair.LimY=[0 GroupPair.SignY*1.2];
   GroupPair.Marker={'o'};
   GroupPair.ViolinLR=[0 0 0 0];
   GroupPair.Test='Ttest';
   Param.Color=[1.0000    0.2000    0.6000;0.3569    0.0784    0.8314;0.4745    0.9686    0.4353;0.6 0.6 0.6];

   GroupPair.GroupName={'L','S','N','0'}

for iWisk=1:2

    SaveP2=[SaveP1 ProcessPar.VolOutLabel{iWisk} '\'] ;
    mkdir(SaveP2);

    tbl_whisktemp = tbl(tbl.Sensory==ProcessPar.VolOut(iWisk),:);
    
    tempCellID=repmat([1:size(OutResult(iWisk).p,1)]',1,4);
    tempFOVID=repmat(CellFOV(:,2),1,4);
    tempFunID=repmat(1:4,size(CellFOV,1),1);
    PvalueTBL=array2table([OutResult(iWisk).p(:) tempCellID(:) tempFunID(:) tempFOVID(:)],'VariableNames',{'Pvalue','Cell','Group','Session'});
    % PvalueTBLShamOpto=array2table(PvalueTBL.Pvalue(PvalueTBL.Group==4),'VariableNames',{'ShamOptoPvalue'});
    PvalueTBLShamOpto=PvalueTBL(PvalueTBL.Group==4,:);
    PvalueTBLShamOpto.Properties.VariableNames{'Pvalue'}='ShamOptoPvalue';
    PvalueTBLShamOpto(:,{'Group'})=[];
    PvalueTBLShamOpto(:,{'Session'})=[];
    PvalueTBL=PvalueTBL(PvalueTBL.Group~=4,:);
    PvalueTBL=innerjoin(PvalueTBL,PvalueTBLShamOpto,'Keys', {'Cell'});




    PvalueTBL=innerjoin(PvalueTBL, unique(ValidSG(:,{'Session','Group'})), ...
                     'Keys', {'Session','Group'});
    tbl_whisktemp=innerjoin(tbl_whisktemp , unique(ValidSG(:,{'Session','Group'})), ...
                     'Keys', {'Session','Group'});




    if size(PvalueTBL,1)~=size(tbl_whisktemp,1)
       disp(['Individual Cell Pvalue table Do NOT match tbl response table size!']);
    else
       tbl_whisktemp=innerjoin(tbl_whisktemp,PvalueTBL,'Keys',{'Cell','Group'});
    end



    % NonTargetCelltbl=tbl_whisktemp(tbl_whisktemp.TargetCell==0,:);
    % subplot(1,2,iWisk)
    figure;
    subplot(1,2,1)
    title('Target Cell')
    TempData={};
    GroupPair.LimY=[0 0.8];
    GroupPair.SignY=0.75;
    
    for iGroup=1:3
        TempData{iGroup}=tbl_whisktemp.Response(tbl_whisktemp.Group==iGroup&tbl_whisktemp.TargetCell==iGroup);
    % [~,pTFromZero(iWisk,iGroup),~,tFromZero(iWisk,iGroup)]=ttest(TempData{iGroup},0);
    % [pRFromZero(iWisk,iGroup),~,RFromZero(iWisk,iGroup)]=signrank(TempData{iGroup},0);
    end
    % TempData{iGroup+1}=a.ShamOpto(a.Group==iGroup);  %%Noted that in the table, ShamOpto value is same for different group.
    % stats=ErrorViolinHalf([1 2 3],TempData(1:3),Param.Color,1,[SaveP1 ProcessPar.VolOutLabel{iWisk} '.txt'],GroupPair,[1 1 1]);
    stats=ErrorViolinHalf([1 2 3],TempData(1:3),[Param.Color],0,[SaveP2 'TargetResponse.txt'],GroupPair,[1 2 3]);
    ylabel('Cell response (Post-Pre ΔF)');
    xlabel('Stim group')
    hold on;
    plot([0 4],[0 0],'k:')
    set(gca,'ylim',[-0.05 0.8],'xlim',[0 4],'xtick',1:3,'XTickLabel',GroupPair.GroupName)
    papersizePX=[0 0 12 12];
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));

    subplot(1,2,2)
    title('Non-target Cell')
    TempData={};
    for iGroup=1:3
        TempData{iGroup}=tbl_whisktemp.Response(tbl_whisktemp.Group==iGroup&tbl_whisktemp.NonTargetCell==1);
    % [~,pTFromZero(iWisk,iGroup),~,tFromZero(iWisk,iGroup)]=ttest(TempData{iGroup},0);
    % [pRFromZero(iWisk,iGroup),~,RFromZero(iWisk,iGroup)]=signrank(TempData{iGroup},0);
    end
    % TempData{iGroup+1}=a.ShamOpto(a.Group==iGroup);  %%Noted that in the table, ShamOpto value is same for different group.
    % stats=ErrorViolinHalf([1 2 3],TempData(1:3),Param.Color,1,[SaveP1 ProcessPar.VolOutLabel{iWisk} '.txt'],GroupPair,[1 1 1]);
    stats=ErrorViolinHalf([1 2 3],TempData(1:3),[Param.Color],0,[SaveP2 'NonTargetResponse.txt'],GroupPair,[1 2 3]);
    ylabel('Cell response (Post-Pre ΔF)');
    xlabel('Stim group')
    hold on;
    plot([0 4],[0 0],'k:')
    set(gca,'ylim',[-0.1 0.2],'xlim',[0 4],'xtick',1:3,'XTickLabel',GroupPair.GroupName)
    papersizePX=[0 0 20 12];
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
    % print(gcf, [SaveP2 'CellResponse.svg'], '-dsvg', '-painters');
    % print(gcf, [SaveP2 'CellResponse.tif'], '-dtiffn', '-painters');

    % tbl_whisktemp.Properties.VariableNames{'Session_tbl_whisktemp'}='Session';



    %%Session Analysis of Non-target cell
    NonTargetTBL=tbl_whisktemp(tbl_whisktemp.NonTargetCell==1,:);
    tbl_whisktemp(tbl_whisktemp.NonTargetCell==1,:);
    tbl_whisktemp.Properties.VariableNames{'Session_tbl_whisktemp'}='Session';
    NonTargetTBL=tbl_whisktemp(tbl_whisktemp.NonTargetCell==1,:);

    % NonTargetTBL.Properties.VariableNames{'Session_tbl_whisktemp'}='Session';




    KeyMerge={'Session','Group'};
    NonTargetY1=groupsummary(NonTargetTBL,KeyMerge, @(x,y) corr(x,y),{"Response","SpeedR"});
    NonTargetY1.Properties.VariableNames{'fun1_Response_SpeedR'}='Response_SpeedR';
    NonTargetY2=groupsummary(NonTargetTBL,KeyMerge, @(x,y) corr(x,y),{"Response","SensoryR"});
    NonTargetY2.Properties.VariableNames{'fun1_Response_SensoryR'}='Response_SensoryR';
    NonTargetY=innerjoin(NonTargetY1,NonTargetY2,"Keys",KeyMerge);
    % A=NonWhiskT(NonWhiskT.TargetCell==1,:);
    TargetX=groupsummary(NonTargetTBL,KeyMerge, 'mean');
    TargetX= groupsummaryBack2OldNames(NonTargetTBL, TargetX, 'mean');
    AveResTBL = innerjoin(NonTargetY,TargetX,"Keys",KeyMerge);


   Param.Color=ProcessPar.GroupColor;
   Param.Marker='o';
   Param.MarkerSize=12;
   Param.Rtype='pearson';
   % Param.xLim=[min(data1) max(data1)];
   % Param.yLim=[min(data2) max(data2)];
   Param.xLabel=[];
   Param.yLabel=[];
   Param.ExcludeColor=[0.5 0.5 0.5];
   Param.PlotExclude=0;
   Param.xLim=[-0.2 0.5];
   Param.yLim=[-0.2 0.5];
 
    Param.xLim=[-0.3 0.6];
    Param.yLim=[-0.3 0.6];
  

    FOVLevelCorrPlots(AveResTBL,Param,SaveP2)
    FOVLevelCorrPlots_GroupSubPlots(AveResTBL, Param, SaveP2)



    GroupParamNet.ScoreLim=[-0.3 0.3];
    GroupParamNet.ResponseLim=[-0.2 0.2];
    GroupParamNet.ScoreMap=slanCM('wildfire',64);
    GroupParamNet.ScoreLabel='Speed Corr.';
    GroupParamNet.ResponseMap=slanCM('seismic',64);
    % GroupParamNet.NodeColor=NodeColor;
    % GroupParamNet.statCellRes=statGroupRes;
    GroupParamNet.CirSubP=[0.02 0.1 0.15 0.12 0.06 0.02];
    GroupParamNet.SessionMap=slanCM('glasbey_light',32);
    GroupParamNet.GroupColor=ProcessPar.GroupColor;

    for iGroup=1:3
    [GgraphOut,tHandels]=tblMeta2CircosPlot(tbl_whisktemp, iGroup, GroupParamNet, 'SpeedR');
     papersizePX=[0 0 22 13];
     set(gcf, 'PaperUnits', 'centimeters');
     set(gcf, 'PaperPosition', papersizePX, 'PaperSize', papersizePX(3:4));
     print(gcf, [SaveP2 'Circosgroup' num2str(iGroup) 'ResposeVsSpeedR.tif'], '-dtiffn', '-painters');
     print(gcf, [SaveP2 'Circosgroup' num2str(iGroup) 'ResposeVsSpeedR.svg'], '-dsvg', '-painters');
     close all
        % print(gcf, [SaveP2 fileName '.svg'], '-dsvg', '-painters');
    [GgraphOut,tHandels]=tblMeta2CircosPlot(tbl_whisktemp, iGroup, GroupParamNet, 'SensoryR');
     papersizePX=[0 0 22 13];
     set(gcf, 'PaperUnits', 'centimeters');
     set(gcf, 'PaperPosition', papersizePX, 'PaperSize', papersizePX(3:4));
     print(gcf, [SaveP2 'Circosgroup' num2str(iGroup) 'ResposeVsSensoryR.tif'], '-dtiffn', '-painters');
     print(gcf, [SaveP2 'Circosgroup' num2str(iGroup) 'ResposeVsSensoryR.svg'], '-dsvg', '-painters');
     close all

    end

end
 




end


