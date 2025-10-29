clear all

SaveFunFOV=[DataSavePath 'GroupSLM' num2str(length(IndexFOVNeed)) 'Sessions\'];
mkdir(SaveFunFOV)



% save([SaveFunFOV 'FOVoutputs.mat'])
load([SaveFunFOV 'FOVoutputs.mat'])

load(['\\nimhlabstore1.nimh.nih.gov\UFNC\FNC2\Zhang\Projects\Project-LocalProcessing\Step3\awakeRefSpon\GroupSLM23-Oct-2025\FOVoutputs.mat'])



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

for iData=1:2

    PSTHparam.iData=iData;

PSTHparam.TestWinNum=5;
tblResponse=[];
for iWin=1:PSTHparam.TestWinNum
    PSTHparam.PostWinN=iWin;
    [tbl,InfoTable]=OfflineSLM_FOVmeta2NeuroDeltaTrial(Output(iData),ProcessPar,PSTHparam);

    % SLMGroupTableTrial=[tbl InfoTable];
    NeedInfo={'Group','VolOut','PowerZero','AwakeState'};
    SLMGroupTableTrial = [tbl InfoTable(:,NeedInfo)];
    SLMGroupTableTrial(isnan(SLMGroupTableTrial.Group),:)=[];
    SLMGroupTableTrial=removevars(SLMGroupTableTrial,{'PointSpeedR','PointStimR'});
    
    SLMGroupTableTrial.Properties.VariableNames{'VolOut'} = 'Whisk';
    
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
    writetable(SLMGroupTableTrial,[SaveFunCon NDataName{iData} 'SLMGroupResTrialDynWin' num2str(iWin) '.csv']);
    tblResponse(:,iWin)=tbl.Response;
end

end

PSTHparam.TestWinNum=1;
for iData=1:2

    PSTHparam.iData=iData;
    tblResponse=[];
    for iWin=1:PSTHparam.TestWinNum
        PSTHparam.PostWinN=iWin;
        [tbl,InfoTable]=OfflineSLM_FOVmeta2NeuroDeltaTrial(Output(iData),ProcessPar,PSTHparam);
        % SLMGroupTableTrial=[tbl InfoTable];
        NeedInfo={'UncagingLaserPower','Point','PointTargetCell','PointTargetCellGroup'};
    
        SLMPointTableTrial = [tbl InfoTable(:,NeedInfo)];
        SLMPointTableTrial(isnan(SLMPointTableTrial.Point),:)=[];
        % SLMPointTableTrial=removevars(SLMPointTableTrial,{'PointSpeedR','PointStimR'});
        
        % SLMPointTableTrial.Properties.VariableNames{'VolOut'} = 'Whisk';
        
        SLMPointTrialID=SLMPointTableTrial.TrialID;
        EndSessInd=[find(diff(SLMPointTrialID)<0);length(SLMPointTrialID)];
        StartSessInd=[1;find(diff(SLMPointTrialID)<0)+1];
        Time=[];
        for iT=1:length(EndSessInd)
            TempSeq=StartSessInd(iT):EndSessInd(iT);
            TempSeq=(TempSeq-min(TempSeq))/(max(TempSeq)-min(TempSeq));
            Time=[Time;TempSeq(:)];
        end
        SLMPointTableTrial.Time=Time;
        writetable(SLMPointTableTrial,[SaveFunCon NDataName{iData} 'SLMPointResTrialDynWin' num2str(iWin) '.csv']);
        tblResponse(:,iWin)=tbl.Response;
    end

end


iData=1;
SLMPointTableTrial=readtable([SaveFunFOV NDataName{iData} 'SLMPointResTrialDynWin1.csv']);

SaveFunCon=[SaveFunFOV 'FunCon\'];
mkdir(SaveFunCon)


DimentionReductionMethod='NonNegMatFac';
SaveFunConDim=[SaveFunCon DimentionReductionMethod '\'];
mkdir(SaveFunConDim)


SponInd=[1:8200];
ExcludeStimInd=[1:2000 2800:4100];  
ExcludeStimInd=union(ExcludeStimInd,ExcludeStimInd+4100);

% AddStimNoise=random('norm',0,0.1,length(SponInd),1);

SaveFunSub=[SaveFunConDim 'FOV\'];
mkdir(SaveFunSub)
% clear FOVres

SaveFunConDim
% iData=2;
% pcaN=40;
ScoreLabel={'SpeedScore','WhiskScore'}
activate_Pth=[0.2;0.1;0.05;0.02;0.01;0.005;0.001];
PthLabel={'20','10','05','02','01','005','001'};
MultiF_Pth=0.05;

% for iPth=1:3
for iPth=1:length(PthLabel)

SaveFunConPth=[SaveFunCon 'FC_pth' PthLabel{iPth} '\'];
mkdir(SaveFunConPth)

r1Total=[];
r2Total=[];
% for iFOV=1:3
rgroup=[];
ngroup=[];
tblFOVclean=[];   %%refers to remove trials with power test when the power is not strong enough to activate one target cell.


    for iFOV=1:length(IndexFOVNeed)
        close all
    tblFOV=SLMPointTableTrial(SLMPointTableTrial.Session==iFOV,:);
    tempNeuroPos3D=Output(1).NeuroPos3DumMeta(Output(1).NeuroPos3DMeta(:,4)==IndexFOVNeed(iFOV),1:3);
    % tempNeuroPos3D(:,1:3)=tempNeuroPos3D(:,1:3)*FOVscalePerPixel;



    tempTargetList=Output(iData).TargetCellListFOV{iFOV};
    tempPowerI=Output(iData).PowerTargetIFOV{iFOV};
    tempPower=unique(tblFOV.UncagingLaserPower);
    tempGroup=Output(iData).TargetCellListFunGroupFOV{iFOV};


    %Cells pass the tested might be more than the final included cell, this
    %part only inlucde the cell later used in SLM group stimuli
        finalTargetList=[];
        for igroup=1:length(Output(iData).GroupTargetCellMeta{iFOV})
            finalTargetList=[finalTargetList;Output(iData).GroupTargetCellMeta{iFOV}{igroup}(:)];
        end
        finalTargetList=sort(finalTargetList);
        i1=ismember(tempTargetList,finalTargetList);
    % tempTargetList=tempTargetList(i1);
    % tempPowerI=tempPowerI(i1);
    % tempGroup=tempGroup(i1);
    %Cells pass the tested might be more than the final included cell, this
    %part only inlucde the cell later used in SLM group stimuli

    %Only significant activated cell were considered
    validI=find(tempPowerI>0&i1>0);
    tempTargetList=tempTargetList(validI);
    tempPowerI=tempPowerI(validI);
    tempPower=tempPower(tempPowerI);
    tempGroup=tempGroup(validI);
    %Only significant activated cell were considered




    for iCell=1:length(validI)
    
        ActivateI=Output(iData).statCellRes{iFOV}(validI(iCell),tempPowerI(iCell)).p<activate_Pth(iPth);
        TrialID=unique(tblFOV.TrialID(tblFOV.PointTargetCell==tempTargetList(iCell)&tblFOV.UncagingLaserPower==tempPower(iCell)));
        tempData=[];
        for itrial=1:length(TrialID)
            tempData(:,itrial)=tblFOV.Response(tblFOV.TrialID==TrialID(itrial));
        end
    
        rgroup(end+1)=tempGroup(iCell);
    
        r1=corr(tempData);
        r1Total(end+1)=mean(r1(triu(true(size(r1)),1)));
        if sum(ActivateI)>5
           r2=corr(tempData(ActivateI,:));
           r2Total(end+1)=mean(r2(triu(true(size(r2)),1)));
        else
           r2Total(end+1)=NaN;
        end
    
        temptbl=tblFOV(ismember(tblFOV.TrialID,TrialID),:);
        temptbl.ActivateI=repmat(ActivateI(:),length(TrialID),1);
        temptbl.VecR=repmat(r1Total(end),size(temptbl,1),1);
        temptbl.ActVecR=repmat(r2Total(end),size(temptbl,1),1);
        

        % Coortbl=[repmat(tempNeuroPos3D,length(TrialID),1) repmat(tempNeuroPos3D(tempTargetList(iCell),:),size(temptbl,1),1)];

        Diff3Dtemp=tempNeuroPos3D-tempNeuroPos3D(tempTargetList(iCell),:);
        Dist3Dtemp=sqrt(sum(Diff3Dtemp.^2,2));
        Dist2Dtemp=sqrt(sum(Diff3Dtemp(:,1:2).^2,2));
        CoLayer=Diff3Dtemp(:,3)==0;

        spaMetrictemp= table(Dist3Dtemp(:),Dist2Dtemp(:),CoLayer(:),'VariableNames',{'Dist','DistXY','CoLayer'});
        
        Coortbl=[tempNeuroPos3D repmat(tempNeuroPos3D(tempTargetList(iCell),:),size(tempNeuroPos3D,1),1)];
        Coortbl=array2table(Coortbl,'VariableNames',{'X','Y','Z','PointX','PointY','PointZ'});

        temptbl=[temptbl repmat([spaMetrictemp Coortbl],length(TrialID),1)];
    
        tblFOVclean=[tblFOVclean;temptbl];
    
        
    end
    
    end




GroupMethod='mean';

illegalSession=[5 6];
tblFOVclean(ismember(tblFOVclean.Session,illegalSession),:)=[];




tblFOVclean.ActivateIPos=tblFOVclean.ActivateI.*tblFOVclean.Response>0;
tblFOVclean.ActivateINeg=tblFOVclean.ActivateI.*tblFOVclean.Response<0;

tblFOVclean.ActSpeedR=tblFOVclean.ActivateI.*tblFOVclean.SpeedR;
tblFOVclean.ActPosSpeedR=tblFOVclean.ActivateIPos.*tblFOVclean.SpeedR;
tblFOVclean.ActNegSpeedR=tblFOVclean.ActivateINeg.*tblFOVclean.SpeedR;


tblFOVclean.ActStimR=tblFOVclean.ActivateI.*tblFOVclean.StimR;
tblFOVclean.ActPosStimR=tblFOVclean.ActivateIPos.*tblFOVclean.StimR;
tblFOVclean.ActNegStimR=tblFOVclean.ActivateINeg.*tblFOVclean.StimR;

GroupMethod='mean';
newtbl1 = groupsummary(tblFOVclean, {'Session', 'Cell'}, GroupMethod);
newtbl1=groupsummaryBack2OldNames(tblFOVclean,newtbl1,GroupMethod);

[a,p]=corr(newtbl1.SpeedR, newtbl1.StimR,'type','Spearman');

% tblFOVcleanavg = groupsummary(tblFOVclean(tblFOVclean.Response>0,:), {'Session', 'PointTargetCell','Cell'}, GroupMethod);
tblFOVcleanavg = groupsummary(tblFOVclean, {'Session', 'PointTargetCell','Cell'}, GroupMethod);
newtbl2=groupsummaryBack2OldNames(tblFOVclean,tblFOVcleanavg,GroupMethod);

newtbl2.ActivateIPos=newtbl2.ActivateI.*newtbl2.Response>0;
newtbl2.ActivateINeg=newtbl2.ActivateI.*newtbl2.Response<0;

% newtbl2.ActSpeedR=newtbl2.ActivateI.*newtbl2.SpeedR;
% newtbl2.ActPosSpeedR=newtbl2.ActivateIPos.*newtbl2.SpeedR;
% newtbl2.ActNegSpeedR=newtbl2.ActivateINeg.*newtbl2.SpeedR;
% 
% 
% newtbl2.ActStimR=newtbl2.ActivateI.*newtbl2.StimR;
% newtbl2.ActPosStimR=newtbl2.ActivateIPos.*newtbl2.StimR;
% newtbl2.ActNegStimR=newtbl2.ActivateINeg.*newtbl2.StimR;
% 




newtbl3 = groupsummary(newtbl2, {'Session', 'PointTargetCell'}, GroupMethod);
newtbl3=groupsummaryBack2OldNames(newtbl2,newtbl3,GroupMethod);

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


AIabs_StimR   = splitapply(myfun, newtbl2.ActivateI,    newtbl2.StimR,  G);
AIP_StimR  = splitapply(myfun, newtbl2.ActivateIPos, newtbl2.StimR,  G);
AIN_StimR  = splitapply(myfun, newtbl2.ActivateINeg, newtbl2.StimR,  G);
AI_StimR =  splitapply(myfun2, newtbl2.ActivateIPos,newtbl2.ActivateINeg, newtbl2.StimR, G);
% Build result table using the grouping labels that match
newtbl4 = table( ...
    sessionVals, TargetcellVals, ...
    AIabs_SpeedR,AI_SpeedR, AIP_SpeedR, AIN_SpeedR, ...
    AIabs_StimR,AI_StimR, AIP_StimR, AIN_StimR, ...
    'VariableNames', {'Session','PointTargetCell', ...
                      'AabsSpeedR','ASpeedR','APSpeedR','ANSpeedR', ...
                      'AabsStimR','AStimR','APStimR','ANStimR'});

sumT1 = join(newtbl3, newtbl4, 'Keys', {'Session','PointTargetCell'});
sumT1(:,{'TargetSpeedR','TargetStimR'})=[];

UniqueG=unique(G);
UniqueS=unique(sessionVals);

tblBridge=table(sessionVals,TargetcellVals,'VariableNames',{'Session','PointTargetCell'});
newtbl5 = groupsummary(newtbl2, {'Session', 'Cell'}, 'mean');

tempVec=cell(length(UniqueS),1);
tempVecAct=cell(length(UniqueS),1);
tempVecActP=cell(length(UniqueS),1);
tempVecActN=cell(length(UniqueS),1);


tempActI=cell(length(UniqueS),1);
for iFOV=1:length(tempActI)
    tempActI{iFOV}={};
end
tempActIPos=tempActI;
tempActINeg=tempActI;
for iG=1:length(UniqueG)
    iFOV=find(UniqueS==sessionVals(iG));
    A1=find(newtbl5.Session==sessionVals(iG));
    B1=find(newtbl5.mean_ActivateI(A1)>CoActITh);
    B2=find(newtbl5.mean_ActivateIPos(A1)>CoActITh);
    B3=find(newtbl5.mean_ActivateINeg(A1)>CoActITh);
    testtbl=newtbl2(G==UniqueG(iG),:);


    tempVec{iFOV}(:,end+1)=newtbl2.Response(G==UniqueG(iG));
    tempActI{iFOV}{end+1} = newtbl2.ActivateI(G==UniqueG(iG));
    tempActIPos{iFOV}{end+1} = newtbl2.ActivateIPos(G==UniqueG(iG));
    tempActINeg{iFOV}{end+1} = newtbl2.ActivateINeg(G==UniqueG(iG));


    tempVecAct{iFOV}(:,end+1)=tempVec{iFOV}(B1,end);
    tempVecActP{iFOV}(:,end+1)=tempVec{iFOV}(B2,end);
    tempVecActN{iFOV}(:,end+1)=tempVec{iFOV}(B3,end);

end

for iFOV=1:length(tempVec)
    rPairVec{iFOV}=zeros(length(tempActI{iFOV}));
    CoPairActIParir{iFOV}=zeros(length(tempActI{iFOV}));

    CoPairActIParirP{iFOV}=zeros(length(tempActI{iFOV}));
    CoPairActIParirN{iFOV}=zeros(length(tempActI{iFOV}));

    CoConPairActIParir{iFOV}=zeros(length(tempActI{iFOV}));


    for i=1:length(tempActI{iFOV})
        for j=i+1:length(tempActI{iFOV})
            sharedI=(tempActI{iFOV}{i}+tempActI{iFOV}{j})>=1;
            rPairVec{iFOV}(i,j)=corr(tempVec{iFOV}(sharedI,i),tempVec{iFOV}(sharedI,j));
            CosharedI1=(tempActI{iFOV}{i}+tempActI{iFOV}{j})==2;
            CoPairActIParir{iFOV}(i,j)=sum(CosharedI1)/length(CosharedI1);

            CosharedI2=(tempActIPos{iFOV}{i}+tempActIPos{iFOV}{j})==2;
            CoPairActIParirP{iFOV}(i,j)=sum(CosharedI2)/length(CosharedI2);


            CosharedI3=(tempActINeg{iFOV}{i}+tempActINeg{iFOV}{j})==2;
            CoPairActIParirN{iFOV}(i,j)=sum(CosharedI3)/length(CosharedI3);

            CoConPairActIParir{iFOV}(i,j)=(sum(CosharedI2)+sum(CosharedI3))/sum(CosharedI1);


        end
    end
    rPairVec{iFOV}=rPairVec{iFOV}+rPairVec{iFOV}';
    CoPairActIParir{iFOV}=CoPairActIParir{iFOV}+CoPairActIParir{iFOV}';
    CoPairActIParirP{iFOV}=CoPairActIParirP{iFOV}+CoPairActIParirP{iFOV}';
    CoPairActIParirN{iFOV}=CoPairActIParirN{iFOV}+CoPairActIParirN{iFOV}';
    CoConPairActIParir{iFOV}=CoConPairActIParir{iFOV}+CoConPairActIParir{iFOV}';


end



for iFOV=1:length(tempVec)

    VecGroupFOV{iFOV}=newtbl3.PointTargetCellGroup(newtbl3.Session==UniqueS(iFOV));


    rVec{iFOV}=corr(tempVec{iFOV},'type','spearman');
    rVecAct{iFOV}=corr(tempVecAct{iFOV},'type','spearman');
    rVecActP{iFOV}=corr(tempVecActP{iFOV});
    if size(tempVecActN{iFOV},1)>1
       rVecActN{iFOV}=corr(tempVecActN{iFOV});
    else
       rVecActN{iFOV}=zeros(size(tempVecActN{iFOV},2))+NaN;
    end


end

CoTargetsLabel={'CoTargets','CoTargetsPos','CoTargetsNeg','CoConTargets'};
CoTargetsVar={CoPairActIParir CoPairActIParirP CoPairActIParirN CoConPairActIParir};

close all

for iCoType=1:length(CoTargetsVar)
figure;
for iFOV=1:length(tempVec)
    subplot(4,4,iFOV)
    AdjComImagesc(CoTargetsVar{iCoType}{iFOV},VecGroupFOV{iFOV},ProcessPar.GroupColor);
    colormap(ResponseMap);
    if iFOV==4
        set(gca,'clim',[-1 1]);
    else
        set(gca,'clim',[-0.5 0.5]);
    end
end
papersizePX=[0 0 16 16];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
print(gcf, [SaveFunConPth  NDataName{iData} CoTargetsLabel{iCoType} '.tif'], '-dtiffn', '-painters');
end

close all
figure;
for iFOV=1:length(tempVec)
    subplot(4,4,iFOV)
    AdjComImagesc(rVec{iFOV},VecGroupFOV{iFOV},ProcessPar.GroupColor);
    colormap(ResponseMap);
    set(gca,'clim',[-1 1]);
end
papersizePX=[0 0 16 16];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
print(gcf, [SaveFunConPth  NDataName{iData} 'VectorStability.tif'], '-dtiffn', '-painters');




figure;
for iFOV=1:length(tempVec)
    subplot(4,4,iFOV)
    AdjComImagesc(rVecAct{iFOV},VecGroupFOV{iFOV},ProcessPar.GroupColor);
    colormap(ResponseMap);
    set(gca,'clim',[-1 1]);
end
papersizePX=[0 0 16 16];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
print(gcf, [SaveFunConPth  NDataName{iData} 'ActVectorStability.tif'], '-dtiffn', '-painters');

close all
figure;
for iFOV=1:length(tempVec)
    subplot(4,4,iFOV)
    AdjComImagesc(rPairVec{iFOV},VecGroupFOV{iFOV},ProcessPar.GroupColor);
    colormap(ResponseMap);
    set(gca,'clim',[-1 1]);
end
papersizePX=[0 0 16 16];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
print(gcf, [SaveFunConPth  NDataName{iData} 'PairwiseSharingActVectorStability.tif'], '-dtiffn', '-painters');



figure;
for iFOV=1:length(tempVec)
    subplot(4,4,iFOV)
    AdjComImagesc(rVecActP{iFOV},VecGroupFOV{iFOV},ProcessPar.GroupColor);
    colormap(ResponseMap);
    set(gca,'clim',[-0.8 0.8]);
end
papersizePX=[0 0 16 16];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
print(gcf, [SaveFunConPth  NDataName{iData} 'PosActVectorStability.tif'], '-dtiffn', '-painters');

figure;
for iFOV=1:length(tempVec)
    subplot(4,4,iFOV)
    AdjComImagesc(rVecActN{iFOV},VecGroupFOV{iFOV},ProcessPar.GroupColor);
    colormap(ResponseMap);
    set(gca,'clim',[-0.8 0.8]);
end
papersizePX=[0 0 16 16];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
print(gcf, [SaveFunConPth  NDataName{iData} 'NegActVectorStability.tif'], '-dtiffn', '-painters');


% intraRgroup{1}=[];
% intraRgroup{2}=[];
% intraRgroup{3}=[];
rMerge={rVec,rVecAct,rPairVec,rVecActP,rVecActN}
rTypeLabel={'Vec','ActVec','PairActVec','ActVecP','ActVecN'};

for irType=1:length(rMerge)
    for igroup=1:3
        intraRgroup{igroup,irType}=[];
    end
end


for iFOV=1:length(tempVec)
    cgroup=newtbl3.PointTargetCellGroup(newtbl3.Session==UniqueS(iFOV));
    for igroup=1:3
        for irType=1:length(rMerge)
            I1=find(cgroup==igroup);
            if length(I1)>2
               r1=rMerge{irType}{iFOV}(I1,I1);
               intraRgroup{igroup,irType}(end+1)=nanmean(r1(triu(true(size(r1)),1)));
            end
        end
    end
end

figure;
for irType=1:length(rMerge)
    subplot(1,length(rMerge),irType)
    stats=ErrorBoxPlotLU(1:3,intraRgroup(:,irType),ProcessPar.GroupColor,[]);
    set(gca,'ylim',[-0.2 0.2])
    yticks([-0.2:0.1:0.2])
    yticklabels([-0.2:0.1:0.2])
    if irType==1
       ylabel('Intra Group Vector Similarity')
    end
    xlabel(rTypeLabel{irType});
end
papersizePX=[0 0 length(rMerge)*6 8];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
print(gcf, [SaveFunConPth  NDataName{iData} 'IntraGroupVectorSimilarity.tif'], '-dtiffn', '-painters');



% G = groupsummary(sumT1,"PointTargetCellGroup",@(x,y) corr(x,y),{["PointSpeedR","PointSpeedR","PointSpeedR","PointStimR","PointStimR","PointStimR"],...
%                                                                 ["ActivateI","ActivateIPos","ActivateINeg","ActivateI","ActivateIPos","ActivateINeg"]})

Predictor={'ActivateI','ActivateIPos','ActivateINeg','AabsSpeedR','ASpeedR','APSpeedR','ANSpeedR','AabsStimR','AStimR','APStimR','ANStimR'};



Factor={'PointSpeedR','SpeedR','PointStimR','StimR','VecR','ActVecR'};
% FactorCov{1}=[newtbl3.StimR newtbl3.SpeedR newtbl3.PointStimR];
% FactorCov{2}=[newtbl3.StimR newtbl3.PointSpeedR newtbl3.PointStimR];
% FactorCov{3}=[newtbl3.StimR newtbl3.SpeedR newtbl3.PointSpeedR];
% FactorCov{4}=[newtbl3.PointStimR newtbl3.PointSpeedR newtbl3.SpeedR];

DisX=[-0.6:0.025:0.6];
% % figure;
% % histPlotLU(sumT1.VecR,DisX,[0.5 0.5 0.5],0.5);
% % hold on;
% % histPlotLU(sumT1.ActVecR,DisX,[0 0 1],0.5);
% % papersizePX=[0 0 12 12];
% % set(gcf, 'PaperUnits', 'centimeters');
% % set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
% % print(gcf, [SaveFunConPth  NDataName{iData} 'TrialStability.svg'], '-dsvg', '-painters');
% % print(gcf, [SaveFunConPth  NDataName{iData} 'TrialStability.tif'], '-dtiffn', '-painters');
% % 
% % 

   GroupPair.CorrName='fdr';
   GroupPair.Q=0.1;
   GroupPair.Pair=[];
   GroupPair.SignY=0.5;
   GroupPair.Plot=1;
   GroupPair.Std=0;      %%%%%%%%%using standard deviation as errorbar
   GroupPair.SamplePlot=1; %%%%%%%%%Plot Individual Sample Point
   GroupPair.SamplePairedPlot=0; %%%%%%%%%Dash line for paired comparison sample
   GroupPair.LimY=[-0.003 GroupPair.SignY*1.2];
   GroupPair.Marker={'o'};   % Datatype=varargin{1};  %%%%%%Datatype=0, non-paired data; 1 for paired data
   GroupPair.GroupName={'VecR','ActVecR'}

figure;
GroupPairViolin=GroupPair;
GroupPairViolin.ViolinLR=[0 0];
GroupPairViolin.SignY=0.8;
GroupPairViolin.LimY=[0 1];
GroupPairViolin.Test='Ranktest';

stats=ErrorViolinHalf(1:2,{sumT1.VecR sumT1.ActVecR},[0.5 0.5 0.5;0 0 1],1,[SaveFunConPth  NDataName{iData} 'TrialStabilityWholeVsAll'])
xticks(1:2)
xticklabels(GroupPairViolin.GroupName)
ylabel('Trial Stability')
yticks(-0.2:0.2:1);
papersizePX=[0 0 8 12];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
print(gcf, [SaveFunConPth  NDataName{iData} 'TrialStability.svg'], '-dsvg', '-painters');
print(gcf, [SaveFunConPth  NDataName{iData} 'TrialStability.tif'], '-dtiffn', '-painters');



FactorVariGroupComp={'VecR','ActVecR'};
figure;
for jFF=1:length(FactorVariGroupComp)
    subplot(1,length(FactorVariGroupComp),jFF)
    clear tempdata
    for igroup=1:3
        tempdata{igroup}=table2array(sumT1(sumT1.PointTargetCellGroup==igroup,FactorVariGroupComp{jFF}));
    %     histPlotLU(sumT1.VecR(sumT1.PointTargetCellGroup==igroup),DisX,ProcessPar.GroupColor(igroup,:),0.5);
    % hold on;
    end
   GroupPair.CorrName='fdr';
   GroupPair.Q=0.1;
   GroupPair.Pair=[];
   GroupPair.SignY=0.2;
   GroupPair.Plot=1;
   GroupPair.Std=0;      %%%%%%%%%using standard deviation as errorbar
   GroupPair.SamplePlot=1; %%%%%%%%%Plot Individual Sample Point
   GroupPair.SamplePairedPlot=0; %%%%%%%%%Dash line for paired comparison sample
   GroupPair.LimY=[-0.003 GroupPair.SignY*1.2];
   GroupPair.Marker={'o'};   % Datatype=varargin{1};  %%%%%%Datatype=0, non-paired data; 1 for paired data
   GroupPair.GroupName={'L','S','N'}
   GroupPairViolin=GroupPair;
   GroupPairViolin.Test='Ranktest';
   GroupPairViolin.ViolinLR=[0 0 0]
   % subplot(1,2,2)
   % stats=ErrorBarPlotLU(1:3,tempdata,[],ProcessPar.GroupColor,2,0,'VecTrialStabilityGroup.txt',GroupPair,[1 2 3]);
   stats=ErrorViolinHalf(1:3,tempdata,ProcessPar.GroupColor,0,[SaveFunConPth  NDataName{iData} FactorVariGroupComp{jFF} 'Group.txt'],GroupPairViolin,[1 2 3])
   xticks(1:3)
   xticklabels(GroupPairViolin.GroupName)
   ylabel(FactorVariGroupComp{jFF})
   set(gca,'ylim',[-0.2 1]);

end

papersizePX=[0 0 16 12];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
print(gcf, [SaveFunConPth  NDataName{iData} 'TrialStabilityGroup.svg'], '-dsvg', '-painters');
print(gcf, [SaveFunConPth  NDataName{iData} 'TrialStabilityGroup.tif'], '-dtiffn', '-painters');



% G = groupsummary(sumT1,"PointTargetCellGroup",@(x,y) corr(x,y,'rows','complete'),{["VecR","VecR","VecR","ActVecR","ActVecR","ActVecR"],...
%                                                                 ["ActivateI","ActivateIPos","ActivateINeg","ActivateI","ActivateIPos","ActivateINeg"]})
% SpaFactor={'Dist','DistXY'}
% SpaPredictor={'Response','ActSpeedR','ActPosSpeedR','ActNegSpeedR','ActStimR','ActPosStimR','ActNegStimR'}

% for iFF=1:length(SpaFactor)
%     figure;
%     for iPP=1:length(SpaPredictor)
%         ParRegress.yLabel=SpaPredictor{iPP};
%         ParRegress.xLabel=SpaFactor{iFF};
%         ParRegress.xLim=[];
%         ParRegress.yLim=[];
% 
%         data1=table2array(newtbl2(:,SpaFactor{iFF}));
%         data2=table2array(newtbl2(:,SpaPredictor{iPP}));
% 
%         dataCov=table2array(newtbl2(:,setdiff(SpaFactor,SpaFactor{iFF})));
% 
%         ValidI=(~isnan(data1))|(~isnan(data2));
%         data1=data1(ValidI);
%         data2=data2(ValidI);
%         dataCov=dataCov(ValidI,:);
% 
%         h(1)=subplotLU(2,length(Predictor),1,iPP,SubplotParam);
%         [OutPut,r,p]=LuPairRegressPlot_Group(data1,data2,newtbl2.PointTargetCellGroup(ValidI),ParRegress)
%         xlabel('');
%         % a=title(Predictor{iPP})
% 
%         h(1)=subplotLU(2,length(Predictor),2,iPP,SubplotParam);
%         [B,BINT,R,RINT,STATS]=regress(data2,[ones(length(dataCov),1) dataCov]);
%         ParRegress.yLim=[];
%         [OutPut,r,p]=LuPairRegressPlot_Group(data1,R,newtbl2.PointTargetCellGroup(ValidI),ParRegress);
% 
% 
%     end
% 
%     papersizePX=[0 0 length(Predictor)*4.5 10];
%     set(gcf, 'PaperUnits', 'centimeters');
%     set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
%     % print(gcf, [SaveFunConPth  NDataName{iData} Factor{iFF} 'MixedGroup.svg'], '-dsvg', '-painters');
%     % print(gcf, [SaveFunConPth  NDataName{iData} Factor{iFF} 'MixedGroup.tif'], '-dtiffn', '-painters');
% 
% 
% end


SpaFactor={'Dist','DistXY'}
SpaPredictor={'Response'};
DistTh=15;

% for iFF=1:length(SpaFactor)
% 
% newtbl2temp=newtbl2(newtbl2.ActivateI==1,:);
% newtbl2tempPos=newtbl2(newtbl2.ActivateIPos==1,:);
% newtbl2tempNeg=newtbl2(newtbl2.ActivateINeg==1,:);
% 
% % newtbl2temp=newtbl2;
Edges=[0 2 10:20:200 250:50:500 900];
% % [Y,E]=discretize(table2array(newtbl2temp(:,SpaFactor{iFF})),100);
% [Y,E]=discretize(table2array(newtbl2temp(:,SpaFactor{iFF})),Edges,'IncludedEdge','right');
% newtbl2temp.Cat=Y;
% newtbl2tempspa = groupsummary(newtbl2temp, {'Cat'}, GroupMethod);
% newtbl2tempspa=groupsummaryBack2OldNames(newtbl2temp,newtbl2tempspa,GroupMethod);
% 
% [Y,E]=discretize(table2array(newtbl2tempPos(:,SpaFactor{iFF})),Edges,'IncludedEdge','right');
% newtbl2tempPos.Cat=Y;
% newtbl2tempspaPos = groupsummary(newtbl2tempPos, {'Cat'}, GroupMethod);
% newtbl2tempspaPos=groupsummaryBack2OldNames(newtbl2tempPos,newtbl2tempspaPos,GroupMethod);
% 
% 
% [Y,E]=discretize(table2array(newtbl2tempNeg(:,SpaFactor{iFF})),Edges,'IncludedEdge','right');
% newtbl2tempNeg.Cat=Y;
% newtbl2tempspaNeg = groupsummary(newtbl2tempNeg, {'Cat'}, GroupMethod);
% newtbl2tempspaNeg=groupsummaryBack2OldNames(newtbl2tempNeg,newtbl2tempspaNeg,GroupMethod);
% 
% 
% figure;
% subplot(1,2,1)
% plot(table2array(newtbl2tempspa(:,SpaFactor{iFF})),newtbl2tempspa.Response,'k')
% hold on;
% plot(table2array(newtbl2tempspaPos(:,SpaFactor{iFF})),newtbl2tempspaPos.Response,'r')
% plot(table2array(newtbl2tempspaNeg(:,SpaFactor{iFF})),newtbl2tempspaNeg.Response,'b')
% xlabel(SpaFactor{iFF})
% ylabel('Response')
% 
% 
% subplot(1,2,2)
% t1=find(table2array(newtbl2tempPos(:,SpaFactor{iFF}))>DistTh);
% t2=find(table2array(newtbl2tempNeg(:,SpaFactor{iFF}))>DistTh);
% % stats=ErrorBoxPlotLU(1:2,{table2array(newtbl2tempPos(t1,SpaFactor{iFF})),table2array(newtbl2tempNeg(t2,SpaFactor{iFF}))},[1 0 0;0 0 1],...
% %     [SaveFunConPth  NDataName{iData} SpaFactor{iFF} "ExVsIn"],0,GroupPair,[1 2]);
% 
% stats=ErrorViolinHalf(1:2,{table2array(newtbl2tempPos(t1,SpaFactor{iFF})),table2array(newtbl2tempNeg(t2,SpaFactor{iFF}))},[1 0 0;0 0 1],0,[SaveFunConPth  NDataName{iData} SpaFactor{iFF} 'ExVsIn.txt'])
% ylabel(SpaFactor{iFF})
% xticks(1:2)
% xticklabels({'Excitation','Inhibition'})
% 
%     papersizePX=[0 0 16 12];
%     set(gcf, 'PaperUnits', 'centimeters');
%     set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
%     print(gcf, [SaveFunConPth  NDataName{iData} SpaFactor{iFF} 'FromPoint.svg'], '-dsvg', '-painters');
%     print(gcf, [SaveFunConPth  NDataName{iData} SpaFactor{iFF} 'FromPoint.tif'], '-dtiffn', '-painters');
% 
% 
% end



groupSepFac={'ActivateI','ActivateIPos','ActivateINeg'};
groupSepLabel={'Response','Excitation','Inhibition'};
groupSepColor=[0 0 0;1 0 0;0 0 1];

close all

for iFF=1:length(SpaFactor)
figure;

clear tempdata
tempdata={};
subplot(1,2,1)

for jFac=1:length(groupSepFac)
    hold on;
    newtbl2tempgroup=newtbl2(table2array(newtbl2(:,groupSepFac{jFac}))==1,:);
    newtbl2tempgroup=newtbl2tempgroup(table2array(newtbl2tempgroup(:,SpaFactor{iFF}))>DistTh,:);

    if jFac>1
       tempdata{end+1}=table2array(newtbl2tempgroup(:,SpaFactor{iFF}));
    end
    [Y,E]=discretize(table2array(newtbl2tempgroup(:,SpaFactor{iFF})),Edges,'IncludedEdge','right');
    newtbl2tempgroup.Cat=Y;
    newtbl2tempgroupspa = groupsummary(newtbl2tempgroup, {'Cat'}, GroupMethod);
    newtbl2tempgroupspa=groupsummaryBack2OldNames(newtbl2tempgroup,newtbl2tempgroupspa,GroupMethod);
    plot(table2array(newtbl2tempgroupspa(:,SpaFactor{iFF})),newtbl2tempgroupspa.Response,'color',groupSepColor(jFac,:));
    hold on;
    xlabel(SpaFactor{iFF})
    ylabel('Neuro Response')

end
subplot(1,2,2)
% stats=ErrorBoxPlotLU(1:3,tempdata,ProcessPar.GroupColor,[],GroupPair,[1 2 3]);
GroupPairViolin=GroupPair;
GroupPairViolin.ViolinLR=[0 0];
GroupPairViolin.SignY=900;
GroupPairViolin.LimY=[0 900];
GroupPairViolin.Test='Ranktest';
stats=ErrorViolinHalf(1:2,tempdata,groupSepColor(2:end,:),0,[SaveFunConPth SpaFactor{iFF} 'FromPoint.txt'],GroupPairViolin,[1 2])
ylabel(SpaFactor{iFF});
yticks([0:300:900])
xticks([1 2])
xticklabels(groupSepLabel(2:end))
    papersizePX=[0 0 16 12];
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
    print(gcf, [SaveFunConPth  NDataName{iData} SpaFactor{iFF} 'FromPoint.svg'], '-dsvg', '-painters');
    print(gcf, [SaveFunConPth  NDataName{iData} SpaFactor{iFF} 'FromPoint.tif'], '-dtiffn', '-painters');




figure;
for jFac=1:length(groupSepFac)
subplot(length(groupSepFac),2,1+(jFac-1)*2)
hold on;
clear tempdata
for igroup=1:3
    newtbl2tempgroup=newtbl2(table2array(newtbl2(:,groupSepFac{jFac}))==1&newtbl2.PointTargetCellGroup==igroup,:);
    newtbl2tempgroup=newtbl2tempgroup(table2array(newtbl2tempgroup(:,SpaFactor{iFF}))>DistTh,:);
    tempdata{igroup}=table2array(newtbl2tempgroup(:,SpaFactor{iFF}));
    [Y,E]=discretize(table2array(newtbl2tempgroup(:,SpaFactor{iFF})),Edges,'IncludedEdge','right');
    newtbl2tempgroup.Cat=Y;
    newtbl2tempgroupspa = groupsummary(newtbl2tempgroup, {'Cat'}, GroupMethod);
    newtbl2tempgroupspa=groupsummaryBack2OldNames(newtbl2tempgroup,newtbl2tempgroupspa,GroupMethod);
    plot(table2array(newtbl2tempgroupspa(:,SpaFactor{iFF})),newtbl2tempgroupspa.Response,'color',ProcessPar.GroupColor(igroup,:));
    hold on;
end
xlabel(SpaFactor{iFF})
ylabel(groupSepLabel{jFac})


subplot(length(groupSepFac),2,2+(jFac-1)*2)
% stats=ErrorBoxPlotLU(1:3,tempdata,ProcessPar.GroupColor,[],GroupPair,[1 2 3]);
GroupPairViolin=GroupPair;
GroupPairViolin.ViolinLR=[0 0 0];
GroupPairViolin.SignY=900;
GroupPairViolin.LimY=[0 900];
GroupPairViolin.Test='Ranktest';
stats=ErrorViolinHalf(1:3,tempdata,ProcessPar.GroupColor,0,[SaveFunConPth SpaFactor{iFF} groupSepFac{jFac} '.txt'],GroupPairViolin,[1 2 3])
ylabel(SpaFactor{iFF});
yticks([0:300:900])
xticks(1:3)
xticklabels(GroupPair.GroupName)

end

    papersizePX=[0 0 16 12];
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
    print(gcf, [SaveFunConPth  NDataName{iData} SpaFactor{iFF} 'FromPointSepGroup.svg'], '-dsvg', '-painters');
    print(gcf, [SaveFunConPth  NDataName{iData} SpaFactor{iFF} 'FromPointSepGroup.tif'], '-dtiffn', '-painters');



end






% figure;
% Ineed=newtbl2temp.ActivateI>0&newtbl2temp.Dist>0
% plot(newtbl2temp.Dist(Ineed),newtbl2temp.Response(Ineed),'.')
% figure;
% Ineed=newtbl2temp.ActivateIPos>0&newtbl2temp.Dist>0
% plot(newtbl2temp.Dist(Ineed),newtbl2temp.Response(Ineed),'.')
% figure;
% Ineed=newtbl2temp.ActivateINeg>0&newtbl2temp.Dist>0
% plot(newtbl2temp.Dist(Ineed),newtbl2temp.Response(Ineed),'.')
% 




% figure;
% plot(Edges,newtbl2tempspa.Response)
% 
% % Edges=[10:50:300]
% newtbl3temp=newtbl2;
% % newtbl3temp.Cat=Y;
% newtbl3temp=newtbl3temp(newtbl3temp.ActivateIPos>0&newtbl3temp.Dist>0,:);
% % figure;hist(newtbl3temp.Response,20)
% % figure;hist(newtbl3temp.Dist,20)
% 



% newtbl3temp=newtbl3temp(newtbl3temp.DistXY>10,:);
% 
% newtbl3tempspa = groupsummary(newtbl3temp, {'Cat'}, GroupMethod);
% newtbl3tempspa=groupsummaryBack2OldNames(newtbl3temp,newtbl3tempspa,GroupMethod);
% 
% 
% figure;
% plot(newtbl3tempspa.Cat,newtbl3tempspa.Response)




% % for iFF=1:length(SpaFactor)
% %     figure;
% %     for iPP=1:length(SpaPredictor)
% %         ParRegress.yLabel=SpaPredictor{iPP};
% %         ParRegress.xLabel=SpaFactor{iFF};
% %         ParRegress.xLim=[];
% %         ParRegress.yLim=[];
% % 
% %         data1=table2array(newtbl2(:,SpaFactor{iFF}));
% %         data2=table2array(newtbl2(:,SpaPredictor{iPP}));
% % 
% %         dataCov=table2array(newtbl2(:,setdiff(SpaFactor,SpaFactor{iFF})));
% % 
% %         ValidI=(~isnan(data1))|(~isnan(data2));
% %         data1=data1(ValidI);
% %         data2=data2(ValidI);
% %         dataCov=dataCov(ValidI,:);
% % 
% %         h(1)=subplotLU(2,length(SpaPredictor),1,iPP,SubplotParam);
% %         [OutPut,r,p]=LuPairRegressPlot_Group(data1,data2,newtbl2.PointTargetCellGroup(ValidI),ParRegress)
% %         xlabel('');
% %         % a=title(Predictor{iPP})
% % 
% %         h(1)=subplotLU(2,length(SpaPredictor),2,iPP,SubplotParam);
% %         [B,BINT,R,RINT,STATS]=regress(data2,[ones(length(dataCov),1) dataCov]);
% %         ParRegress.yLim=[];
% %         [OutPut,r,p]=LuPairRegressPlot_Group(data1,R,newtbl2.PointTargetCellGroup(ValidI),ParRegress);
% % 
% % 
% %     end
% % 
% %     papersizePX=[0 0 length(Predictor)*4.5 10];
% %     set(gcf, 'PaperUnits', 'centimeters');
% %     set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
% %     print(gcf, [SaveFunConPth  NDataName{iData} Factor{iFF} 'MixedGroup.svg'], '-dsvg', '-painters');
% %     print(gcf, [SaveFunConPth  NDataName{iData} Factor{iFF} 'MixedGroup.tif'], '-dtiffn', '-painters');
% % 
% % 
% % end







close all
ParRegress.Color=ProcessPar.GroupColor;
ParRegress.Marker='s';
ParRegress.MarkerSize=10;
ParRegress.Rtype='Spearman';
ParRegress.xLim=[0 0.1];
ParRegress.yLim=[0 0.2];
% ParRegress.xLabel='StimR';
% ParRegress.yLabel=Predictor{iPP};
% [OutPut,r,p]=LuPairRegressPlot_Group(newtbl3.StimR,newtbl3.ActivateI,newtbl3.PointTargetCellGroup,ParRegress)    
ParRegress.xLim=[];
ParRegress.yLim=[];


SubplotParam=[0.05 0.03 0.1 0.2 0.03 0.1]

for iFF=1:length(Factor)
    figure;
    for iPP=1:length(Predictor)
        ParRegress.yLabel=Predictor{iPP};
        ParRegress.xLabel=Factor{iFF};
        ParRegress.xLim=[];
        ParRegress.yLim=[];

        data1=table2array(sumT1(:,Factor{iFF}));
        data2=table2array(sumT1(:,Predictor{iPP}));

        dataCov=table2array(sumT1(:,setdiff(Factor,Factor{iFF})));

        ValidI=(~isnan(data1))|(~isnan(data2));
        data1=data1(ValidI);
        data2=data2(ValidI);
        dataCov=dataCov(ValidI,:);

        h(1)=subplotLU(2,length(Predictor),1,iPP,SubplotParam);
        [OutPut,r,p]=LuPairRegressPlot_Group(data1,data2,sumT1.PointTargetCellGroup(ValidI),ParRegress)
        xlabel('');
        % a=title(Predictor{iPP})

        h(1)=subplotLU(2,length(Predictor),2,iPP,SubplotParam);
        [B,BINT,R,RINT,STATS]=regress(data2,[ones(length(dataCov),1) dataCov]);
        ParRegress.yLim=[];
        [OutPut,r,p]=LuPairRegressPlot_Group(data1,R,sumT1.PointTargetCellGroup(ValidI),ParRegress);

            
    end

    papersizePX=[0 0 length(Predictor)*4.5 10];
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
    print(gcf, [SaveFunConPth  NDataName{iData} Factor{iFF} 'MixedGroup.svg'], '-dsvg', '-painters');
    print(gcf, [SaveFunConPth  NDataName{iData} Factor{iFF} 'MixedGroup.tif'], '-dtiffn', '-painters');


end



SubplotParam=[0.05 0.03 0.025 0.05 0.025 0.025]
ParRegressSub=ParRegress;
for iFF=1:length(Factor)
    figure;
    for iPP=1:length(Predictor)
        ParRegressSub.yLabel=Predictor{iPP};
        ParRegressSub.xLabel=Factor{iFF};

        for igroup=1:3
        ParRegressSub.Color=ParRegress.Color(igroup,:);
        ParRegressSub.xLim=[];
        ParRegressSub.yLim=[];
        subI=find(sumT1.PointTargetCellGroup==igroup);
        data1=table2array(sumT1(subI,Factor{iFF}));
        data2=table2array(sumT1(subI,Predictor{iPP}));

        dataCov=table2array(sumT1(subI,setdiff(Factor,Factor{iFF})));

        ValidI=(~isnan(data1))|(~isnan(data2));
        data1=data1(ValidI);
        data2=data2(ValidI);
        dataCov=dataCov(ValidI,:);


        h(1)=subplotLU(6,length(Predictor),igroup,iPP,SubplotParam);
        [OutPut,r,p]=LuPairRegressPlot_Group(data1,data2,ones(size(data1)),ParRegressSub);

        xlabel('');

        % a=title(Predictor{iPP})

        h(1)=subplotLU(6,length(Predictor),igroup+3,iPP,SubplotParam);
        [B,BINT,R,RINT,STATS]=regress(data2,[ones(length(dataCov),1) dataCov]);
        ParRegressSub.yLim=[min(R) max(R)];
        [OutPut,r,p]=LuPairRegressPlot_Group(data1,R,ones(size(data1)),ParRegressSub);
        if igroup~=3
        xlabel('');
        end
        end
    end

    papersizePX=[0 0 length(Predictor)*4.5 5*6];
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
    print(gcf, [SaveFunConPth  NDataName{iData} Factor{iFF} 'SepGroup.svg'], '-dsvg', '-painters');
    print(gcf, [SaveFunConPth  NDataName{iData} Factor{iFF} 'SepdGroup.tif'], '-dtiffn', '-painters');


end
close all




SubplotParam=[0.1 0.03 0.05 0.2 0.04 0.06]
close all
for iPP=1:length(Predictor)
     ParRegress.xLim=[];
     ParRegress.yLim=[];

    ParRegress.yLabel=Predictor{iPP};
    data1=table2array(sumT1(:,Factor));
    data2=table2array(sumT1(:,Predictor{iPP}));
    
    STATS=regstats(data2,nanzscore(data1),'linear');
    Predictor{iPP} 
    Ifactor=find(STATS.tstat.pval(2:end)'<MultiF_Pth);
    figure;
    papersizePX=[0 0 10 10];
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));

    if length(Ifactor)==1
       ParRegress.xLabel=Factor{Ifactor};
       ParRegress.yLabel=Predictor{iPP};

       [OutPut,r,p]=LuPairRegressPlot_Group(data1(:,Ifactor),data2,sumT1.PointTargetCellGroup,ParRegress);
    elseif length(Ifactor)==2
       ParRegress.xLabel=Factor{Ifactor(1)};
       ParRegress.yLabel=Factor{Ifactor(2)};
       ParRegress.zLabel=Predictor{iPP};
       ParRegress.View=[45 30]; % default 3D view
       [OutPut,r,p]=LuPairRegressPlot3D_Group(data1(:,Ifactor),data2,sumT1.PointTargetCellGroup,ParRegress)
    elseif length(Ifactor)>=3
       for iFF=1:length(Ifactor)
        h(1)=subplotLU(2,length(Ifactor),1,iFF,SubplotParam);
        ParRegress.xLabel=Factor{Ifactor(iFF)};
        ParRegress.yLabel=Predictor{iPP};

        datatemp1=data1(:,Ifactor(iFF));
        dataCov=data1(:,setdiff(Ifactor,Ifactor(iFF)));
        [OutPut,r,p]=LuPairRegressPlot_Group(datatemp1,data2,sumT1.PointTargetCellGroup,ParRegress);
        xlabel('');
        % xticks([])

        if iFF~=1
        ylabel('');

        end
        [B,BINT,R,RINT,STATS]=regress(data2,[ones(length(dataCov),1) dataCov]);
        ParRegress.yLim=[];

        h(1)=subplotLU(2,length(Ifactor),2,iFF,SubplotParam);
        [OutPut,r,p]=LuPairRegressPlot_Group(datatemp1,R,sumT1.PointTargetCellGroup,ParRegress);
        if iFF~=1
        ylabel('');

        end
       end
    papersizePX=[0 0 length(Ifactor)*4.5 10];
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));

    else
        Ifactor=1:length(Factor);
        for iFF=1:length(Ifactor)
        h(1)=subplotLU(2,length(Ifactor),1,iFF,SubplotParam);
        ParRegress.xLabel=Factor{Ifactor(iFF)};
        ParRegress.yLabel=Predictor{iPP};

        datatemp1=data1(:,Ifactor(iFF));
        dataCov=data1(:,setdiff(Ifactor,Ifactor(iFF)));
        [OutPut,r,p]=LuPairRegressPlot_Group(datatemp1,data2,sumT1.PointTargetCellGroup,ParRegress);
        xlabel('');
        % xticks([])

        if iFF~=1
        ylabel('');

        end
        [B,BINT,R,RINT,STATS]=regress(data2,[ones(length(dataCov),1) dataCov]);
        ParRegress.yLim=[];

        h(1)=subplotLU(2,length(Ifactor),2,iFF,SubplotParam);
        [OutPut,r,p]=LuPairRegressPlot_Group(datatemp1,R,sumT1.PointTargetCellGroup,ParRegress);
        if iFF~=1
        ylabel('');

        end
    papersizePX=[0 0 length(Ifactor)*4.5 10];
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));

       end

    end

    print(gcf, [SaveFunConPth  NDataName{iData} 'MultiF' Predictor{iPP} 'MixedGroup.svg'], '-dsvg', '-painters');
    print(gcf, [SaveFunConPth  NDataName{iData} 'MultiF' Predictor{iPP} 'MixedGroup.tif'], '-dtiffn', '-painters');




end

SubplotParam=[0.1 0.03 0.05 0.2 0.04 0.06]
for iPP=1:length(Predictor)
        % ParRegress.xLim=[];
        % ParRegress.yLim=[];
        % 
    % ParRegress.yLabel=Predictor{iPP};
    % data1=table2array(sumT1(:,Factor));
    % data2=table2array(sumT1(:,Predictor{iPP}));
    
    close all
        for igroup=1:3
        ParRegressSub.Color=ParRegress.Color(igroup,:);
        ParRegressSub.xLim=[];
        ParRegressSub.yLim=[];
        subI=find(sumT1.PointTargetCellGroup==igroup);

        data1=table2array(sumT1(subI,Factor));
        data2=table2array(sumT1(subI,Predictor{iPP}));
        dataCov=table2array(sumT1(subI,setdiff(Factor,Factor{iFF})));

        ValidI=(sum(~isnan(data1),2)==size(data1,2))|(~isnan(data2));
        data1=data1(ValidI,:);
        data2=data2(ValidI);
        dataCov=dataCov(ValidI,:);

    STATS=regstats(data2,nanzscore(data1),'linear');
    Predictor{iPP} 
    Ifactor=find(STATS.tstat.pval(2:end)'<MultiF_Pth);
    figure;
    papersizePX=[0 0 10 10];
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));

    if length(Ifactor)==1
       ParRegressSub.xLabel=Factor{Ifactor};
       ParRegressSub.yLabel=Predictor{iPP};
        ParRegressSub.xLim=[];
        ParRegressSub.yLim=[];

       [OutPut,r,p]=LuPairRegressPlot_Group(data1(:,Ifactor),data2,ones(size(data2)),ParRegressSub);
    elseif length(Ifactor)==2
       ParRegressSub.xLabel=Factor{Ifactor(1)};
       ParRegressSub.yLabel=Factor{Ifactor(2)};
       ParRegressSub.zLabel=Predictor{iPP};
       ParRegressSub.View=[45 30]; % default 3D view
       [OutPut,r,p]=LuPairRegressPlot3D_Group(data1(:,Ifactor),data2,ones(size(data2)),ParRegressSub)
    elseif length(Ifactor)>=3
       for iFF=1:length(Ifactor)
        h(1)=subplotLU(2,length(Ifactor),1,iFF,SubplotParam);
        ParRegressSub.xLabel=Factor{Ifactor(iFF)};
        ParRegressSub.yLabel=Predictor{iPP};
        ParRegressSub.xLim=[];
        ParRegressSub.yLim=[];

        datatemp1=data1(:,Ifactor(iFF));
        dataCov=data1(:,setdiff(Ifactor,Ifactor(iFF)));
        [OutPut,r,p]=LuPairRegressPlot_Group(datatemp1,data2,ones(size(data2)),ParRegressSub);
        xlabel('');
        % xticks([])

        if iFF~=1
        ylabel('');

        end
        [B,BINT,R,RINT,STATS]=regress(data2,[ones(length(dataCov),1) dataCov]);
        ParRegressSub.yLim=[];

        h(1)=subplotLU(2,length(Ifactor),2,iFF,SubplotParam);
        [OutPut,r,p]=LuPairRegressPlot_Group(datatemp1,R,ones(size(data2)),ParRegressSub);

        if iFF~=1
        ylabel('');

        end

       end
    papersizePX=[0 0 length(Ifactor)*4.5 10];
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));

    else
        Ifactor=1:length(Factor);
        for iFF=1:length(Ifactor)
        h(1)=subplotLU(2,length(Ifactor),1,iFF,SubplotParam);
        ParRegressSub.xLabel=Factor{Ifactor(iFF)};
        ParRegressSub.yLabel=Predictor{iPP};
        ParRegressSub.xLim=[];
        ParRegressSub.yLim=[];
        
        datatemp1=data1(:,Ifactor(iFF));
        dataCov=data1(:,setdiff(Ifactor,Ifactor(iFF)));
        [OutPut,r,p]=LuPairRegressPlot_Group(datatemp1,data2,ones(size(data2)),ParRegressSub);
        xlabel('');
        % xticks([])
        if iFF~=1
        ylabel('');

        end
        [B,BINT,R,RINT,STATS]=regress(data2,[ones(length(dataCov),1) dataCov]);
        ParRegressSub.yLim=[];

        h(1)=subplotLU(2,length(Ifactor),2,iFF,SubplotParam);
        [OutPut,r,p]=LuPairRegressPlot_Group(datatemp1,R,ones(size(data2)),ParRegressSub);
        if iFF~=1
        ylabel('');

        end
    papersizePX=[0 0 length(Ifactor)*4.5 10];
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
       end

    print(gcf, [SaveFunConPth  NDataName{iData} 'MultiF' Predictor{iPP} 'SepGroup' num2str(igroup) '.svg'], '-dsvg', '-painters');
    print(gcf, [SaveFunConPth  NDataName{iData} 'MultiF' Predictor{iPP} 'SepGroup' num2str(igroup) '.tif'], '-dtiffn', '-painters');

    end


   end


end


close all





CellPopVec=tempVec;
CellPopActVec=tempVecAct;
CellPopActVecP=tempVecActP;
CellPopActVecN=tempVecActN;

CellPopActI=tempActI;
CellPopActPI=tempActIPos;
CellPopActNI=tempActINeg;



save([SaveFunConPth NDataName{iData} 'VecTable.mat'],'tblFOVclean','tblFOVcleanavg','sumT1','SpaFactor','SpaPredictor', 'Predictor','Factor','VecGroupFOV',...
    'CellPopVec','CellPopActVec','CellPopActVecP','CellPopActVecN','CellPopActI','CellPopActPI','CellPopActNI','ProcessPar','PSTHparam','-v7.3')




end

G = groupsummary(sumT1,"PointTargetCellGroup",@(x,y) corr(x,y),{["PointSpeedR","PointSpeedR","PointSpeedR","PointStimR","PointStimR","PointStimR"],...
                                                                ["ActivateI","ActivateIPos","ActivateINeg","ActivateI","ActivateIPos","ActivateINeg"]})

% GroupMethod='corr';
% newtblcorr = groupsummary(newtbl2, {'Session', 'PointTargetCellGroup'}, GroupMethod);



% Define grouping
% [G, sessionVals, cellVals] = findgroups(newtbl2.Session, newtbl2.PointTargetCellGroup);


% G = groupsummary(newtbl2,"PointTargetCellGroup",@(x,y) corr(x,y),{["PointSpeedR","PointSpeedR","PointSpeedR","PointStimR","PointStimR","PointStimR"],...
%                                                                  ["ActivateI","ActivateIPos","ActivateINeg","ActivateI","ActivateIPos","ActivateINeg"]})


[G, sessionVals, TargetcellVals] = findgroups(newtbl2.Session, newtbl2.PointTargetCell);

% Define custom function
myfun = @(x,y) sum(x.*y) / sum(x);

% Apply for each pair
% AI_SpeedR  = splitapply(myfun, newtbl2.ActivateI,    newtbl2.SpeedR, G);
% temp = splitapply(myfun, newtbl2.Response,G);



AIN_SpeedR = splitapply(myfun, newtbl2.ActivateINeg, newtbl2.SpeedR, G);
AI_SpeedR = AIP_SpeedR - AIN_SpeedR;

% AI_StimR   = splitapply(myfun, newtbl2.ActivateI,    newtbl2.StimR,  G);
AIP_StimR  = splitapply(myfun, newtbl2.ActivateIPos, newtbl2.StimR,  G);
AIN_StimR  = splitapply(myfun, newtbl2.ActivateINeg, newtbl2.StimR,  G);
AI_StimR = AIP_StimR - AIN_StimR;
% Build result table using the grouping labels that match
newtbl4 = table( ...
    sessionVals, cellVals, ...
    AI_SpeedR, AIP_SpeedR, AIN_SpeedR, ...
    AI_StimR, AIP_StimR, AIN_StimR, ...
    'VariableNames', {'Session','PointTargetCell', ...
                      'ASpeedR','APSpeedR','ANSpeedR', ...
                      'AStimR','APStimR','ANStimR'});

sumT1 = join(newtbl3, newtbl4, 'Keys', {'Session','PointTargetCell'});
sumT1(:,{'TargetSpeedR','TargetStimR'})=[];


[G, sessionVals, TargetcellVals] = findgroups(newtbl2.Session, newtbl2.PointTargetCell);
% for igroup=1:3
% for iFOV=1:length(tempVec)
%     rVecActGroup{igroup}{iFOV}=corr(tempVecActGroup{igroup}{iFOV});
% end
% end

% for igroup=1:3
% for iFOV=1:length(tempVec)
%     rVecActGroup{igroup}{iFOV}=corr(tempVecActGroup{igroup}{iFOV});
% end
% end

% for igroup=1:3
% for iFOV=1:length(tempVec)
%     rVecActGroup{igroup}{iFOV}=corr(tempVecActGroup{igroup}{iFOV});
% end
% end




%%%%%%Datatype=varargin{1};  %%%%%%default 0, non-paired data; 1 for paired data
%%%%%%SavePath=varargin{2};  %%%%%%Path to save statis to a text file,[] for no saving
%%%%%%GroupPair=varargin{3};  %%%%%%pairs of groups for t-test

stats=ErrorBarPlotLU(1:3,intraRgroup,[],ProcessPar.GroupColor,2,1,[])

        r1Total(end+1)=mean(r1(triu(true(size(r1)),1)));

figure;
for iFOV=1:length(tempVec)
    subplot(4,4,iFOV)

    [~,sortI]=sort(newtbl3.PointSpeedR(newtbl3.Session==UniqueS(iFOV)),'descend');
    imagesc(rVecAct{iFOV}(sortI,sortI));
    colormap(ResponseMap);
    set(gca,'clim',[-0.8 0.8]);
end


figure;
for iFOV=1:length(tempVec)
    subplot(4,4,iFOV)
    [~,sortI]=sort(newtbl3.PointSpeedR(newtbl3.Session==UniqueS(iFOV)),'descend');
    imagesc(rVec{iFOV}(sortI,sortI));
    colormap(ResponseMap);
    set(gca,'clim',[-0.8 0.8]);
end

figure;
for iFOV=1:length(tempVec)
    subplot(4,4,iFOV)
    [~,sortI]=sort(newtbl3.PointSpeedR(newtbl3.Session==UniqueS(iFOV)),'descend');
    imagesc(rVecAct{iFOV}(sortI,sortI));
    colormap(ResponseMap);
    set(gca,'clim',[-0.8 0.8]);
end

figure;
for iFOV=1:length(tempVec)
    subplot(4,4,iFOV)
    [~,sortI]=sort(newtbl3.PointStimR(newtbl3.Session==UniqueS(iFOV)),'descend');
    imagesc(rVecAct{iFOV}(sortI,sortI));
    colormap(ResponseMap);
    set(gca,'clim',[-0.8 0.8]);
end




for iSS=1:length(UniqueS)
    IS1=find(sessionVals==UniqueS(iSS))
    vecTemp=[];
    tempTargetList=TargetcellVals(IS1);
    for iT=1:length(tempTargetList)
        vecTemp(:,end+1)=newtbl2.Response(:,)
        find(sessionVals==UniqueS(iSS))
    end

for i=1:length(UniqueG)
    recVec{i}=newtbl2.Response(G==UniqueG(i));
    recVecAct{i}=newtbl2.Response(G==UniqueG(i)&newtbl2.ActivateI==1);


end
end










for iFOV=1:length(IndexFOVNeed)
    close all
    tbltemp=newtbl2(newtbl2.Session==iFOV,:);
    tblgrouping=tbltemp.PointTargetCell


tempTargetList=Output(iData).TargetCellListFOV{iFOV};
tempPowerI=Output(iData).PowerTargetIFOV{iFOV};
tempPower=unique(tblFOV.UncagingLaserPower);
tempGroup=Output(iData).TargetCellListFunGroupFOV{iFOV};

validI=find(tempPowerI>0);
tempTargetList=tempTargetList(validI);
tempPowerI=tempPowerI(validI);
tempPower=tempPower(tempPowerI);
tempGroup=tempGroup(validI);
for iCell=1:length(validI)

    ActivateI=Output(iData).statCellRes{iFOV}(validI(iCell)).p<activate_Pth(iPth);
    TrialID=unique(tblFOV.TrialID(tblFOV.PointTargetCell==tempTargetList(iCell)&tblFOV.UncagingLaserPower==tempPower(iCell)));
    tempData=[];
    for itrial=1:length(TrialID)
        tempData(:,itrial)=tblFOV.Response(tblFOV.TrialID==TrialID(itrial));
    end

    rgroup(end+1)=tempGroup(iCell);

    r1=corr(tempData);
    r1Total(end+1)=mean(r1(triu(true(size(r1)),1)));
    if sum(ActivateI)>5
       r2=corr(tempData(ActivateI,:));
       r2Total(end+1)=mean(r2(triu(true(size(r2)),1)));
    else
       r2Total(end+1)=NaN;
    end

    temptbl=tblFOV(ismember(tblFOV.TrialID,TrialID),:);
    temptbl.ActivateI=repmat(ActivateI(:),length(TrialID),1);
    temptbl.VecR=repmat(r1Total(end),size(temptbl,1),1);
    temptbl.ActVecR=repmat(r2Total(end),size(temptbl,1),1);
    

    tblFOVclean=[tblFOVclean;temptbl];


end

end













for iFF=1:length(Factor)
    figure;
    for iPP=1:length(Predictor)
        ParRegress.yLabel=Predictor{iPP};
        ParRegress.xLabel=Factor{iFF};
        ParRegress.xLim=[];
        ParRegress.yLim=[];

        data1=table2array(sumT1(:,Factor{iFF}));
        data2=table2array(sumT1(:,Predictor{iPP}));

        dataCov=table2array(sumT1(:,setdiff(Factor,Factor{iFF})));

        ValidI=(~isnan(data1))|(~isnan(data2));
        data1=data1(ValidI);
        data2=data2(ValidI);
        dataCov=dataCov(ValidI,:);

        h(1)=subplotLU(2,length(Predictor),1,iPP,SubplotParam);
        [OutPut,r,p]=LuPairRegressPlot_Group(data1,data2,sumT1.PointTargetCellGroup(ValidI),ParRegress)
        xlabel('');
        % a=title(Predictor{iPP})

        h(1)=subplotLU(2,length(Predictor),2,iPP,SubplotParam);
        [B,BINT,R,RINT,STATS]=regress(data2,[ones(length(dataCov),1) dataCov]);
        ParRegress.yLim=[];
        [OutPut,r,p]=LuPairRegressPlot_Group(data1,R,sumT1.PointTargetCellGroup(ValidI),ParRegress);

            
    end

    papersizePX=[0 0 length(Predictor)*4.5 10];
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
    % print(gcf, [SaveFunConPth  NDataName{iData} Factor{iFF} 'MixedGroup.svg'], '-dsvg', '-painters');
    % print(gcf, [SaveFunConPth  NDataName{iData} Factor{iFF} 'MixedGroup.tif'], '-dtiffn', '-painters');


end





% DisX=0:0.005:0.3;
% figure;
% for igroup=1:3
%     histPlotLUPerc(newtbl2.ActivateI(newtbl2.PointTargetCellGroup==igroup),DisX,ProcessPar.GroupColor(igroup,:),0.2);
%     hold on;
% end
close all
for iPP=1:length(Factor)
ParRegress.Color=ProcessPar.GroupColor;
ParRegress.Marker='s';
ParRegress.MarkerSize=10;
ParRegress.Rtype='pearson';
ParRegress.xLim=[0 0.1];
ParRegress.yLim=[0 0.2];
% ParRegress.xLabel='StimR';
ParRegress.yLabel=Factor{iPP};
% [OutPut,r,p]=LuPairRegressPlot_Group(newtbl3.StimR,newtbl3.ActivateI,newtbl3.PointTargetCellGroup,ParRegress)    

ParRegress.xLim=[-0.1 0.3];
ParRegress.yLim=[0 0.2];
% [OutPut,r,p]=LuPairRegressPlot_Group(newtbl3.PointStimR,newtbl3[:,Factor{iPP}],newtbl3.PointTargetCellGroup,ParRegress)    

figure;
ParRegress.xLim=[-0.1 0.3];
ParRegress.yLim=[0 0.2];
ParRegress.xLabel='TargetStimR';
ParRegress.yLabel=Factor{iPP};
[OutPut,r,p,h]=LuPairRegressPlot_Group_Cov(newtbl3.PointStimR,table2array(newtbl3(:,Factor{iPP})),[newtbl3.StimR newtbl3.SpeedR newtbl3.PointSpeedR],newtbl3.PointTargetCellGroup,ParRegress);   
figure;
ParRegress.yLabel='ActivateStimR';
ParRegress.yLim=[0 0.01];
[OutPut,r,p,h]=LuPairRegressPlot_Group_Cov(newtbl3.PointStimR,table2array(newtbl3(:,Factor{iPP})).*newtbl3.StimR,[newtbl3.SpeedR newtbl3.PointSpeedR],newtbl3.PointTargetCellGroup,ParRegress);    
figure;
ParRegress.yLabel='ActivateSpeedR';
ParRegress.yLim=[0 0.01];
[OutPut,r,p,h]=LuPairRegressPlot_Group_Cov(newtbl3.PointStimR,table2array(newtbl3(:,Factor{iPP})).*newtbl3.SpeedR,[newtbl3.StimR newtbl3.PointSpeedR],newtbl3.PointTargetCellGroup,ParRegress);   






figure;
ParRegress.xLim=[-0.1 0.5];
ParRegress.xLabel='TargetSpeedR';
ParRegress.yLim=[0 0.2];
ParRegress.yLabel=Factor{iPP};
[OutPut,r,p,h]=LuPairRegressPlot_Group_Cov(newtbl3.PointSpeedR,table2array(newtbl3(:,Factor{iPP})),[newtbl3.StimR newtbl3.SpeedR newtbl3.PointStimR],newtbl3.PointTargetCellGroup,ParRegress)    

figure;
ParRegress.xLim=[-0.3 0.5];
ParRegress.yLabel='ActivateSpeedR';
ParRegress.yLim=[0 0.01];
[OutPut,r,p,h]=LuPairRegressPlot_Group_Cov(newtbl3.PointSpeedR,table2array(newtbl3(:,Factor{iPP})).*newtbl3.SpeedR,[newtbl3.StimR newtbl3.PointStimR],newtbl3.PointTargetCellGroup,ParRegress)    

figure;
ParRegress.yLabel='ActivateStimR';
ParRegress.yLim=[0 0.01];
[OutPut,r,p,h]=LuPairRegressPlot_Group_Cov(newtbl3.PointSpeedR,table2array(newtbl3(:,Factor{iPP})).*newtbl3.StimR,[newtbl3.SpeedR newtbl3.PointStimR],newtbl3.PointTargetCellGroup,ParRegress)    




figure;
ParRegress.xLim=[0 0.1];
ParRegress.xLabel='NeuronSpeedR';
[OutPut,r,p,h]=LuPairRegressPlot_Group_Cov(newtbl3.SpeedR,table2array(newtbl3(:,Factor{iPP})),[newtbl3.StimR newtbl3.PointSpeedR newtbl3.PointSpeedR],newtbl3.PointTargetCellGroup,ParRegress)    


figure;
ParRegress.xLim=[0 0.1];
ParRegress.xLabel='NeuronStimR';
[OutPut,r,p,h]=LuPairRegressPlot_Group_Cov(newtbl3.StimR,table2array(newtbl3(:,Factor{iPP})),[newtbl3.PointStimR newtbl3.PointSpeedR newtbl3.SpeedR],newtbl3.PointTargetCellGroup,ParRegress)    

end
% figure;
% ParRegress.xLabel='SpeedR';
% [OutPut,r,p]=LuPairRegressPlot_Group(newtbl3.SpeedR,newtbl3.ActivateI,newtbl3.PointTargetCellGroup,ParRegress)    


figure;
for igroup=1:3
    histPlotLUPerc(newtbl2.ActivateI(newtbl2.PointTargetCellGroup==igroup),DisX,ProcessPar.GroupColor(igroup,:),0.2);
    hold on;
end

unique(tblFOVclean.Point(tblFOVclean.Session==1))

DisX=[-0.2:0.05:0.4];
figure;
histPlotLU(r1Total,DisX,[0.5 0.5 0.5],0.5);
hold on;
histPlotLU(r2Total,DisX,[0 0 1],0.5);

figure;
for igroup=1:3
    histPlotLUPerc(r1Total(rgroup==igroup),DisX,ProcessPar.GroupColor(igroup,:),0.2);
    hold on;
end
figure;
for igroup=1:3
    histPlotLUPerc(r2Total(rgroup==igroup),DisX,ProcessPar.GroupColor(igroup,:),0.2);
    hold on;
end


tbltemp=SLMPointTableTrial(SLMPointTableTrial.Session==iFOV&SLMPointTableTrial.TrialID==TrialID(1),:);
% [W,H] = nnmf(NeuroTrace{iFOV}{iData}(tbltemp.TargetCell==0,SponInd),pcaN);


[W,H] = nnmf(NeuroTrace{iFOV}{iData}(:,SponInd),pcaN,'algorithm','als','replicates',20);

I1=find(tblFOV.TrialID==TrialID(1));
NonTarget=find(tblFOV.TargetCell(I1)==0);
Target=find(tblFOV.TargetCell(I1)==1);

% Wsub=W(NonTarget,:);
% norms = sqrt( sum( Wsub.^2 , 1 ) );  
% WNontarget = bsxfun(@rdivide, Wsub, norms);  
% Asub=NeuroTrace{iFOV}{iData}(NonTarget,SponInd);

WNontarget=W(NonTarget,:);
% end
pinvW = pinv(W);
subNonTargetData=NeuroTrace{iFOV}{iData}(:,SponInd);
% subNonTargetData(FOVres(iFOV).Target,:)=0;
subNonTargetData(Target,:)=0;
HNontarget = (pinvW * subNonTargetData); 
FOVres(iFOV).W=W;
FOVres(iFOV).H=H;
FOVres(iFOV).WNontarget=WNontarget;
FOVres(iFOV).HNontarget=HNontarget;
FOVres(iFOV).Nontarget=NonTarget;
FOVres(iFOV).Target=Target;




SubData=NeuroTrace{iFOV}{iData}(NonTarget,SponInd);
[W_sub,H_sub] = nnmf(SubData,pcaN,'algorithm','als','replicates',20);
FOVres(iFOV).W_sub=W_sub;
FOVres(iFOV).H_sub=H_sub;

% HNontarget=WNontarget'*NeuroTrace{iFOV}{iData}(NonTarget,SponInd);

% for t=1:length(SponInd)
%     HNontarget(:,t) = lsqnonneg(WNontarget, SubData(:,t));



% nComp = size(WNontarget,2);
% nCols = size(Asub,2);
% Hnorm = zeros(nComp,nCols);
% for j = 1:nCols
%   Hnorm(:,j) = lsqnonneg(WNontarget, Asub(:,j) );
% end
% HNontarget = bsxfun(@times, Hnorm, norms(:));
% HNontarget=Hnorm;
% 
% H=H';
% HNontarget=HNontarget';
ProcessPar.GroupColor=[255 51 153;91 20 212;121 247 111]/255;

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
% 
HH=H';

rSpeed=corr(HH(ExcludeStimInd,1:pcaN),nanmean(BehTrace(iFOV).Speed(ExcludeStimInd,:),2),'rows','complete');


% rSpeed=corr(H(:,1:pcaN),nanmean(BehTrace(iFOV).Speed(:,:),2),'rows','complete');

% rSpeedNon=corr(H(ExcludeStimInd,1:pcaN),nanmean(BehTrace(iFOV).Speed(ExcludeStimInd,:),2),'rows','complete');
maxLag=10;
clear rStim c rStim_sub c_sub


for PCI = 1:pcaN
    [c(:,PCI), lags] = xcorr(AmpNormalize(HH(SponInd,PCI),[0 100]), AmpNormalize(nanmean(BehTrace(iFOV).Stim(SponInd,:),2),[0 100]), maxLag, 'coeff');
   
    PostI = find(lags >= 0&lags <= 10);
    PreI=find(lags<0);
    % rStim(PCI, 1)=max(c(PostI,PCI));
    % [~, i1] = max(c(PostI,PCI));
    
    % [~, i1] = max(abs(c(PostI,PCI)-mean(c(PreI,PCI))));

    %%Current used
    [~, i1] = max(abs(c(PostI,PCI))-mean(c(PreI,PCI)));
    rStim(PCI, 1) = c(PostI(i1),PCI)-mean(AmpNormalize(HH(SponInd,PCI)));
    %%Current used


end
[~,pcaIspeed]=max(rSpeed)
[~,pcaIstim]=max(rStim)

FOVres(iFOV).SpeedScore=rSpeed([pcaIspeed pcaIstim]);
FOVres(iFOV).StimScore=rStim([pcaIspeed pcaIstim]);


% for PCI = 1:pcaN
%     figure
%     plot(AmpNormalize(HH_sub(SponInd,PCI)));
%     hold on;
%     plot(nanmean(BehTrace(iFOV).Stim(SponInd,:),2)>1,'color',ProcessPar.GroupColor(2,:))
% end


HH_sub=H_sub';
SpeedSpon=nanmean(BehTrace(iFOV).Speed(ExcludeStimInd,:),2);
rSpeed_sub=corr(HH_sub(ExcludeStimInd,1:pcaN),SpeedSpon,'rows','complete');

for PCI = 1:pcaN
    [c_sub(:,PCI), lags] = xcorr(AmpNormalize(HH_sub(SponInd,PCI),[0 100]), AmpNormalize(nanmean(BehTrace(iFOV).Stim(SponInd,:),2),[0 100]), maxLag, 'coeff');
    PostI = find(lags >= 0&lags <= 10);
    PreI=find(lags<0);
  %%Current used
    [~, i1] = max(abs(c_sub(PostI,PCI))-mean(c_sub(PreI,PCI)));
    rStim_sub(PCI, 1) = c_sub(PostI(i1),PCI)-mean(AmpNormalize(HH_sub(SponInd,PCI)));
    %%Current used


end


PP=[0.1 0.03 0.05 0.1];
      % P.xLeft=PP(1);
      % P.xRight=PP(2);
      % P.yTop=PP(3);
      % P.yBottom=PP(4);



TBLall=[];
TBLallNontarget=[];
TBLall_sub=[];
    
for iFOV=1:length(IndexFOVNeed)
close all
tblFOV=SLMGroupTableTrial(SLMGroupTableTrial.Session==iFOV,:);
TrialID=unique(tblFOV.TrialID);


DataVecInfo=[];
DataVec=[];
ScoreNontarget=[];
Score=[];

Score_sub=[];

% DataVecTBL=[];
clear DataVecTBL
tempScore=[];
for iTrial=1:length(TrialID)
    DataVec(:,iTrial)= tblFOV.Response(tblFOV.TrialID==TrialID(iTrial),:);
    I1=find(tblFOV.TrialID==TrialID(iTrial));

    % I1=find(tblFOV.TrialID==TrialID(1));
    NonTarget=find(tblFOV.TargetCell(I1)==0);
    Target=find(tblFOV.TargetCell(I1)==1);

    temp=tblFOV(tblFOV.TrialID==TrialID(iTrial),:);
    % DataVecInfo(iTrial,:)=[tblFOV.Group(I1(1)) tblFOV.Whisk(I1(1)) tblFOV.PowerZero(I1(1)) tblFOV.Time(I1(1))  tblFOV.Speed(I1(1))];
    % DataVec(tblFOV.TargetCell(I1)==1,iTrial)=0;
    DataVecTBL(iTrial,:)=tblFOV(I1(1),:);
    DataVecTBL.SpeedR(iTrial)=mean(tblFOV.SpeedR(I1(NonTarget),:));
    DataVecTBL.StimR(iTrial)=mean(tblFOV.StimR(I1(NonTarget),:));
    % DataVecTBL(iTrial,:).TargetSpeedR=mean(tblFOV.SpeedR(I1(Target),:));
    % DataVecTBL(iTrial,:).TargetStimR=mean(tblFOV.StimR(I1(Target),:));
    if temp.Group(1)~=3
        tempScore(iTrial,:)=[FOVres(iFOV).SpeedScore(temp.Group(1)) FOVres(iFOV).StimScore(temp.Group(1))];
        % DataVecTBL.TargetStimScore(iTrial)=FOVres.StimScore(temp.Group(1));     
    else
        tempScore(iTrial,:)=[NaN NaN];
    end


    Score(:,iTrial)=lsqnonneg(FOVres(iFOV).W(:,FOVres(iFOV).ComI), DataVec(:,iTrial));
    ScoreNontarget(:,iTrial)=lsqnonneg(FOVres(iFOV).WNontarget(:,FOVres(iFOV).ComI), DataVec(NonTarget,iTrial));
    Score_sub(:,iTrial)=lsqnonneg(FOVres(iFOV).W_sub(:,FOVres(iFOV).ComI_sub), DataVec(NonTarget,iTrial));

end

DataVecTBL=[DataVecTBL array2table(tempScore,"VariableNames",{'TargetSpeedScore','TargetWhiskerScore'})]
Score=Score';
ScoreNontarget=ScoreNontarget';
Score_sub=Score_sub';


DataVecRaw=DataVec;
DataVecNontarget=DataVec(NonTarget,:);

% Score=DataVec'*WSub;
% ScoreNontarget=DataVecNontarget'*WNontargetSub;

FOVres(iFOV).DataVec=DataVec;
FOVres(iFOV).InfoTBL=DataVecTBL;
% FOVres(iFOV).NonTarget=NonTarget;
FOVres(iFOV).DataVecNonTarget=DataVecNontarget;


    pinvW = pinv(FOVres(iFOV).W(:,FOVres(iFOV).ComI));
    Score = (pinvW * FOVres(iFOV).DataVec)'; 
    Score(:,1) = Score(:,1) - [ones(length(Score(:,1)),1) DataVecTBL.Speed]*FOVres(iFOV).SpeedReg.B;


    % pinvW = pinv(FOVres(iFOV).W(:,FOVres(iFOV).ComI));
    subNonTargetData=FOVres(iFOV).DataVec;
    % subNonTargetData(FOVres(iFOV).Target,:)=0;
    subNonTargetData(Target,:)=0;
    % ScoreNontarget = (pinvW * subNonTargetData)'; 
    ScoreNontarget(:,1) = ScoreNontarget(:,1) - [ones(length(ScoreNontarget(:,1)),1) DataVecTBL.Speed]*FOVres(iFOV).SpeedRegNontarget.B;

    pinvW_sub = pinv(FOVres(iFOV).W_sub(:,FOVres(iFOV).ComI_sub));
    % Score_sub = (pinvW_sub * FOVres(iFOV).DataVecNonTarget)'; 
    Score_sub(:,1) = Score_sub(:,1) - [ones(length(Score_sub(:,1)),1) DataVecTBL.Speed]*FOVres(iFOV).SpeedReg_sub.B;

maxScore=max(Score);
% maxScore=max(ScoreNontarget)

% SpeedTh=0.5;
clear ScoreGroup SpeedGroup ;
% maxScore
for iWhisk=1:2
for iGroup=1:3
    I1=find(DataVecTBL.Group==iGroup&DataVecTBL.Whisk==ProcessPar.VolOut(iWhisk)&DataVecTBL.PowerZero==0&abs(DataVecTBL.Speed)<1000);
    ScoreGroup(iGroup,iWhisk).Score=Score(I1,:);
    ScoreGroup(iGroup,iWhisk).ScoreNonTarget=ScoreNontarget(I1,:);
    ScoreGroup(iGroup,iWhisk).Time=DataVecTBL.Time(I1);
    SpeedGroup{iGroup,iWhisk}=DataVecTBL.Speed(I1);
    ScoreGroup(iGroup,iWhisk).Score_sub=Score_sub(I1,:);
end

% figure;hold on;



figure;
for iScore=1:2
TempScore={};

    for iGroup=1:3
        % plot(ScoreGroup(iGroup).Time,ScoreGroup(iGroup).Score(:,1),'Color',ProcessPar.GroupColor(iGroup,:));
        NeuroScore{iGroup,iScore}=ScoreGroup(iGroup,iWhisk).Score(iScore,:);
        NeuroScoreNontarget{iGroup,iScore}=ScoreGroup(iGroup,iWhisk).ScoreNonTarget(iScore,:);
        NeuroScore_sub{iGroup,iScore}=ScoreGroup(iGroup,iWhisk).Score_sub(iScore,:);


    end


   GroupPair.CorrName='fdr';
   GroupPair.Q=0.1;
   GroupPair.Pair=[];
   GroupPair.SignY=maxScore(iScore);
   GroupPair.Plot=1;
   GroupPair.Std=0;      %%%%%%%%%using standard deviation as errorbar
   GroupPair.SamplePlot=1; %%%%%%%%%Plot Individual Sample Point
   GroupPair.SamplePairedPlot=0; %%%%%%%%%Dash line for paired comparison sample
   GroupPair.LimY=[-0.003 GroupPair.SignY*1.6];
   GroupPair.Marker={'o'};   % Datatype=varargin{1};  %%%%%%Datatype=0, non-paired data; 1 for paired data
   % SavePath=varargin{2};  %%%%%%Save the stats to txt file
   % GroupPair=varargin{3};
   % RGroupID=varargin{4};


   GroupPairSpeed=GroupPair;
   GroupPairSpeed.SignY=3;
   GroupPairSpeed.LimY=[-0.2 GroupPairSpeed.SignY*1.2];

% stats=ErrorBarPlotLU(1:3,NeuroScore,[],ProcessPar.GroupColor,2,0,[],GroupPair)
   % SavePath=varargin{4};
   % GroupPair=varargin{5};
   %
    subplotLU(2,4,iScore,1)
    stats=ErrorBoxPlotLU(1:3,NeuroScore(:,iScore),ProcessPar.GroupColor,[],GroupPair,[1 2 3]);
    set(gca,'ylim',GroupPair.LimY);

    if iScore==2
       xlabel('All')
    end

     ylabel(ScoreLabel{iScore})



    subplotLU(2,4,iScore,2)
    stats=ErrorBoxPlotLU(1:3,NeuroScoreNontarget(:,iScore),ProcessPar.GroupColor,[],GroupPair,[1 2 3])
    set(gca,'ylim',GroupPair.LimY);

    if iScore==2
       xlabel('NonTarget')
    end

    subplotLU(2,4,iScore,3)
    stats=ErrorBoxPlotLU(1:3,NeuroScore_sub(:,iScore),ProcessPar.GroupColor,[],GroupPair,[1 2 3])
    set(gca,'ylim',GroupPair.LimY);

    if iScore==2
       xlabel('NonTarget')
    end


    if iScore==2
    subplotLU(2,4,iScore,4)
    stats=ErrorBoxPlotLU(1:3,SpeedGroup(:,iWhisk),ProcessPar.GroupColor,[],GroupPairSpeed,[1 2 3])
    set(gca,'ylim',GroupPairSpeed.LimY);

    % if iScore==2
       xlabel('Speed')
    % end
    end
    % subplotLU(2,2,1,1)
    % stats=ErrorBoxPlotLU(1:3,NeuroScore(:,iScore),ProcessPar.GroupColor,[],GroupPair,[1 2 3]);
    % subplotLU(2,2,1,2)
    % stats=ErrorBoxPlotLU(1:3,NeuroScoreNontarget(:,iScore),ProcessPar.GroupColor,[],GroupPair,[1 2 3])

end


    papersizePX=[0 0 15 12];
    set(gcf, 'PaperUnits', 'centimeters');
    set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
    print(gcf, [SaveFunSub ProcessPar.VolOutLabel{iWhisk} 'FOV' num2str(iFOV) 'NeuroCom' '.svg'], '-dsvg', '-painters');
    print(gcf, [SaveFunSub ProcessPar.VolOutLabel{iWhisk} 'FOV' num2str(iFOV) 'NeuroCom' '.tif'], '-dtiffn', '-painters');



end



    ScoreNontarget=array2table(ScoreNontarget,"VariableNames",{'SpeedScore','WhiskerScore'});
    Score=array2table(Score,"VariableNames",{'SpeedScore','WhiskerScore'});
    Score_sub=array2table(Score_sub,"VariableNames",{'SpeedScore','WhiskerScore'});

    TBLallNontarget=[TBLallNontarget;[ScoreNontarget FOVres(iFOV).InfoTBL(:,2:end)]];
    TBLall=[TBLall;[Score FOVres(iFOV).InfoTBL(:,2:end)]];
    TBLall_sub=[TBLall_sub;[Score_sub FOVres(iFOV).InfoTBL(:,2:end)]];





clear DataVecTBL DataVec DataVecRaw DataVecNontarget

end



save([SaveFunSub NDataName{iData} date 'nnmfResult.mat']);



TBLall=[];
TBLallNontarget=[];
TBLall_sub=[];

for iFOV=1:length(IndexFOVNeed)
    % if iFOV==6
    %    continue;
    % end

    tblFOV=SLMGroupTableTrial(SLMGroupTableTrial.Session==iFOV,:);
    TrialID=unique(tblFOV.TrialID);
    ScoreNontarget=[];
    Score=[];
    Score_sub=[];
    DataVec=[];
    tempScore=[];
    clear DataVecTBL
    for iTrial=1:size(FOVres(iFOV).DataVec,2)
        DataVec(:,iTrial)= tblFOV.Response(tblFOV.TrialID==TrialID(iTrial),:);
        I1=find(tblFOV.TrialID==TrialID(iTrial));
    
        % I1=find(tblFOV.TrialID==TrialID(1));
        NonTarget=find(tblFOV.TargetCell(I1)==0);
        Target=find(tblFOV.TargetCell(I1)==1);
    
        % temp=tblFOV(tblFOV.TrialID==TrialID(iTrial),:);
        % DataVecInfo(iTrial,:)=[tblFOV.Group(I1(1)) tblFOV.Whisk(I1(1)) tblFOV.PowerZero(I1(1)) tblFOV.Time(I1(1))  tblFOV.Speed(I1(1))];
        % DataVec(tblFOV.TargetCell(I1)==1,iTrial)=0;
        DataVecTBL(iTrial,:)=tblFOV(I1(1),:);
        DataVecTBL.SpeedR(iTrial)=mean(tblFOV.SpeedR(I1(NonTarget),:));
        DataVecTBL.StimR(iTrial)=mean(tblFOV.StimR(I1(NonTarget),:));
        % DataVecTBL(iTrial,:).TargetSpeedR=mean(tblFOV.SpeedR(I1(Target),:));
        % DataVecTBL(iTrial,:).TargetStimR=mean(tblFOV.StimR(I1(Target),:));
        if temp.Group(1)~=3
            tempScore(iTrial,:)=[FOVres(iFOV).SpeedScore(temp.Group(1)) FOVres(iFOV).StimScore(temp.Group(1))];
            % DataVecTBL.TargetStimScore(iTrial)=FOVres.StimScore(temp.Group(1));     
        else
            tempScore(iTrial,:)=[NaN NaN];
        end
        % 
        % Keep Score positive
        Score(:,iTrial)=lsqnonneg(FOVres(iFOV).W(:,FOVres(iFOV).ComI), FOVres(iFOV).DataVec(:,iTrial));
        ScoreNontarget(:,iTrial)=lsqnonneg(FOVres(iFOV).WNontarget(:,FOVres(iFOV).ComI), FOVres(iFOV).DataVecNonTarget(:,iTrial));
        Score_sub(:,iTrial)=lsqnonneg(FOVres(iFOV).W_sub(:,FOVres(iFOV).ComI_sub), FOVres(iFOV).DataVecNonTarget(:,iTrial));

        % %%Keep Score positive

    end

    DataVecTBL=[DataVecTBL array2table(tempScore,"VariableNames",{'TargetSpeedScore','TargetWhiskerScore'})];

    Score=Score';
    ScoreNontarget=ScoreNontarget';
    Score_sub=Score_sub';
    % 
    % Score=zscore(Score);
    % ScoreNontarget=zscore(ScoreNontarget);
    % Score_sub=zscore(Score_sub);
    % 
    % pinvW = pinv(FOVres(iFOV).W(:,FOVres(iFOV).ComI));
    % Score = (pinvW * FOVres(iFOV).DataVec)'; 
    Score(:,1) = Score(:,1) - [ones(length(Score(:,1)),1) DataVecTBL.Speed]*FOVres(iFOV).SpeedReg.B;

% x=[min(data1(:)) max(data1(:))]';
% y=[[1 1]' x]*B;
% close all
% plot(lags,c_sub(:,12))
    % subNonTargetData=FOVres(iFOV).DataVec;
    % subNonTargetData(Target,:)=0;
    % ScoreNontarget = (pinvW * subNonTargetData)'; 
    ScoreNontarget(:,1) =  ScoreNontarget(:,1) - [ones(length(Score(:,1)),1) DataVecTBL.Speed]*FOVres(iFOV).SpeedRegNontarget.B;

    % pinvW_sub = pinv(FOVres(iFOV).W_sub(:,FOVres(iFOV).ComI_sub));
    % Score_sub = (pinvW_sub * FOVres(iFOV).DataVecNonTarget)'; 
    % 
    Score_sub(:,1) = Score_sub(:,1) - [ones(length(Score(:,1)),1) DataVecTBL.Speed]*FOVres(iFOV).SpeedReg_sub.B;

    % ScoreNontarget=array2table(ScoreNontarget,"VariableNames",{'SpeedScore','WhiskerScore'});
    % Score=array2table(Score,"VariableNames",{'SpeedScore','WhiskerScore'});

    % ScoreNontarget=zscore(FOVres(iFOV).DataVecNonTarget'*FOVres(iFOV).WNontarget(:,FOVres(iFOV).ComI));
    % ScoreNontarget=array2table(ScoreNontarget',"VariableNames",{'SpeedScore','WhiskerScore'});
    % Score=array2table(Score',"VariableNames",{'SpeedScore','WhiskerScore'});
    % 
    ScoreNontarget=array2table(ScoreNontarget,"VariableNames",{'SpeedScore','WhiskerScore'});
    Score=array2table(Score,"VariableNames",{'SpeedScore','WhiskerScore'});
    Score_sub=array2table(Score_sub,"VariableNames",{'SpeedScore','WhiskerScore'});


    % ScoreNontarget=array2table(zscore(ScoreNontarget),"VariableNames",{'SpeedScore','WhiskerScore'});
    % Score=array2table(zscore(Score),"VariableNames",{'SpeedScore','WhiskerScore'});
    % Score_sub=array2table(zscore(Score_sub),"VariableNames",{'SpeedScore','WhiskerScore'});

    % ScoreNontarget=array2table(ScoreNontarget,"VariableNames",{'SpeedScore','WhiskerScore'});
    % Score=array2table(Score,"VariableNames",{'SpeedScore','WhiskerScore'});

    TBLallNontarget=[TBLallNontarget;[ScoreNontarget FOVres(iFOV).InfoTBL(:,2:end)]];
    TBLall=[TBLall;[Score FOVres(iFOV).InfoTBL(:,2:end)]];
    TBLall_sub=[TBLall_sub;[Score_sub FOVres(iFOV).InfoTBL(:,2:end)]];

    
end

iWin=1;
writetable(TBLallNontarget,[SaveFunConDim 'NonTargetSLMScoreTrialDynWin' num2str(iWin) '.csv']);
writetable(TBLall,[SaveFunConDim 'AllSLMScoreTrialDynWin' num2str(iWin) '.csv']);


for iFOV=1:length(IndexFOVNeed)
    % if iFOV==6
    %    continue;
    % end

    tblFOV=SLMGroupTableTrial(SLMGroupTableTrial.Session==iFOV,:);
    TrialID=unique(tblFOV.TrialID);
    ScoreNontarget=[];
    Score=[];
    Score_sub=[];
    DataVec=[];
    clear DataVecTBL
    for iTrial=1:size(FOVres(iFOV).DataVec,2)
        DataVec(:,iTrial)= tblFOV.Response(tblFOV.TrialID==TrialID(iTrial),:);
        I1=find(tblFOV.TrialID==TrialID(iTrial));
    
        % I1=find(tblFOV.TrialID==TrialID(1));
        NonTarget=find(tblFOV.TargetCell(I1)==0);
        Target=find(tblFOV.TargetCell(I1)==1);
    
        % temp=tblFOV(tblFOV.TrialID==TrialID(iTrial),:);
        % DataVecInfo(iTrial,:)=[tblFOV.Group(I1(1)) tblFOV.Whisk(I1(1)) tblFOV.PowerZero(I1(1)) tblFOV.Time(I1(1))  tblFOV.Speed(I1(1))];
        % DataVec(tblFOV.TargetCell(I1)==1,iTrial)=0;
        DataVecTBL(iTrial,:)=tblFOV(I1(1),:);
        DataVecTBL(iTrial,:).SpeedR=mean(tblFOV.SpeedR(I1(NonTarget),:));
        DataVecTBL(iTrial,:).StimR=mean(tblFOV.StimR(I1(NonTarget),:));
        % DataVecTBL(iTrial,:).TargetSpeedR=mean(tblFOV.SpeedR(I1(Target),:));
        % DataVecTBL(iTrial,:).TargetStimR=mean(tblFOV.StimR(I1(Target),:));
    
        % 
        % %Keep Score positive
        Score(:,iTrial)=lsqnonneg(FOVres(iFOV).W(:,FOVres(iFOV).ComI), FOVres(iFOV).DataVec(:,iTrial));
        ScoreNontarget(:,iTrial)=lsqnonneg(FOVres(iFOV).WNontarget(:,FOVres(iFOV).ComI), FOVres(iFOV).DataVecNonTarget(:,iTrial));
        Score_sub(:,iTrial)=lsqnonneg(FOVres(iFOV).W_sub(:,FOVres(iFOV).ComI_sub), FOVres(iFOV).DataVecNonTarget(:,iTrial));
        % 
        % % %%Keep Score positive

    end
    Score=Score';
    ScoreNontarget=ScoreNontarget';
    Score_sub=Score_sub';

    % Score=zscore(Score);
    % ScoreNontarget=zscore(ScoreNontarget);
    % Score_sub=zscore(Score_sub);

    pinvW = pinv(FOVres(iFOV).W(:,FOVres(iFOV).ComI));
    % Score = (pinvW * FOVres(iFOV).DataVec)'; 
    % Score(:,1) = Score(:,1) - [ones(length(Score(:,1)),1) DataVecTBL.Speed]*FOVres(iFOV).SpeedReg.B;

% x=[min(data1(:)) max(data1(:))]';
% y=[[1 1]' x]*B;
% close all
% plot(lags,c_sub(:,12))
    subNonTargetData=FOVres(iFOV).DataVec;
    subNonTargetData(Target,:)=0;
    % ScoreNontarget = (pinvW * subNonTargetData)'; 
    % ScoreNontarget(:,1) =  ScoreNontarget(:,1) - [ones(length(Score(:,1)),1) DataVecTBL.Speed]*FOVres(iFOV).SpeedRegNontarget.B;

    pinvW_sub = pinv(FOVres(iFOV).W_sub(:,FOVres(iFOV).ComI_sub));
    % Score_sub = (pinvW_sub * FOVres(iFOV).DataVecNonTarget)'; 

    % Score_sub(:,1) = Score_sub(:,1) - [ones(length(Score(:,1)),1) DataVecTBL.Speed]*FOVres(iFOV).SpeedReg_sub.B;

    % ScoreNontarget=array2table(ScoreNontarget,"VariableNames",{'SpeedScore','WhiskerScore'});
    % Score=array2table(Score,"VariableNames",{'SpeedScore','WhiskerScore'});

    % ScoreNontarget=zscore(FOVres(iFOV).DataVecNonTarget'*FOVres(iFOV).WNontarget(:,FOVres(iFOV).ComI));
    % ScoreNontarget=array2table(ScoreNontarget',"VariableNames",{'SpeedScore','WhiskerScore'});
    % Score=array2table(Score',"VariableNames",{'SpeedScore','WhiskerScore'});
    % 
    ScoreNontarget=array2table(ScoreNontarget,"VariableNames",{'SpeedScore','WhiskerScore'});
    Score=array2table(Score,"VariableNames",{'SpeedScore','WhiskerScore'});
    Score_sub=array2table(Score_sub,"VariableNames",{'SpeedScore','WhiskerScore'});


    % ScoreNontarget=array2table(zscore(ScoreNontarget),"VariableNames",{'SpeedScore','WhiskerScore'});
    % Score=array2table(zscore(Score),"VariableNames",{'SpeedScore','WhiskerScore'});
    % Score_sub=array2table(zscore(Score_sub),"VariableNames",{'SpeedScore','WhiskerScore'});

    % ScoreNontarget=array2table(ScoreNontarget,"VariableNames",{'SpeedScore','WhiskerScore'});
    % Score=array2table(Score,"VariableNames",{'SpeedScore','WhiskerScore'});

    TBLallNontarget=[TBLallNontarget;[ScoreNontarget FOVres(iFOV).InfoTBL(:,2:end)]];
    TBLall=[TBLall;[Score FOVres(iFOV).InfoTBL(:,2:end)]];
    TBLall_sub=[TBLall_sub;[Score_sub FOVres(iFOV).InfoTBL(:,2:end)]];

    
end



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
SpeedTh=1000
TBLallavg = groupsummary(TBLall(TBLall.Whisk==0&TBLall.PowerZero==0&abs(TBLall.Speed)<SpeedTh,:), {'Session', 'Group'}, GroupMethod);
TBL_subavg= groupsummary(TBLall_sub(TBLall_sub.Whisk==0&TBLall_sub.PowerZero==0&abs(TBLall_sub.Speed)<SpeedTh,:),{'Session', 'Group'}, GroupMethod);
% 
% TBLallavg = groupsummary(TBLall(TBLall.Whisk==0&TBLall.PowerZero==0&TBLall.Group~=3&abs(TBLall.Speed)<SpeedTh,:), {'Session', 'Group'}, GroupMethod);
% TBL_subavg= groupsummary(TBLall_sub(TBLall.Whisk==0&TBLall.PowerZero==0&TBLall.Group~=3&abs(TBLall.Speed)<SpeedTh,:),{'Session', 'Group'}, GroupMethod);
% 
% SpeedTh=1000;
% TBLallavg = groupsummary(TBLall(TBLall.Whisk==0&TBLall.PowerZero==0&TBLall.Group~=3&TBLall.Speed<SpeedTh,:), {'Session', 'Group'}, GroupMethod);
% TBL_subavg= groupsummary(TBLall_sub(TBLall.Whisk==0&TBLall.PowerZero==0&TBLall.Group~=3&TBLall.Speed<SpeedTh,:),{'Session', 'Group'}, GroupMethod);
% 
oldNames = TBLallavg.Properties.VariableNames;
% Find which names start with 'mean_'
meanVars = startsWith(oldNames, [GroupMethod '_']);
% Remove 'mean_' prefix
newNames = oldNames;
newNames(meanVars) = erase(oldNames(meanVars), [GroupMethod '_']);
TBLallavg.Properties.VariableNames=newNames;
TBL_subavg.Properties.VariableNames=newNames;



% for igroup=1:3
%     tempData{igroup,1}=TBLallavg.mean_SpeedScore(TBLallavg.Group==igroup);
% end
% a=ErrorBarPlotLU(1:3,tempData,[],ProcessPar.GroupColor,2,1,[],GroupPair,[1 1 1])
% 

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






SpeedTh=1;
% A=TBLallNontarget(TBLallNontarget.Whisk==0&TBLallNontarget.PowerZero==0&TBLallNontarget.Speed<SpeedTh,:);
% B=TBLall(TBLall.Whisk==0&TBLall.PowerZero==0&TBLall.Speed<SpeedTh,:);
figure;
clear dataAll dataNonT dataAll_sub statsSpeed
for iWhisk=1:length(ProcessPar.VolOut)

tempNonT=TBLallNontarget(TBLallNontarget.Whisk==ProcessPar.VolOut(iWhisk)&TBLallNontarget.PowerZero==0&TBLallNontarget.Speed<SpeedTh,:);
tempAll=TBLall(TBLall.Whisk==ProcessPar.VolOut(iWhisk)&TBLall.PowerZero==0&TBLall.Speed<SpeedTh,:);
tempAll_sub=TBLall_sub(TBLall_sub.Whisk==ProcessPar.VolOut(iWhisk)&TBLall_sub.PowerZero==0&TBLall_sub.Speed<SpeedTh,:);
% 
% tempNonT=TBLallNontarget(TBLallNontarget.Whisk==ProcessPar.VolOut(iWhisk)&TBLallNontarget.PowerZero==0&abs(TBLallNontarget.Speed)<SpeedTh,:);
% tempAll=TBLall(TBLall.Whisk==ProcessPar.VolOut(iWhisk)&TBLall.PowerZero==0&abs(TBLall.Speed)<SpeedTh,:);
% tempAll_sub=TBLall_sub(TBLall_sub.Whisk==ProcessPar.VolOut(iWhisk)&TBLall_sub.PowerZero==0&abs(TBLall_sub.Speed<SpeedTh),:);


for iGroup=1:3
    dataAll{iGroup,iWhisk}=tempAll.SpeedScore(tempAll.Group==iGroup);
    dataNonT{iGroup,iWhisk}=tempNonT.SpeedScore(tempNonT.Group==iGroup);
    dataAll_sub{iGroup,iWhisk}=tempAll_sub.SpeedScore(tempAll_sub.Group==iGroup);

end

subplotLU(2,3,iWhisk,1)
% statsSpeed(iWhisk,1)=ErrorBoxPlotLU(1:3,dataAll(:,iWhisk),ProcessPar.GroupColor,[],GroupPair)
statsSpeed(iWhisk,1)=ErrorBarPlotLU(1:3,dataAll(:,iWhisk),[],ProcessPar.GroupColor,2,0,[],GroupPair)

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
clear dataAll dataNonT dataAll_sub
for iWhisk=1:length(ProcessPar.VolOut)

tempNonT=TBLallNontarget(TBLallNontarget.Whisk==ProcessPar.VolOut(iWhisk)&TBLallNontarget.PowerZero==0&TBLallNontarget.Speed<SpeedTh,:);
tempAll=TBLall(TBLall.Whisk==ProcessPar.VolOut(iWhisk)&TBLall.PowerZero==0&TBLall.Speed<SpeedTh,:);
tempAll_sub=TBLall_sub(TBLall_sub.Whisk==ProcessPar.VolOut(iWhisk)&TBLall_sub.PowerZero==0&TBLall_sub.Speed<SpeedTh,:);
% 
tempNonT=TBLallNontarget(TBLallNontarget.Whisk==ProcessPar.VolOut(iWhisk)&TBLallNontarget.PowerZero==0&abs(TBLallNontarget.Speed)<SpeedTh,:);
tempAll=TBLall(TBLall.Whisk==ProcessPar.VolOut(iWhisk)&TBLall.PowerZero==0&abs(TBLall.Speed)<SpeedTh,:);
tempAll_sub=TBLall_sub(TBLall_sub.Whisk==ProcessPar.VolOut(iWhisk)&TBLall_sub.PowerZero==0&abs(TBLall_sub.Speed<SpeedTh),:);


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

tempNonT=TBLallNontarget(TBLallNontarget.Whisk==ProcessPar.VolOut(iWhisk)&TBLallNontarget.PowerZero==0&TBLallNontarget.Speed<SpeedTh,:);
tempAll=TBLall(TBLall.Whisk==ProcessPar.VolOut(iWhisk)&TBLall.PowerZero==0&TBLall.Speed<SpeedTh,:);
tempAll_sub=TBLall_sub(TBLall_sub.Whisk==ProcessPar.VolOut(iWhisk)&TBLall_sub.PowerZero==0&TBLall_sub.Speed<SpeedTh,:);




for iGroup=1:3
    dataAll{iGroup,iWhisk}=tempAll.WhiskerScore(tempAll.Group==iGroup);
    dataNonT{iGroup,iWhisk}=tempNonT.WhiskerScore(tempNonT.Group==iGroup);
    dataAll_sub{iGroup,iWhisk}=tempAll_sub.WhiskerScore(tempAll_sub.Group==iGroup);

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
A=TBLallNontarget(TBLallNontarget.Whisk==1&TBLallNontarget.PowerZero==0&abs(TBLallNontarget.Speed)<SpeedTh,:);
B=TBLall(TBLall.Whisk==1&TBLall.PowerZero==0&abs(TBLall.Speed)<SpeedTh,:);

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
pcaN=20;
ScoreLabel={'SpeedScore','WhiskScore'}
DataTable=[];
for iFOV=1:length(IndexFOVNeed)
tblFOV=SLMGroupTableTrial(SLMGroupTableTrial.Session==iFOV,:);
TrialID=unique(tblFOV.TrialID);
tbltemp=SLMGroupTableTrial(SLMGroupTableTrial.Session==iFOV&SLMGroupTableTrial.TrialID==TrialID(1),:);
% [W,H] = nnmf(NeuroTrace{iFOV}{iData}(tbltemp.TargetCell==0,SponInd),pcaN);

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


pcaN=10;
rSpeed=corr(coeff(ExcludeStimInd,1:pcaN),nanmean(BehTrace(iFOV).Speed(ExcludeStimInd,:),2),'rows','complete')

maxLag=10;
clear rStim c
for PCI = 1:pcaN
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
























