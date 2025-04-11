clear all

WorkFolder='E:\LuSLMOnlineTest\SL0855-Emx1G6CII-AAV9CAMKII\03042025\';
% ConfigFolder='C:\Users\zhangl33\Projects\Project-SLMonlineControl\config\';
% ConfigFile='CurrentSLMsetting.yml';%<----------------------------------------------------------------------------------Edit, configuration file
% [~,~,~,CaData,CaDataPlane,stat,yaml,confSet]=ROIToXYZ(WorkFolder,ConfigFile);
ProcessFolder = Get_ExpDataFolder(WorkFolder, 'SpeedStimEdgeExc', {'Data','AllIncluded','DataSum','.gpl','.xml'})
DataFolder=[ProcessFolder 'Data\'];
load('C:\Users\zhangl33\Projects\Project-SLMonlineControl\CodeAndPackages\subfun\Color\colorMapPN3.mat','colorMapPN1')

ResultFolder=[ProcessFolder 'Results\'];
mkdir(ResultFolder)
ResultFolder=[ResultFolder 'Step1\'];
mkdir(ResultFolder)
load([DataFolder 'TableForSuite2p.mat'])
motionTh=5;

confSet=SLMPosInfo.confSetFinal;
confSet.save_path0=DataFolder;
Zdepth=confSet.scan_Z+confSet.ETL;
nPlane=length(Zdepth);


[FileIDList,i1]=unique(Suite2pTable.FileID);
for i=1:length(FileIDList)
    [fSpeed{i},fStim{i},timeStampCa{i},FrameTS{i},fVideo{i},VideoStartFrameTime{i}]=PV_VolExtract_MultiCyc(confSet,FileIDList(i));
end

SessFileTable=Suite2pTable(i1,:);
tempTiffmark=[];
suite2pInd=[];

CumTiffNum=cumsum(SessFileTable.Suite2pTiffNum);
CumFrameNum=[0;CumTiffNum(1:end-1)]/nPlane;
clear i1
for i=1:length(FileIDList)
    i2=find(Suite2pTable.FileID==FileIDList(i));
    for j=1:length(i2)
        AddTemp=Suite2pTable.markCycle(i2(j))-(j-1)*RemoveFrame;
        tempTiffmark=[tempTiffmark;AddTemp];
        suite2pInd=[suite2pInd;AddTemp+CumFrameNum(i)];
    end
end
clear i2

Suite2pTable.PostSLMTiffCycle=tempTiffmark;
Suite2pTable.suite2pInd=suite2pInd;

[NeuronPos3D,NeuronPos3DRaw,CaData,CaDataPlane,Neuronstat]=Extract_Suite2p(confSet);

SpeedAll=[];
StimAll=[];
for i=1:length(FileIDList)
    i2=find(Suite2pTable.FileID==FileIDList(i));

    frameN=SessFileTable.Suite2pTiffNum(i)/nPlane;
    if ~isempty(fSpeed{i})
        InValid=[];
        for j=1:RemoveFrame  
            InValid=[InValid;Suite2pTable.markCycle(i2)+j-1];
        end
        InValid(isnan(InValid))=[];
        tempSpeed=fSpeed{i};
        tempStim=fStim{i};
        tempSpeed(InValid,:)=[];
        tempStim(InValid+1,:)=[];
    else
        tempSpeed=zeros(frameN,nPlane)+NaN;
        tempStim=tempSpeed;

    end
    SpeedAll=[SpeedAll;tempSpeed];
    StimAll=[StimAll;tempStim];
   
    clear tempSpeed tempStim
end





SLMPos3D=SLMTestInfo.Pos3Dneed;
SLMGroup=SLMTestInfo.FunScore(:,1);





DistTh=6;
[SLMtarget,SLMtargetCellDist]=SLMtargetMatchCell(SLMPos3D,NeuronPos3D,DistTh)
PointList1=find(SLMtarget>0);

SLMtargetTable=zeros(size(Suite2pTable,1),1);
SLMtargetTableGroup=zeros(size(Suite2pTable,1),1);

for iTarget=1:length(SLMtarget)
    SLMtargetTable(find(Suite2pTable.Point==iTarget))=SLMtarget(iTarget);
    SLMtargetTableGroup(find(Suite2pTable.Point==iTarget))=SLMGroup(iTarget);
end
Suite2pTable.PointTargetCell=SLMtargetTable;
Suite2pTable.PointTargetCellGroup=SLMtargetTableGroup;


ismember(SLMPosInfo.FinalPos3D,SLMPos3D)

[SLMFinalInSLMtest,~]=SLMtargetMatchCell(SLMPosInfo.FinalPos3D,SLMPos3D,0.1);
% [SLMPosInfo.FinalPos3D(1:2,:) SLMPos3D(SLMFinalInSLMtest(1:2),:)]


[TargetCellList,ia]=unique(SLMtargetTable(SLMtargetTable>0));
temp=SLMtargetTableGroup(SLMtargetTable>0);
TargetCellListFunGroup=temp(ia);

GroupTargetCell={};
for iFun=1:length(SLMPosInfo.Group)
    [temp,~]=SLMtargetMatchCell(SLMPosInfo.FinalPos3D(SLMPosInfo.Group(iFun).Indices,:),NeuronPos3D,DistTh);
    GroupTargetCell{iFun}=temp(temp>0);
end


deltaFoF=double(F2deltaFoF(CaData.F,CaData.Fneu,CaData.ops.fs));
iscell=find(CaData.iscell(:,1)>0.99);
deltaFoF=deltaFoF(:,iscell);
spks=double(CaData.spks(iscell,:));



% preSLM=10;
% postSLM=3;

PSTHparam.PreSLMCal=10;        %<----------------------------------------------------------------------------------Edit,Frame # before SLM to calculate baseline map
PSTHparam.PostSLMCal=15;        %<----------------------------------------------------------------------------------Edit,Frame # after SLM to calculate responsive map
% PSTHparam.YLim=[-50 600];       % %<-----------------------------------------For Suite2p based ROI signal only method
PSTHparam.pTh=0.05;             % %<-----------------------------------------For Suite2p based ROI signal only method
PSTHparam.TestMethod='ranksum'; % %<-----------------------------------------For Suite2p based ROI signal only method
PSTHparam.MPFrameJump=2; % %<-----------------------------------------For Suite2p based ROI signal only method


TimBinFrame=-PSTHparam.PreSLMCal:PSTHparam.PostSLMCal-1;        %<----------------------------------------------------------------------------------Edit,Frame # before SLM to calculate baseline map



NData={deltaFoF' spks};
Nlabel={'DeltaF','Spks'};


% rSpeedInitial=corr(deltaFoF,nanmean(SpeedAll,2),'rows','complete');
% rSpeedInitial=corr(deltaFoF,nanmean(SpeedAll,2),'rows','complete');

Suite2pTable.Suite2pTiffNum(isnan(Suite2pTable.markCycle))/nPlane;
InitialInd=find(isnan(Suite2pTable.markCycle)==1);

InitialDataI=[];
frameGap=100;             %%%Frame before/after 1st/last whiskerStim were consider to exclude to calculate neural-speed correlation.
i=1
tempF=double(fStim{InitialInd(i)}(:,1));
t1=min(find(tempF>1))-frameGap;
t2=max(find(tempF>1))+frameGap;
ExcludeStimIndraw=[1:t1 t2:Suite2pTable.Suite2pTiffNum(InitialInd(i))/nPlane];  

ExcludeStimInd=ExcludeStimIndraw;
for i=2:length(InitialInd)
    ExcludeStimInd=[ExcludeStimInd,ExcludeStimInd+Suite2pTable.Suite2pTiffNum(InitialInd(i))/nPlane];
end
InitialInd=1:sum(Suite2pTable.Suite2pTiffNum(InitialInd))/nPlane;


XTimesStdTh = 3;
MinInterVal = 20;
maxLag = 10;
clear rSpeed pSpeed rStim;

for iData =1:length(NData)

for iPlane=1:nPlane
    I1=find(CaData.CellPlaneID==iPlane);
    [rSpeed(I1,1,iData),pSpeed(I1,1,iData)]=corr(NData{iData}(I1,ExcludeStimInd)',SpeedAll(ExcludeStimInd,iPlane),'type','Spearman','rows','pairwise');
end

clear c
StimAllT=double(StimAll>1);


for iCell = 1:length(iscell)
    iCell;
    iPlane=CaData.CellPlaneID(iCell);
    [c(:,iCell), lags] = xcorr(NData{iData}(iCell,InitialInd)', StimAllT(InitialInd,iPlane), maxLag, 'coeff');
    PostI = find(lags >= 0);
    [~, i1] = max(abs(c(PostI,iCell)));
    rStim(iCell, 1,iData) = c(PostI(i1),iCell);
end

%%Stim TTL for during anesthesia
StimAllTAne=double(StimAllT);


[~,rankStimTemp]=sort(rStim(:, 1,1),'descend')
figure;
imagesc(NData{1}(rankStimTemp,InitialInd))


% AneI=Suite2pTable.suite2pInd(Suite2pTable.AwakeState==2&Suite2pTable.PowerZero==0)
% StimAllTAne(AneI,:)=3;
% StimAllTAne=double(StimAllTAne==3);
StimSLM0TAne=double(StimAllT);
t1=Suite2pTable.suite2pInd(Suite2pTable.AwakeState==2&Suite2pTable.PowerZero==1);
StimSLM0TAne(t1,:)=3;
StimSLM0TAne=double(StimSLM0TAne==3);

%%Stim TTL for during awake
StimSLM0TWake=double(StimAllT);
t1=Suite2pTable.suite2pInd(Suite2pTable.AwakeState==1&Suite2pTable.PowerZero==1);
StimSLM0TWake(t1,:)=3;
StimSLM0TWake=double(StimSLM0TWake==3);




FunExpStartI=min(find(SessFileTable.Group>0))-1;
FunExpEndI=max(find(SessFileTable.AwakeState==1));
FunExpInd=sum(SessFileTable.Suite2pTiffNum(1:FunExpStartI))/nPlane+1:sum(SessFileTable.Suite2pTiffNum(1:FunExpEndI))/nPlane;
for iPlane=1:nPlane
    I1=find(CaData.CellPlaneID==iPlane);
    [rSpeed(I1,2,iData),pSpeed(I1,2,iData)]=corr(NData{iData}(I1,FunExpInd)',SpeedAll(FunExpInd,iPlane),'type','Spearman','rows','pairwise');
end

StimAllT=StimAll>1;
clear c;
for iCell = 1:length(iscell)
    iCell;
    iPlane=CaData.CellPlaneID(iCell);
    [c(:,iCell), lags] = xcorr(NData{iData}(iCell,FunExpInd)', StimSLM0TWake(FunExpInd,iPlane), maxLag, 'coeff');
    PostI = find(lags >= 0);
    [~, i1] = max(abs(c(PostI,iCell)));
    rStim(iCell, 2,iData) = c(PostI(i1),iCell);
end



AneExpStartI=min(find(SessFileTable.AwakeState==2))-1;
AneExpEndI=max(find(SessFileTable.AwakeState==2));
AneExpInd=sum(SessFileTable.Suite2pTiffNum(1:AneExpStartI))/nPlane+1:sum(SessFileTable.Suite2pTiffNum(1:AneExpEndI))/nPlane;

clear c;
for iCell = 1:length(iscell)
    iCell;
    iPlane=CaData.CellPlaneID(iCell);
    [c(:,iCell), lags] = xcorr(NData{iData}(iCell,AneExpInd)', StimSLM0TAne(AneExpInd,iPlane), maxLag, 'coeff');
    PostI = find(lags >= 0);
    [~, i1] = max(abs(c(PostI,iCell)));
    rStim(iCell,3,iData) = c(PostI(i1),iCell);
end





   P.xLeft=0.08;        %%%%%%Left Margin
   P.xRight=0.02;       %%%%%%Right Margin
   P.yTop=0.04;         %%%%%%Top Margin
   P.yBottom=0.12;      %%%%%%Bottom Margin
   P.xInt=0.08;         %%%%%%Width-interval between subplots
   P.yInt=0.02;         %%%%%%Height-interval between subplots


figure;
regParam.Color=[0.1 0.1 0.1];
regParam.Marker='o';
regParam.MarkerSize=3;
regParam.Rtype='pearson';
regParam.xLim=[-0.3 0.6];
regParam.yLim=[-0.3 0.6];
regParam.xLabel='SpeedCorr Pre';
regParam.yLabel='SpeedCorr Post';
subplotLU(1,3,1,1,P)
[OutPut,r,p]=LuPairRegressPlot(rSpeed(:,1,iData),rSpeed(:,2,iData),regParam);
plot([-0.3 0.6],[-0.3 0.6],'k:');hold on;
% set(gca,'xlim',[-0.3 0.6],'ylim',[-0.3 0.6])

regParam.xLabel='StimCorr Pre';
regParam.yLabel='StimCorr Awake';
regParam.xLim=[-0.05 0.2];
regParam.yLim=[-0.05 0.2];
subplotLU(1,3,1,2,P)
[OutPut,r,p]=LuPairRegressPlot(rStim(:,1,iData),rStim(:,2,iData),regParam)    
plot([-0.3 0.6],[-0.3 0.6],'k:');hold on;

subplotLU(1,3,1,3,P)
regParam.xLabel='StimCorr Pre';
regParam.yLabel='StimCorr Ane';

[OutPut,r,p]=LuPairRegressPlot(rStim(:,1,iData),rStim(:,3,iData),regParam)    
plot([-0.3 0.6],[-0.3 0.6],'k:');hold on;
papersizePX=[0 0 30 10];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
saveas(gcf,[ResultFolder 'Beh' Nlabel{iData} 'CorrNoSLM'],'png');
close all
end


figure;
regParam.Color=[0.1 0.1 0.1];
regParam.Marker='o';
regParam.MarkerSize=3;
regParam.Rtype='pearson';
regParam.xLim=[-0.3 0.6];
regParam.yLim=[-0.3 0.6];
regParam.xLabel='SpeedCorr Pre DeltaF';
regParam.yLabel='SpeedCorr Pre spks';
subplotLU(1,2,1,1,P)
[OutPut,r,p]=LuPairRegressPlot(rSpeed(:,1,1),rSpeed(:,1,2),regParam);
plot([-0.3 0.6],[-0.3 0.6],'k:');hold on;
% set(gca,'xlim',[-0.3 0.6],'ylim',[-0.3 0.6])

regParam.xLabel='StimCorr Pre DeltaF';
regParam.yLabel='StimCorr Pre spks';
regParam.xLim=[-0.05 0.2];
regParam.yLim=[-0.05 0.2];
subplotLU(1,2,1,2,P)
[OutPut,r,p]=LuPairRegressPlot(rStim(:,1,1),rStim(:,1,2),regParam)    
plot([-0.3 0.6],[-0.3 0.6],'k:');hold on;


papersizePX=[0 0 20 10];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
saveas(gcf,[ResultFolder 'Beh' Nlabel{1} Nlabel{2} 'CorrNoSLM'],'png');
close all




PVpower=xmlPower2PVpower(SLMTestInfo.confSet.UncagingLaserPower);
PVpower=intersect(round(PVpower),SLMInfoTable.UncagingLaserPower);


AlignedSpeed=[];
AlignedStim=[];

for iData =1:length(NData)

tempNData=NData{iData};
tempNData=AmpNormalizeRow(tempNData,[0 100]);



SLMInfoTable=Suite2pTable(Suite2pTable.suite2pInd>0,:);
AlignedtempNData=[];

for i=1:size(SLMInfoTable,1)
    I0=SLMInfoTable.suite2pInd(i);
    s1=I0-PSTHparam.PreSLMCal:I0-1;
    s2=I0:I0+PSTHparam.PostSLMCal-1;
    AlignedtempNData(:,:,i)=tempNData(:,[s1 s2])-repmat(nanmean(tempNData(:,s1),2),1,length(s1)+length(s2));

    if iData==1
    AlignedSpeed(:,i)=mean(SpeedAll([s1 s2],:),2);
    AlignedStim(:,i)=mean(StimAll([s1 s2],:),2);
    end
end

TrialThNum=3;

PSTHparam.TestStepFrame=3;
[AlignedtempNData,TargetInfoTable,CellResponse,statCellRes,TargetResponse,TargetCellResP,TargetCellResR,CellSampleN]=Aligned_FromSuite2p(NData{iData},TargetCellList,iscell,Suite2pTable,PVpower,PSTHparam);
if iData==1

FDR=0.05;
[h, crit_p, adj_p]=fdr_bh(TargetCellResP((~isnan(TargetCellResP))&CellSampleN>=TrialThNum),FDR,'pdep');crit_p


PowerTargetI=zeros(length(TargetCellList),1);
SuccTarget=PowerTargetI;
for iCell=1:length(TargetCellList)
    for iPower=1:length(PVpower)
        if TargetCellResP(iCell,iPower)<=crit_p&&TargetCellResR(iCell,iPower)>0&&CellSampleN(iCell,iPower)>=TrialThNum
           PowerTargetI(iCell)=iPower;
           SuccTarget(iCell)=1;
           SuccAmp(iCell)=TargetCellResR(iCell,iPower);
        end
    end
end
sum(SuccTarget)

pAll=[];
for iCell=1:length(TargetCellList)
    if SuccTarget(iCell)
       temp=statCellRes(iCell,PowerTargetI(iCell)).p(:);
       % temp(TargetCellList(iCell))=[];
       pAll=[pAll;temp];
    end
end
FDR=0.1;

[h, crit_pAll, adj_p]=fdr_bh(pAll,FDR,'pdep');crit_pAll

end



close all

Param.PlotType=3
Param.statisP=0;
Param.LegendShow=0;
Param.Legend=[]


   P.xLeft=0.06;        %%%%%%Left Margin
   P.xRight=0.02;       %%%%%%Right Margin
   P.yTop=0.02;         %%%%%%Top Margin
   P.yBottom=0.1;      %%%%%%Bottom Margin
   P.xInt=0.02;         %%%%%%Width-interval between subplots
   P.yInt=0.01;         %%%%%%Height-interval between subplots

figure;
for iCell = 1:length(TargetCellList)
    for iPower=1:length(PVpower)
        if ~isempty(TargetResponse(iCell,iPower))&CellSampleN(iCell,iPower)>=TrialThNum

            % TargetCellList(iCell)
            subplotLU(length(TargetCellList),length(PVpower),iCell,iPower,P)
            RateHist_GroupPlot(TimBinFrame+0.5,TargetResponse(iCell,iPower),[0.1 0.1 0.1],Param)
            hold on;
            text(-10,0.1,['C' num2str(iCell) 'n = ' num2str(CellSampleN(iCell,iPower))])
            set(gca,'xlim',[-PSTHparam.PreSLMCal PSTHparam.PostSLMCal],'xtick',[-PSTHparam.PreSLMCal:5:PSTHparam.PostSLMCal])
            ylabel(['Target' num2str(iCell)]);
        end
    end

end
papersizePX=[0 0 length(PVpower)*4 5*length(TargetCellList) ];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
saveas(gcf,[ResultFolder 'AllSLMTargetResponse' Nlabel{iData}],'png');
close all



figure;
for iCell = 1:length(TargetCellList)
    for iPower=1:length(PVpower)
        if ~isempty(CellResponse{iCell,iPower})&CellSampleN(iCell,iPower)>=TrialThNum
            % TargetCellList(iCell)
            subplotLU(length(TargetCellList),length(PVpower),iCell,iPower,P)
            imagesc(TimBinFrame+0.5,1:length(iscell),CellResponse{iCell,iPower});
            colormap(colorMapPN1);
            set(gca,'clim',[-0.1 0.1],'ylim',[0 length(iscell)+1]);
            hold on;
            plot(TimBinFrame(1),TargetCellList(iCell)+0.5,'g>');
            set(gca,'xlim',[-PSTHparam.PreSLMCal PSTHparam.PostSLMCal],'xtick',[-PSTHparam.PreSLMCal:5:PSTHparam.PostSLMCal])
        end
    end

end
papersizePX=[0 0 length(PVpower)*3 5*length(TargetCellList) ];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
saveas(gcf,[ResultFolder 'AllSLMTargetResponseMap' Nlabel{iData}],'png');
close all




AlignedNData{iData}=AlignedtempNData;




end







save([ResultFolder 'Step1Meta.mat']);

for iData=1:length(NData)
clear Data1 Data2
for iFun=1:length(SLMPosInfo.Group)
    Data1{iFun}=rSpeed(TargetCellList(TargetCellListFunGroup==iFun),1,iData);
    Data1{iFun+length(SLMPosInfo.Group)}=rSpeed(TargetCellList(TargetCellListFunGroup==iFun),2,iData);

    Data2{iFun}=rStim(TargetCellList(TargetCellListFunGroup==iFun),1,iData);
    Data2{iFun+length(SLMPosInfo.Group)}=rStim(TargetCellList(TargetCellListFunGroup==iFun),2,iData);
    Data2{iFun+length(SLMPosInfo.Group)*2}=rStim(TargetCellList(TargetCellListFunGroup==iFun),3,iData);

end

x1=1:length(SLMPosInfo.Group);x1=[x1 x1+length(SLMPosInfo.Group)+1];
x2=1:length(SLMPosInfo.Group);x2=[x2 x2+length(SLMPosInfo.Group)+1 x2+length(SLMPosInfo.Group)*2+2];
GroupLabel={'L.','S.','N.'};

nGroup=length(SLMPosInfo.Group);



% % subplotLU(2,1,1,1)
% % stats=ErrorBoxPlotLU(x1,Data1,repmat(GroupColor,2,1),[ResultFolder 'SpeedGroup']);
% % LuLegend([9 9 9;0.5 0.4 0.3;10 10 10;0.5 0.4 0.3],0,GroupLabel,GroupColor,8);
% % set(gca,'xlim',[0 12],'ylim',[-0.2 1]);
% % subplotLU(2,1,2,1)
% % stats=ErrorBoxPlotLU(x2,Data2,repmat(GroupColor,3,1),[ResultFolder 'SpeedGroup']);
% % set(gca,'xlim',[0 12]);

GroupColor=[247 150 111;239 109 249;121 247 111]/255;

   GroupPair.CorrName='fdr';
   GroupPair.Q=0.1;

   
   GroupPair.SignY=1;
   GroupPair.Plot=1;
   GroupPair.Std=1;      %%%%%%%%%using standard deviation as errorbar
   GroupPair.SamplePlot=1; %%%%%%%%%Plot Individual Sample Point
   GroupPair.SamplePairedPlot=1; %%%%%%%%%Dash line for paired comparison sample
   GroupPair.LimY=[0 GroupPair.SignY*1.2];
   GroupPair.Marker={'o'};
   GroupPair.GroupName=repmat(GroupLabel,1,2);
   % GroupID=repmat(1:nGroup,1,2);

P.xInt=0.1;
P.xLeft=0.1;

figure;
subplotLU(1,2,1,1,P)
   p1=[1 1 2;2 3 3];p1=[p1 p1+nGroup];
   p2=[1:nGroup];p2=[p2;p2+nGroup];
   GroupPair.Pair=[p1 p2];
GroupPair.GroupName=repmat(GroupLabel,1,2);
   GroupPair.SignY=0.5;

ErrorBarPlotLU(x1,Data1,[],repmat(GroupColor,2,1),2,1,[ResultFolder Nlabel{iData} 'SpeedGroup.txt'],GroupPair,repmat(1:nGroup,1,2));
LuLegend([9 9 9;0.5 0.4 0.3;10 10 10;0.5 0.4 0.3],0,GroupLabel,GroupColor,8);
set(gca,'xlim',[0 12],'ylim',[-0.2 0.6],'ytick',[-0.2:0.2:0.6],'xtick',[2 6],'xticklabel',{'Spon 1','Spon 2'});
ylabel('Speed Corr')

subplotLU(1,2,1,2,P)
   p1=[1 1 2;2 3 3];p1=[p1 p1+nGroup p1+nGroup*2];
   p2=[1:nGroup];p2=[p2 p2 p2+nGroup;p2+nGroup p2+nGroup*2 p2+nGroup*2];
   GroupPair.Pair=[p1 p2];
   GroupPair.SignY=0.25;
   GroupID=repmat(1:nGroup,1,3);
   GroupPair.GroupName=repmat(GroupLabel,1,3);

ErrorBarPlotLU(x2,Data2,[],repmat(GroupColor,3,1),2,1,[ResultFolder Nlabel{iData} 'StimGroup.txt'],GroupPair,repmat(1:nGroup,1,3));
set(gca,'xlim',[0 12],'ylim',[-0.05 0.3],'ytick',[-0.05:0.05:0.3],'xtick',[2 6 10],'xticklabel',{'Spon 1','Spon 2','Ane'});
ylabel('Stim Corr')

papersizePX=[0 0 22 10];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
saveas(gcf,[ResultFolder 'BehCorrCompare' Nlabel{iData}],'png');
close all

end

