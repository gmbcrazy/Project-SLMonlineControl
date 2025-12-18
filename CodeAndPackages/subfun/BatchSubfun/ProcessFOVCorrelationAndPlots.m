function CorrResults = ProcessFOVCorrelationAndPlots(CaData, SpeedAll, StimAll, Suite2pTable, SessFileTable, ResultFolder)
% Handles all correlation, cross-correlation, sorting, and plotting tasks for one FOV

% Prepare neural data
iscell = find(CaData.iscell(:,1)>0.99);
deltaFoF = double(F2deltaFoF(CaData.F, CaData.Fneu, CaData.ops.fs));
deltaFoF = deltaFoF(:, iscell);
spks = double(CaData.spks(iscell,:));
NData = {deltaFoF', spks};
Nlabel = {'DeltaF', 'Spks'};

% Correlation/xcorr settings
XTimesStdTh = 3;
MinInterVal = 20;
maxLag = 10;

rSpeed = [];
pSpeed = [];
rStim = [];

nPlane = length(unique(CaData.CellPlaneID));

% --- Exclude stim periods for spontaneous analysis ---
InitialSessInd = find(isnan(Suite2pTable.markCycle)==1);
nFrameInitialSession=Suite2pTable.Suite2pTiffNum(InitialSessInd(1))/nPlane;
InitialInd=1:nFrameInitialSession;

InitialDataI = [];
frameGap = 100;

ExcludeStimInd = [];
for i=1:length(InitialSessInd)
    tempF = double(StimAll(InitialInd,1));
    t1 = min(find(tempF>1))-frameGap;
    t2 = max(find(tempF>1))+frameGap;
    ExcludeStimIndraw = [1:t1 t2:nFrameInitialSession]+nFrameInitialSession*(i-1);
    ExcludeStimInd = [ExcludeStimInd, ExcludeStimIndraw];
end

ExcludeStimInd=sort(ExcludeStimInd);
% plot(ExcludeStimInd)


InitialInd=1:sum(Suite2pTable.Suite2pTiffNum(InitialSessInd))/nPlane;

% --- Main calculation for each data type ---
for iData = 1:length(NData)

    %First spontonouse behavior period
    for iPlane = 1:nPlane
        I1 = find(CaData.CellPlaneID==iPlane);
        [rSpeed(I1,1,iData), pSpeed(I1,1,iData)] = corr(NData{iData}(I1,ExcludeStimInd)', SpeedAll(ExcludeStimInd,iPlane), 'type','Spearman','rows','pairwise');
    end
    
    StimAllT = double(StimAll > 1);
    for iCell = 1:length(iscell)
        iPlane = CaData.CellPlaneID(iCell);
        [c(:,iCell), lags] = xcorr(NData{iData}(iCell,InitialInd)', StimAllT(InitialInd,iPlane), maxLag, 'coeff');
        PostI = find(lags >= 0);
        [~, i1] = max(abs(c(PostI,iCell)));
        rStim(iCell, 1,iData) = c(PostI(i1),iCell);
    end
    % --- Plotting ---
    SaveCorrFigure(NData{iData}, rStim(:, 1,iData), mean(StimAllT(InitialInd,:),2)>0.1,InitialInd, Nlabel{iData}, ResultFolder,'Stim',[91 20 212]/255);
    SaveCorrFigure(NData{iData}, rSpeed(:, 1,iData), mean(SpeedAll(InitialInd,:),2), InitialInd, Nlabel{iData}, ResultFolder,'Speed',[255 51 153]/255);

    %later spontonouse behavior + SLM group stimuli period
    %Stim TTL for during awake
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

   P.xLeft=0.08;        %%%%%%Left Margin
   P.xRight=0.02;       %%%%%%Right Margin
   P.yTop=0.04;         %%%%%%Top Margin
   P.yBottom=0.12;      %%%%%%Bottom Margin
   P.xInt=0.08;         %%%%%%Width-interval between subplots
   P.yInt=0.02;         %%%%%%Height-interval between subplots

% iData=1;
figure;
regParam.Color=[0.1 0.1 0.1];
regParam.Marker='o';
regParam.MarkerSize=3;
regParam.Rtype='pearson';
regParam.xLim=[-0.3 0.6];
regParam.yLim=[-0.3 0.6];
regParam.xLabel='SpeedCorr Pre';
regParam.yLabel='SpeedCorr Post';
subplotLU(1,2,1,1,P);
[OutPut,r,p]=LuPairRegressPlot(rSpeed(:,1,iData),rSpeed(:,2,iData),regParam);
plot([-0.3 0.6],[-0.3 0.6],'k:');hold on;
% set(gca,'xlim',[-0.3 0.6],'ylim',[-0.3 0.6])

regParam.xLabel='StimCorr Pre';
regParam.yLabel='StimCorr Awake';
regParam.xLim=[-0.05 0.2];
regParam.yLim=[-0.05 0.2];
subplotLU(1,2,1,2,P);
[OutPut,r,p]=LuPairRegressPlot(rStim(:,1,iData),rStim(:,2,iData),regParam);  
plot([-0.3 0.6],[-0.3 0.6],'k:');hold on;

% % subplotLU(1,3,1,3,P)
% % regParam.xLabel='StimCorr Pre';
% % regParam.yLabel='StimCorr Ane';
% % 
% % [OutPut,r,p]=LuPairRegressPlot(rStim(:,1,iData),rStim(:,3,iData),regParam)    
% % plot([-0.3 0.6],[-0.3 0.6],'k:');hold on;
papersizePX=[0 0 20 10];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
saveas(gcf,[ResultFolder 'Beh' Nlabel{iData} 'CorrPrePosSLM'],'png');
close all



    

end

CorrResults.rSpeed = rSpeed;
CorrResults.rStim = rStim;



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


% papersizePX=[0 0 20 10];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
saveas(gcf,[ResultFolder 'Beh' Nlabel{1} Nlabel{2} 'SponSLM'],'png');
close all




end
