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
InitialInd = find(isnan(Suite2pTable.markCycle)==1);
InitialDataI = [];
frameGap = 100;

ExcludeStimInd = [];
for i=1:length(InitialInd)
    tempF = double(StimAll(:,1));
    t1 = min(find(tempF>1))-frameGap;
    t2 = max(find(tempF>1))+frameGap;
    ExcludeStimIndraw = [1:t1 t2:Suite2pTable.Suite2pTiffNum(InitialInd(i))/nPlane];
    ExcludeStimInd = [ExcludeStimInd, ExcludeStimIndraw];
end
InitialInd=1:sum(Suite2pTable.Suite2pTiffNum(InitialInd))/nPlane;

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
    SaveCorrFigure(NData{iData}, rStim(:, 1,iData), mean(StimAllT(InitialInd,:),2)>0.1,InitialInd, Nlabel{iData}, ResultFolder,'Stim');
    SaveCorrFigure(NData{iData}, rSpeed(:, 1,iData), mean(SpeedAll(InitialInd,:),2), InitialInd, Nlabel{iData}, ResultFolder,'Speed');

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



    

end

CorrResults.rSpeed = rSpeed;
CorrResults.rStim = rStim;

end
