function OfflineCombineProcess_OneFOV(FOVtemp, Suite2pDataKeywords, suite2pFOVPathLocal,PSTHparam)
% Full modularized pipeline for all analysis/plotting on one FOV.
% Input:
%   FOVtemp: one element of FOVUpdate struct array
%   Suite2pDataKeywords: string, e.g. 'awakeRefSpon'
%   suite2pFOVPathLocal: cell array, same as in the original script
%   PSTHparam.PreSLMCal = 10; PSTHparam.PostSLMCal = 3;
%   PSTHparam.pTh = 0.05; PSTHparam.TestMethod = 'ranksum';
%   PSTHparam.MPFrameJump = 2;
%% --- Setup paths and load data ---
DataFolder = FOVtemp.DataFolder;
ProcessFolder = [fileparts(FOVtemp.DataFolder(1:end-1)) '\'];
Temp1 = [fileparts(ProcessFolder(1:end-1)) '\'];
WorkFolder = [fileparts(Temp1(1:end-1)) '\'];

SLMPosInfo = load([ProcessFolder 'SLMFunGroup.mat']);
SLMTestInfo = load([ProcessFolder 'SLMIncludedIndFromIscell.mat']);
confSet = SLMPosInfo.confSetFinal;
Zdepth = confSet.scan_Z + confSet.ETL;
Pos3DFun = SLMPosInfo.FinalPos3D;
Group = SLMPosInfo.Group;
for iGroup = 1:length(Group)
    Pos3DGroup{iGroup} = Pos3DFun(Group(iGroup).Indices,:);
end
confSet.save_path0 = DataFolder;
nPlane = length(Zdepth);


resultPaths = findAllFoldersKeyWords(suite2pFOVPathLocal, Suite2pDataKeywords);
if length(resultPaths) > 1
    disp('more than one suite2p data folder!');
elseif length(resultPaths) == 1
    suite2pDataFolder = resultPaths{1};
else
    suite2pDataFolder = DataFolder;
end
ResultFolder = [suite2pDataFolder 'result\'];
mkdir(ResultFolder); ResultFolder = [ResultFolder 'Step1Basic\']; mkdir(ResultFolder);

Suite2pTable = FOVtemp.Suite2pTable;
Suite2pTable = Suite2pTable(Suite2pTable.AwakeState<2,:);  % Exclude anesthesia data
[FileIDList,i1] = unique(Suite2pTable.FileID);




RemoveFrame = PSTHparam.MPFrameJump;
SessFileTable = Suite2pTable(i1,:);
tempTiffmark = [];
suite2pInd = [];
CumTiffNum = cumsum(SessFileTable.Suite2pTiffNum);
CumFrameNum = [0; CumTiffNum(1:end-1)]/nPlane;
for i = 1:length(FileIDList)
    i2 = find(Suite2pTable.FileID == FileIDList(i));
    for j = 1:length(i2)
        AddTemp = Suite2pTable.markCycle(i2(j))-(j-1)*RemoveFrame;
        tempTiffmark = [tempTiffmark; AddTemp];
        suite2pInd = [suite2pInd; AddTemp + CumFrameNum(i)];
    end
end
Suite2pTable.PostSLMTiffCycle = tempTiffmark;
Suite2pTable.suite2pInd = suite2pInd;

if ~exist([ResultFolder 'BehAll.mat'])
    for i = 1:length(FileIDList)
        [fSpeed{i}, fStim{i}, timeStampCa{i}, FrameTS{i}, fVideo{i}, VideoStartFrameTime{i}] = PV_VolExtract_MultiCyc(confSet, FileIDList(i));
    end
    SpeedAll = [];
    StimAll = [];
    for i = 1:length(FileIDList)
        i2 = find(Suite2pTable.FileID == FileIDList(i));
        frameN = SessFileTable.Suite2pTiffNum(i)/nPlane;
        if ~isempty(fSpeed{i})
            InValid = [];
            for j = 1:RemoveFrame
                InValid = [InValid; Suite2pTable.markCycle(i2) + j - 1];
            end
            InValid(isnan(InValid)) = [];
            tempSpeed = fSpeed{i}; tempStim = fStim{i};
            tempSpeed(InValid,:) = [];
            tempStim(InValid+1,:) = [];
        else
            tempSpeed = zeros(frameN, nPlane) + NaN;
            tempStim = tempSpeed;
        end
        SpeedAll = [SpeedAll; tempSpeed];
        StimAll = [StimAll; tempStim];
    end
    save([ResultFolder 'BehAll.mat'],'SpeedAll','StimAll');
else
    load([ResultFolder 'BehAll.mat']);
end

confSet.save_path0 = suite2pDataFolder;

[NeuronPos3D, NeuronPos3DRaw, CaData, CaDataPlane, Neuronstat] = Extract_Suite2p(confSet);
TempconfSet = confSet; TempconfSet.save_path0 = WorkFolder;
[~,~,OnlineCaData,~,~] = Extract_Suite2p(TempconfSet);
[~, ~, ~, cellBoundary, ~] = Suite2pCellIDMapFromStat(CaData.statCell, [512 512]);



%% --- 1. Correlation & Speed/Stim/Neural plots ---
CorrResults = ProcessFOVCorrelationAndPlots(CaData, SpeedAll, StimAll, Suite2pTable, SessFileTable, ResultFolder);

%% --- 2. SLM/Cell Mapping ---
DistTh = 6;
[Suite2pTable, SLMtarget, SLMtargetTable, GroupTargetCell, TargetCellList, TargetCellListFunGroup] = ...
    ProcessFOVSLMTargetMapping(CaData, SLMPosInfo, SLMTestInfo, Suite2pTable, NeuronPos3D, DistTh);

%% --- 3. Group Comparison/Boxplots ---
Nlabel = {'DeltaF','Spks'};
ProcessFOVGroupComparison(CorrResults.rSpeed, CorrResults.rStim, TargetCellList, TargetCellListFunGroup, SLMPosInfo, Nlabel, ResultFolder);

%% --- 4. Cell Response/Heatmap Plots (Fill args as needed) ---
% You will likely need to compute PSTHparam, PVpower, etc. as in your original script
% Example placeholders below:
TimBinFrame = -PSTHparam.PreSLMCal:PSTHparam.PostSLMCal-1;
iscell = find(CaData.iscell(:,1)>0.99);
NData = {
    double(F2deltaFoF(CaData.F, CaData.Fneu, CaData.ops.fs));
    double(CaData.spks)
};
NData{1} = NData{1}(:, iscell)';   % deltaFoF: (frames, cells) -> (cells, frames)
NData{2} = NData{2}(iscell, :);    % spks: (cells, frames)

SLMInfoTable=Suite2pTable(Suite2pTable.suite2pInd>0,:);

PVpower=xmlPower2PVpower(SLMTestInfo.confSet.UncagingLaserPower);
PVpower=intersect(round(PVpower),SLMInfoTable.UncagingLaserPower);

% PVpower = 1; % TODO: Replace with proper calculation as in your script
% TargetResponse = {}; CellResponse = {}; CellSampleN = [];
% GroupTargetCell = {}; % TODO: Compute as needed
% For demonstration, call with placeholder/empty where needed
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

[AlignedtempNData,TargetInfoTable,CellResponse,statCellRes,TargetResponse,TargetCellResP,TargetCellResR,CellSampleN]=Aligned_FromSuite2p(NData{iData},TargetCellList,Suite2pTable,PVpower,PSTHparam);




if iData==1

FDR=0.1;
[h, crit_p, adj_p]=fdr_bh(TargetCellResP((~isnan(TargetCellResP))&CellSampleN>=TrialThNum),FDR,'pdep');
crit_p=min([crit_p 0.05]);


PowerTargetI=zeros(length(TargetCellList),1);
SuccTarget=PowerTargetI;
SuccAmp=PowerTargetI;

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
[h, crit_pAll, adj_p]=fdr_bh(pAll,FDR,'pdep');
crit_pAll=min([crit_pAll 0.05]);

end




Param.PlotType=3
Param.statisP=0;
Param.LegendShow=0;
Param.Legend=[]







AlignedNData{iData}=AlignedtempNData;




end

close all
save([ResultFolder 'Step1Meta.mat'],'Suite2pTable','CaData','CorrResults','SLMPosInfo','SLMTestInfo','SLMInfoTable','SessFileTable',...
    'NData','Neuronstat','NeuronPos3D','NeuronPos3DRaw','SuccAmp','SuccTarget','PowerTargetI',...
    'SuccAmp','TargetCellList','TargetCellListFunGroup','GroupTargetCell','AlignedSpeed','AlignedStim');


ProcessFOVCellResponsePlots(NData{1}, TargetCellList, iscell, Suite2pTable, PVpower, PSTHparam, GroupTargetCell, TargetResponse, CellResponse, CellSampleN, TimBinFrame, ResultFolder, Nlabel{1});
ProcessFOVCellResponsePlots(NData{2}, TargetCellList, iscell, Suite2pTable, PVpower, PSTHparam, GroupTargetCell, TargetResponse, CellResponse, CellSampleN, TimBinFrame, ResultFolder, Nlabel{2});

close all


end
















function [TargetResponse, statCellRes, TargetCellResP, TargetCellResR, CellSampleN] = SLMTargetResponseStats(NData, TargetCellList, Suite2pTable, PVpower, PSTHparam)
% Computes SLM target response statistics using ranksum/t-test for pre/post SLM on each target
% Inputs:
%   NData: (nCells x nFrames) neural data
%   TargetCellList: indices of target cells
%   Suite2pTable: trial/cycle metadata
%   PVpower: list of SLM/uncaging power settings
%   PSTHparam: parameters struct (PreSLMCal, PostSLMCal, pTh, TestMethod)

nTarget = length(TargetCellList);
nPower = length(PVpower);
TrialThNum = 3;  % Minimal number of trials per power for stat
TargetResponse = cell(nTarget, nPower);
statCellRes = struct([]);
TargetCellResP = NaN(nTarget, nPower);
TargetCellResR = NaN(nTarget, nPower);
CellSampleN = zeros(nTarget, nPower);

for iCell = 1:nTarget
    for iPower = 1:nPower
        % For each cell/power, get trial indices
        % Example: (replace with your method to select trial indices per power)
        trials = find(Suite2pTable.Point == iCell & Suite2pTable.UncagingLaserPower == PVpower(iPower));
        if length(trials) < TrialThNum
            continue; % skip if too few
        end
        responses = [];
        pvals = [];
        stats = [];
        for iTrial = 1:length(trials)
            trialInd = trials(iTrial);
            frameIdx = Suite2pTable.suite2pInd(trialInd);
            preIdx = frameIdx - PSTHparam.PreSLMCal:frameIdx-1;
            postIdx = frameIdx:frameIdx+PSTHparam.PostSLMCal-1;
            datPre = NData(TargetCellList(iCell), preIdx);
            datPost = NData(TargetCellList(iCell), postIdx);
            % Test (ranksum or ttest)
            if strcmpi(PSTHparam.TestMethod, 'ranksum')
                [p, ~, statsTmp] = ranksum(datPost, datPre);
            else
                [~, p, ~, statsTmp] = ttest2(datPost, datPre);
            end
            pvals = [pvals; p];
            stats = [stats; statsTmp];
            responses = [responses; mean(datPost) - mean(datPre)];
        end
        TargetResponse{iCell, iPower} = responses;
        statCellRes(iCell,iPower).p = pvals;
        statCellRes(iCell,iPower).stat = stats;
        TargetCellResP(iCell,iPower) = median(pvals);
        TargetCellResR(iCell,iPower) = mean(responses);
        CellSampleN(iCell,iPower) = length(responses);
    end
end

end













function SaveFOVResults(ResultFolder, varargin)
% Save relevant analysis results to Step1Meta.mat (or custom .mat)
% varargin: name-value pairs of variables to save

if mod(length(varargin), 2) ~= 0
    error('Arguments must be name-value pairs');
end

vars = struct();
for k = 1:2:length(varargin)
    name = varargin{k};
    value = varargin{k+1};
    vars.(name) = value;
end

save(fullfile(ResultFolder, 'Step1Meta.mat'), '-struct', 'vars');
end


