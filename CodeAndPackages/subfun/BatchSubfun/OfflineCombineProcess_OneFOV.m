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
        FileIDList(i)
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


TestGroup = [];
for iGroup = 1:length(Group)
    TestGroup = [TestGroup;zeros(length(Group(iGroup).Indices),1)+iGroup];
end

MeanImg=AmpNormalizeDim(double(CaData.PlaneMeanImg),3,[0.5 99.5]);
ImgClim=[0 1];

figure;
GroupColor=[255 51 153;91 20 212;121 247 111]/255;

H=MultiPlanes2DShow(permute(MeanImg,[2 1 3]), [], Pos3DFun, [], SLMPosInfo.confSetFinal.ETL+SLMPosInfo.confSetFinal.scan_Z(1), GroupColor(TestGroup,:), ImgClim);
bar=colorbar(H(3));
bar.Location='eastoutside';
bar.Position=[0.95,0.3,0.01,0.3];
bar.Label.String='F.';
bar.Ticks=[0 1];
papersizePX=[0 0 30 9];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));
print(gcf, [ResultFolder 'FinalPowertarget.svg'], '-dsvg', '-painters');
print(gcf, [ResultFolder 'FinalPowertarget.tif'], '-dtiffn', '-painters');          



%% --- 1. Correlation & Speed/Stim/Neural plots ---
CorrResults = ProcessFOVCorrelationAndPlots(CaData, SpeedAll, StimAll, Suite2pTable, SessFileTable, ResultFolder);


%% --- 2. SLM/Cell Mapping ---
DistTh = 10;
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

[DeltaF,DeltaF1]=F2deltaFoF(CaData.F, CaData.Fneu, CaData.ops.fs);

NData = {
    double(DeltaF);
    % double(DeltaF1);
    double(CaData.spks)
};
NData{1} = NData{1}(:, iscell)';   % deltaFoF: (frames, cells) -> (cells, frames)
% NData{2} = NData{2}(:, iscell)';
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

TrialThNum=2;

[AlignedtempNData,TargetInfoTable,CellResponse,statCellRes,TargetResponse,TargetCellResP,TargetCellResR,CellSampleN]=Aligned_FromSuite2p(NData{iData},TargetCellList,Suite2pTable,PVpower,PSTHparam);




if iData==1

FDR=0.1;
[h, crit_p, adj_p]=fdr_bh(TargetCellResP((~isnan(TargetCellResP))&CellSampleN>=TrialThNum),FDR,'pdep');
crit_p=0.05;
% crit_p=min([crit_p 0.05]);
% crit_p=0.1;

PowerTargetI=zeros(length(TargetCellList),1);
SuccTarget=PowerTargetI;
SuccAmp=PowerTargetI;
% SuccTargetPvalue=ones(length(TargetCellList),length(PVpower));
for iCell=1:length(TargetCellList)
    for iPower=1:length(PVpower)
        if TargetCellResP(iCell,iPower)<=crit_p&&TargetCellResR(iCell,iPower)>0&&CellSampleN(iCell,iPower)>=TrialThNum
           PowerTargetI(iCell)=iPower;
           SuccTarget(iCell)=1;
           SuccAmp(iCell)=TargetCellResR(iCell,iPower);
        end
        % if TargetCellResR(iCell,iPower)>0&&CellSampleN(iCell,iPower)>=TrialThNum
        %    SuccTargetPvalue(iCell,iPower)=TargetCellResP(iCell,iPower);
        % end

    end
end
sum(SuccTarget);

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
Param.Legend=[];



[~,~,ip]=intersect(TargetCellList,SLMtarget);
if length(ip)~=length(TargetCellList)
   disp(['Check TargetCellList and PointList, not match']);
   return 

end


AlignedNData{iData}=AlignedtempNData;

ProcessFOVCellResponsePlots(TargetCellList, ip ,iscell, PVpower, PSTHparam, TargetCellResP,TargetResponse, CellResponse, CellSampleN, TimBinFrame, ResultFolder, Nlabel{iData});



end

close all
% save([ResultFolder 'Step1Meta.mat'],'Suite2pTable','CaData','CorrResults','SLMPosInfo','SLMTestInfo','SLMInfoTable','SessFileTable',...
%     'NData','Neuronstat','NeuronPos3D','NeuronPos3DRaw','SuccAmp','SuccTarget','PowerTargetI',...
%     'SuccAmp','TargetCellList','TargetCellListFunGroup','GroupTargetCell','AlignedSpeed','AlignedStim');


    tempTargetList=TargetCellList;
    tempPowerI=PowerTargetI;
    tempPower=unique(SLMInfoTable.UncagingLaserPower);
    tempPower=tempPower(2:end);
    tempGroup=TargetCellListFunGroup;

    % targetListLabel={};
    % for i=1:length(TargetCellList)
    %     targetListLabel
    % 
    % end

        finalTargetList=[];
        finalTargetListGroup=[];
        for igroup=1:length(GroupTargetCell)
            finalTargetList=[finalTargetList;GroupTargetCell{igroup}(:)];
            finalTargetListGroup=[finalTargetListGroup;zeros(length(GroupTargetCell{igroup}),1)+igroup];
        end
        [finalTargetList,i1]=sort(finalTargetList);
        finalTargetListGroup=finalTargetListGroup(i1);

        i1=ismember(TargetCellList,finalTargetList);
        validI=find(tempPowerI>0&i1>0);
        tempTargetList=tempTargetList(validI);
        tempPowerI=tempPowerI(validI);
        tempPower=tempPower(tempPowerI);
        tempGroup=tempGroup(validI);

        % NeuronPos3D(tempTargetList,:);

    FinalActCellList=tempTargetList;
    FinalActPowerI=tempPowerI;
    FinalActCellFunGroup=tempGroup;
    FinalActPower=tempPower;

NeuronPos3Dum=NeuronPos3D;
NeuronPos3Dum(:,1:2)=NeuronPos3Dum(:,1:2)*SLMPosInfo.yaml.umPerlPixelX;

Pos3DFunum=Pos3DFun;
Pos3DFunum(:,1:2)=Pos3DFunum(:,1:2)*SLMPosInfo.yaml.umPerlPixelX;



CellDist=squareform(pdist(NeuronPos3Dum));
CellDistXY=squareform(pdist(NeuronPos3Dum(:,1:2)));


MinTargetDist=min(pdist2(NeuronPos3Dum,Pos3DFunum),[],2);
MinTargetDistXY=min(pdist2(NeuronPos3Dum(:,1:2),Pos3DFunum(:,1:2)),[],2);

for iGroup=1:length(Group)
    ActGroupPos3Dum{iGroup}=NeuronPos3Dum(FinalActCellList(FinalActCellFunGroup==iGroup),:);
    CellGroupDist(:,iGroup) = mean(pdist2(NeuronPos3Dum,NeuronPos3Dum(FinalActCellList(FinalActCellFunGroup==iGroup),:)),2);
    CellGroupDistXY(:,iGroup) = mean(pdist2(NeuronPos3Dum(:,1:2),NeuronPos3Dum(FinalActCellList(FinalActCellFunGroup==iGroup),1:2)),2);

end

   

save([ResultFolder 'Step1Meta.mat'],'Suite2pTable','CaData','CorrResults','SLMPosInfo','SLMTestInfo','SLMInfoTable','SessFileTable',...
    'NData','Neuronstat','NeuronPos3D','NeuronPos3Dum','NeuronPos3DRaw','SuccAmp','SuccTarget','PowerTargetI',...
    'SuccAmp','TargetCellList','TargetCellListFunGroup','GroupTargetCell','AlignedSpeed','AlignedStim',...
    'FinalActCellList','FinalActPowerI','FinalActCellFunGroup','FinalActPower','MinTargetDist','MinTargetDistXY','CellGroupDistXY','CellGroupDist','ActGroupPos3Dum');

ProcessFOVGroupComparison(CorrResults.rSpeed, CorrResults.rStim, FinalActCellList, FinalActCellFunGroup, SLMPosInfo, Nlabel, [ResultFolder 'Act']);


colorCell=repmat([1 1 0],length(tempTargetList),1);
colorCellIncluded=repmat([1 1 0],length(finalTargetList),1);

GroupLabel={'L','S','N'};
nGroup=length(GroupLabel);
% NodeColor=repmat([0.9 0.9 0.9],length(iscell),1);
for iGroup=1:length(Group)
    I1=find(tempGroup==iGroup);
    colorCell(I1,:)=repmat(GroupColor(iGroup,:),length(I1),1);

    I1=find(finalTargetListGroup==iGroup);
    colorCellIncluded(I1,:)=repmat(GroupColor(iGroup,:),length(I1),1);
end

finalTargetListLabel={};
for iCell=1:length(finalTargetList)
    finalTargetListLabel{iCell}=num2str(finalTargetList(iCell));
end

figure;
H=MultiPlanes2DShow(permute(MeanImg,[2 1 3]), cellBoundary(finalTargetList), NeuronPos3D(finalTargetList,:), finalTargetListLabel, SLMPosInfo.confSetFinal.ETL+SLMPosInfo.confSetFinal.scan_Z(1), colorCellIncluded, ImgClim);
bar=colorbar(H(3));
bar.Location='eastoutside';
bar.Position=[0.95,0.3,0.01,0.3];
bar.Label.String='F.';
bar.Ticks=[0 1];
papersizePX=[0 0 30 9];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));

print(gcf, [ResultFolder 'FinalUsedSLMtarget.svg'], '-dsvg', '-painters');
print(gcf, [ResultFolder 'FinalUsedSLMtarget.tif'], '-dtiffn', '-painters');          
% close all

tempTargetListLabel={};
for iCell=1:length(tempTargetList)
    tempTargetListLabel{iCell}=num2str(tempTargetList(iCell));
end

figure;
H=MultiPlanes2DShow(permute(MeanImg,[2 1 3]), cellBoundary(tempTargetList), NeuronPos3D(tempTargetList,:), tempTargetListLabel, SLMPosInfo.confSetFinal.ETL+SLMPosInfo.confSetFinal.scan_Z(1), colorCell, ImgClim);
bar=colorbar(H(3));
bar.Location='eastoutside';
bar.Position=[0.95,0.3,0.01,0.3];
bar.Label.String='F.';
bar.Ticks=[0 1];
papersizePX=[0 0 30 9];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));

print(gcf, [ResultFolder 'FinalActSLMtarget.svg'], '-dsvg', '-painters');
print(gcf, [ResultFolder 'FinalActSLMtarget.tif'], '-dtiffn', '-painters');          
close all

[AlignedtempNData,TargetInfoTable,CellResponse,statCellRes,TargetResponse,TargetCellResP,TargetCellResR,CellSampleN]=Aligned_FromSuite2p(NData{iData},TargetCellList,Suite2pTable,PVpower,PSTHparam);

% ProcessFOVCellResponsePlots(NData{2}, TargetCellList, iscell, Suite2pTable, PVpower, PSTHparam, GroupTargetCell, TargetResponse, CellResponse, CellSampleN, TimBinFrame, ResultFolder, Nlabel{2});

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


