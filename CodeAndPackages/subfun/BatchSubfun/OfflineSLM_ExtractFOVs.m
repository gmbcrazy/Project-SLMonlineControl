function [Output,NeuroTrace,BehTrace]=OfflineSLM_ExtractFOVs(FOV, Suite2pDataKeywords,suite2pFOVPathLocal,PSTHparam)



AlignedInfoTable=[];
CellResponseMeta={};
TargetResponseMeta={};
AlignedSpeedMeta=[];
AlignedStimMeta=[];
TableSessVec=[];
NeuroPos3DMeta=[];
NeuroPos3DumMeta=[];


rSpeedMeta=[];
rStimMeta=[];

MinTargetDistFOV={};
MinTargetDistXYFOV={};

TargetDistFOV={};
TargetDistXYFOV={};

ActGrouopPos3DumFOV={};


for iFOV=1:length(FOV)
suite2pFOVPathLocalTemp=suite2pFOVPathLocal{iFOV};
% GroupLabel={'L','S','N'};
% nGroup=length(GroupLabel);
% GroupColor=[255 51 153;91 20 212;121 247 111]/255;
% NodeColor=repmat([0.9 0.9 0.9],length(iscell),1);
% 
papersizePX=[0 0 33 9];
papersizePXX=[0 0 22 16];
TestStepFrame=PSTHparam.TestStepFrame;
% suite2pFOVPathLocalTemp=suite2pFOVPathLocal{iFOV};
resultPaths = findAllFoldersKeyWords(suite2pFOVPathLocalTemp, Suite2pDataKeywords);
metaDataPath=Get_ExpDataFolder(resultPaths{1},'Step1Basic',{'Step1Meta.mat','BehAll.mat','.png'});
load([metaDataPath 'Step1Meta.mat'])
load([metaDataPath 'BehAll.mat'])

PowerTargetIFOV{iFOV}=PowerTargetI;
PowerTestCellListFOV{iFOV}=TargetCellList;
PowerTestCellListFunGroupFOV{iFOV}=TargetCellListFunGroup;

confSet=SLMPosInfo.confSetFinal;
rSpeedMeta=cat(1,rSpeedMeta,CorrResults.rSpeed);
rStimMeta=cat(1,rStimMeta,CorrResults.rStim);

Zdepth = confSet.scan_Z + confSet.ETL;
suite2pPath = findAllFoldersKeyWords(suite2pFOVPathLocalTemp, Suite2pDataKeywords);
[~,Session,~]=fileparts(suite2pFOVPathLocalTemp(1:end-1));

confSet.save_path0=suite2pPath{1};
[NeuronPos3D,NeuronPos3DRaw,CaData,CaDataPlane,Neuronstat]=Extract_Suite2p(confSet);
% load('C:\Users\zhangl33\Projects\Project-SLMonlineControl\CodeAndPackages\subfun\Color\colorMapPN3.mat','colorMapPN1')
[~, ~, CellMedCenter, cellBoundary, ~] = Suite2pCellIDMapFromStat(CaData.statCell, [confSet.SLM_Pixels_Y confSet.SLM_Pixels_X]);
Cell3DPos=[CellMedCenter Zdepth(CaData.CellPlaneID)'];
CellDistMeta{iFOV}=squareform(pdist(Cell3DPos));

iData=PSTHparam.iData;
% PSTHparam.TestStepFrame=3;
PVpower=xmlPower2PVpower(SLMTestInfo.confSet.UncagingLaserPower);
PVpower=intersect(round(PVpower),SLMInfoTable.UncagingLaserPower);
PowerFOV{iFOV}=PVpower;


NeuroTrace{iFOV}=NData;
BehTrace(iFOV).Speed=SpeedAll;
BehTrace(iFOV).Stim=StimAll;
BehTrace(iFOV).AlignedSpeed=AlignedSpeed;
BehTrace(iFOV).AlignedStim=AlignedStim;



GroupTargetCellMeta{iFOV}=GroupTargetCell;


ActCellListFOV{iFOV}=FinalActCellList;
ActPowerIFOV{iFOV}=FinalActPowerI;
ActCellFunGroupFOV{iFOV}=FinalActCellFunGroup;
ActCellPowerFOV{iFOV}=FinalActPower;


MinTargetDistFOV{iFOV}=MinTargetDist;
MinTargetDistXYFOV{iFOV}=MinTargetDistXY;

TargetDistFOV{iFOV}=CellGroupDist;
TargetDistXYFOV{iFOV}=CellGroupDistXY;

ActGrouopPos3DumFOV{iFOV}=ActGroupPos3Dum;

% [AlignedtempNData,AlignedInfoTable,statCellRes,TargetResponse,TargetCellResP,TargetCellResR,CellSampleN]=Aligned_FromSuite2p(NData{iData},PowerTestCellList,iscell,Suite2pTable,PVpower,PSTHparam);
[AlignedNData{iFOV},AlignedInfoTableFOV{iFOV},CellResponseMeta{iFOV},statCellRes{iFOV},TargetResponseMeta{iFOV},TargetCellResPMeta{iFOV},TargetCellResR{iFOV},CellSampleN{iFOV}]=Aligned_FromSuite2p(NData{iData},TargetCellList,Suite2pTable,PVpower,PSTHparam);

AlignedInfoTable=[AlignedInfoTable;AlignedInfoTableFOV{iFOV}];
AlignedSpeedMeta=[AlignedSpeedMeta AlignedSpeed];
AlignedStimMeta=[AlignedStimMeta AlignedStim];




GroupTargetCellMeta{iFOV}=GroupTargetCell;
NeuroPos3DMeta=[NeuroPos3DMeta;[NeuronPos3D zeros(size(NeuronPos3D,1),1)+iFOV]];
NeuroPos3DumMeta=[NeuroPos3DumMeta;[NeuronPos3Dum zeros(size(NeuronPos3Dum,1),1)+iFOV]];

SuccAmpMeta{iFOV}=SuccAmp(:);
SuccTargetMeta{iFOV}=SuccTarget(:);
confSetMeta(iFOV)=confSet;

SLMTestInfoMeta(iFOV)=SLMTestInfo;
SLMPosInfoMeta(iFOV)=SLMPosInfo;

TableSessVec=[TableSessVec;zeros(size(AlignedInfoTableFOV{iFOV},1),1)+iFOV];

end

TableSessVec=table(TableSessVec,'VariableNames',{'iFOV'});

Output.AlignedNData=AlignedNData;
Output.AlignedInfoTableFOV=AlignedInfoTableFOV;
Output.AlignedInfoTable=[AlignedInfoTable TableSessVec];
Output.AlignedSpeedMeta=AlignedSpeedMeta;
Output.AlignedStimMeta=AlignedStimMeta;
Output.GroupTargetCellMeta=GroupTargetCellMeta;
Output.NeuroPos3DMeta=NeuroPos3DMeta;
Output.NeuroPos3DumMeta=NeuroPos3DumMeta;

Output.SuccAmpMeta=SuccAmpMeta;
Output.SuccTargetMeta=SuccTargetMeta;
Output.confSetMeta=confSetMeta;
Output.SLMTestInfoMeta=SLMTestInfoMeta;
Output.SLMPosInfoMeta=SLMPosInfoMeta;
Output.rSpeed=rSpeedMeta;
Output.rStim=rStimMeta;
Output.CellDistMeta=CellDistMeta;
Output.TargetCellResR=TargetCellResR;
Output.statCellRes=statCellRes;
Output.TargetResponseMeta=TargetResponseMeta;
Output.TargetCellResPMeta=TargetCellResPMeta;
Output.CellSampleN=CellSampleN;
Output.PowerTargetIFOV=PowerTargetIFOV;
Output.PowerTestCellListFOV=PowerTestCellListFOV;
Output.PowerTestCellListFunGroupFOV=PowerTestCellListFunGroupFOV;
Output.MinTargetDistFOV=MinTargetDistFOV;
Output.MinTargetDistXYFOV=MinTargetDistXYFOV;
Output.TargetDistFOV=TargetDistFOV;
Output.TargetDistXYFOV=TargetDistXYFOV;
Output.ActGrouopPos3DumFOV=ActGrouopPos3DumFOV;


GroupList=1:length(GroupTargetCell);

for iFun=1:length(GroupList)
    GroupTargetCell{iFun}=[];
    GroupTargetCellTemp=[];
    Cnum=0;
    for iFOV = 1:length(Output.GroupTargetCellMeta)
         GroupTargetCellTemp=[GroupTargetCellTemp;Output.GroupTargetCellMeta{iFOV}{iFun}(:)+Cnum];
         Cnum = Cnum+sum(abs(Output.NeuroPos3DMeta(:,4)-iFOV)<0.1);
    end
    GroupTargetCell{iFun}=GroupTargetCellTemp;
end

GroupTargetCellAll=[];
for iFun=1:length(GroupList)
    GroupTargetCellAll=[GroupTargetCellAll;[GroupTargetCell{iFun}(:) zeros(size(GroupTargetCell{iFun}(:)))+iFun]];
end
GroupTargetCellMeta=[GroupTargetCell {[]}];

Output.GroupTargetCellMerge=GroupTargetCell;
Output.GroupTargetCellAll=GroupTargetCellAll;



ActGroupTargetCellAll=[];
for iFOV = 1:length(ActCellListFOV)
    ActGroupTargetCellAll=[ActGroupTargetCellAll;[ActCellListFOV{iFOV}(:) ActCellFunGroupFOV{iFOV}(:)]];
end

for iFun=1:length(GroupList)
    ActGroupTargetCell{iFun}=[];
    % ActGroupTargetCellTemp=[];
    Cnum=0;
    for iFOV = 1:length(ActCellListFOV)
         % Ind=ActCellFunGroupFOV{iFOV}==iFun
         ActGroupTargetCell{iFun}=[ActGroupTargetCell{iFun}(:);ActCellListFOV{iFOV}(ActCellFunGroupFOV{iFOV}==iFun)+Cnum];
         Cnum = Cnum+sum(abs(Output.NeuroPos3DMeta(:,4)-iFOV)<0.1);
    end
    % GroupTargetCell{iFun}=GroupTargetCellTemp;
end


Output.ActGroupTargetCellMerge=ActGroupTargetCell;
Output.ActGroupTargetCellAll=ActGroupTargetCellAll;



Output.ActCellListFOV=ActCellListFOV;
Output.ActPowerIFOV=ActPowerIFOV;
Output.ActCellFunGroupFOV=ActCellFunGroupFOV;
Output.ActCellPowerFOV=ActCellPowerFOV;


