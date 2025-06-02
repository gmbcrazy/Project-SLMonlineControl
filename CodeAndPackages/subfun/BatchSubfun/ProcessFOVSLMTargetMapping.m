function [Suite2pTable, SLMtarget, SLMtargetTable, GroupTargetCell, TargetCellList, TargetCellListFunGroup] = ProcessFOVSLMTargetMapping(CaData, SLMPosInfo, SLMTestInfo, Suite2pTable, NeuronPos3D, DistTh)
% Handles SLM target mapping, group assignment, and related table fields

% SLM Pos from SLMTestInfo
SLMPos3D = SLMTestInfo.Pos3Dneed;
SLMGroup = SLMTestInfo.FunScore(:,1);

[SLMtarget, SLMtargetCellDist] = SLMtargetMatchCell(SLMPos3D, NeuronPos3D, DistTh);
PointList1 = find(SLMtarget > 0);

SLMtargetTable = zeros(size(Suite2pTable,1),1);
SLMtargetTableGroup = zeros(size(Suite2pTable,1),1);
for iTarget = 1:length(SLMtarget)
    SLMtargetTable(Suite2pTable.Point == iTarget) = SLMtarget(iTarget);
    SLMtargetTableGroup(Suite2pTable.Point == iTarget) = SLMGroup(iTarget);
end
Suite2pTable.PointTargetCell = SLMtargetTable;
Suite2pTable.PointTargetCellGroup = SLMtargetTableGroup;

[SLMFinalInSLMtest,~] = SLMtargetMatchCell(SLMPosInfo.FinalPos3D, SLMPos3D, 0.1);

[TargetCellList, ia] = unique(SLMtargetTable(SLMtargetTable > 0));
temp = SLMtargetTableGroup(SLMtargetTable > 0);
TargetCellListFunGroup = temp(ia);

GroupTargetCell = cell(1, length(SLMPosInfo.Group));
for iFun = 1:length(SLMPosInfo.Group)
    [tempTarget,~] = SLMtargetMatchCell(SLMPosInfo.FinalPos3D(SLMPosInfo.Group(iFun).Indices,:), NeuronPos3D, DistTh);
    GroupTargetCell{iFun} = tempTarget(tempTarget > 0);
end






end
