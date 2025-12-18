function [Suite2pTable, SLMtarget, SLMtargetTable, GroupTargetCell, TargetCellList, TargetCellListFunGroup] = ProcessFOVSLMTarget_Online_Mapping(SLMPosInfo, SLMTestInfo, Suite2pTable, NeuronPos3D, Mapping)
% Handles SLM target mapping, group assignment, and related table fields

% SLM Pos from SLMTestInfo
SLMPos3D = SLMTestInfo.Pos3Dneed;
SLMGroup = SLMTestInfo.FunScore(:,1);
InitialCell=SLMTestInfo.SLMIncludedIndFromIscell;



% [SLMtarget, SLMtargetCellDist] = SLMtargetMatchCell(SLMPos3D, NeuronPos3D, DistTh);
% PointList1 = find(SLMtarget > 0);
SLMtarget=zeros(size(SLMPos3D,1),1)-1;

[CellFound,MapI]=ismember(InitialCell(:),Mapping(:,1));
% MapI=MapI(MapI>0);

for iCell=1:length(InitialCell)
    if CellFound(iCell)
       SLMtarget(iCell)=Mapping(MapI(iCell),2);
    end
end





SLMtargetTable = zeros(size(Suite2pTable,1),1);
SLMtargetTableGroup = zeros(size(Suite2pTable,1),1);
for iTarget = 1:length(SLMtarget)
    % [iTarget sum(Suite2pTable.Point == iTarget)];
    SLMtargetTable(Suite2pTable.Point == iTarget) = SLMtarget(iTarget);
    SLMtargetTableGroup(Suite2pTable.Point == iTarget) = SLMGroup(iTarget);
end
Suite2pTable.PointTargetCell = SLMtargetTable;
Suite2pTable.PointTargetCellGroup = SLMtargetTableGroup;

% [SLMFinalInSLMtest,~] = SLMtargetMatchCell(SLMPosInfo.FinalPos3D, SLMPos3D, 0.1);

[TargetCellList, ia] = unique(SLMtargetTable(SLMtargetTable > 0));
temp = SLMtargetTableGroup(SLMtargetTable > 0);
TargetCellListFunGroup = temp(ia);


%%Needs to be edited
GroupTargetCell = cell(1, length(SLMPosInfo.Group));

for iFun = 1:length(SLMPosInfo.Group)
    tempCellstat=SLMPosInfo.FinalCellstat(SLMPosInfo.Group(iFun).Indices);
    [map_idx,overlap] = mapCellStat1_To_CellStat2(SLMPosInfo.Cellstat, tempCellstat);

    map_idx=map_idx(map_idx > 0);
    [~,InLater]=ismember(map_idx,Mapping(:,1)); 
    InLater=InLater(InLater>0);
    GroupTargetCell{iFun} = Mapping(InLater,2);

end


% for iFun = 1:length(SLMPosInfo.Group)
%     [tempTarget,~] = SLMtargetMatchCell(SLMPosInfo.FinalPos3D(SLMPosInfo.Group(iFun).Indices,:), NeuronPos3D, DistTh);
%     GroupTargetCell{iFun} = tempTarget(tempTarget > 0);
% end
%%Needs to be edited





end
