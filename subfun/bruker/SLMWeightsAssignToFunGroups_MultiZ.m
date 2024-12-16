function [Group, FinalPos3D, FinalCellstat, FinalFunScore, confSetFinal] = SLMWeightsAssignToFunGroups_MultiZ(FunScore, CellPerGroup, Pos3Dneed, Cellstat, SLMIncludedIndFromIscell, SLMTable, NonTargets, refPVpower, confSet)
% ASSIGNCELLSTOGROUPS Assigns cells to functional groups based on criteria.
% Assign power weights based on SLM power settings in SLMTable

% Inputs:
%   FunScore - Functional scores of cells (matrix)
%   FunType - Types of functions for the cells (vector)
%   CellPerGroup - Number of cells per group (scalar)
%   Pos3Dneed - 3D positions of cells that are needed (matrix)
%   Cellstat - Cell statistics (cell array)
%   SLMIncludedIndFromIscell - SLM indices for included cells (vector)
%   SLMTable - SLM data table (matrix), col 1 is MarkPoints ID; col 2 is the xml PV power to stimulate the corresponding MarkPoint ID. 
%   NonTargets - Positions for non-target cells (matrix) % NonTargets
%   coordinates in case there is not enougch cells within functional group, adding non-cell targets in SLM group. 
%   refPVpower - Reference power for uncaging laser (scalar)
%   confSet - Configuration settings (structure)
% 
% Outputs:
%   FinalPos3D - Final 3D positions of cells (matrix)
%   FinalCellstat - Final cell statistics (cell array)
%   FinalFunScore - Final functional scores of cells (matrix)
%   confSetFinal - Updated configuration settings (structure)

FunType=unique(FunScore(:,1));
numFun=length(FunType);

% Initialize output variables


FinalPos3D = [];
FinalCellstat = {};
FinalFunScore = [];
numFun = length(FunType);


% Iterate through each functional group
for iGroup = 1:numFun
    % Determine indices for cells in the group
    Group(iGroup).Indices = [1:CellPerGroup] + (iGroup - 1) * CellPerGroup;
    
    % Assign functional group score
    FinalFunScore(Group(iGroup).Indices, 1) = iGroup;
    
    % Set default power weight for each cell in the group
    Group(iGroup).PowerWeight = ones(CellPerGroup, 1);

    % Find cells that match the functional group criteria
    I1 = find(FunScore(:, 1) == FunType(iGroup) & SLMTable(:, 2) > 0);

    % If there are more cells than needed
    if length(I1) >= CellPerGroup
        if iGroup < numFun
            % For functional groups, choose cells with the highest functional score
            [~, I2] = sort(FunScore(I1, iGroup + 1), 'descend');
            I1 = I1(I2(1:CellPerGroup));
        else
            % For non-functional groups, choose the first CellPerGroup cells
            I1 = I1(1:CellPerGroup);
        end
        AddedTargetN = 0;
        AddedPos = [];
        AddedCell = cell(1, 0);
        AddedFromNonTarget{iGroup} = [];
    else
        % If there are fewer cells than needed, add non-targets as fake cells
        AddedTargetN = CellPerGroup - length(I1);
        AddedPos = NonTargets(1:AddedTargetN, :);
        AddedCell = cell(1, AddedTargetN);
        AddedFromNonTarget{iGroup} = 1:AddedTargetN;
    end

    % Extract 3D positions and cell statistics for selected cells
    Pos3Dtemp = Pos3Dneed(I1, :);
    cellstattemp = Cellstat(SLMIncludedIndFromIscell(I1));

    % Assign power weights based on SLM power settings
    Group(iGroup).PowerWeight(1:length(I1)) = PVpower2PVweight(xmlPower2PVpower(SLMTable(I1, 2)), refPVpower);

    % Update functional scores for the group
    FinalFunScore(Group(iGroup).Indices(1:length(I1)), 2:numFun) = FunScore(I1, 2:numFun);
    FinalFunScore(Group(iGroup).Indices(length(I1) + 1:CellPerGroup), 2:numFun) = nan;

    % Append positions and cell statistics to final output
    FinalPos3D = [FinalPos3D; Pos3Dtemp; AddedPos];
    FinalCellstat = [FinalCellstat cellstattemp AddedCell];
end

% Update configuration settings
confSetFinal = confSet;
confSetFinal.UncagingLaserPower = PVpower2xmlPower(refPVpower);

end