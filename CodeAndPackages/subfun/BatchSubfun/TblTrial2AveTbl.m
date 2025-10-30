
function AveTable = TblTrial2AveTbl(SLMGroupTableTrial)
% AGGREGATOR for STRATEGY B
% Input:
%   SLMGroupTableTrial : tall table with at least columns
%       Response, PowerZero, Group, Session, Sensory, and
%       either AwakeState or Whisk (Awake surrogate), plus Cell.
%
% Output:
%   AveTable : per (Cell,Group,Sensory,AwakeState,Session) mean Response,
%              sham baseline ('ShamOpto'), GroupCount, ShamCount, and
%              ResponseNorm = Response - ShamOpto.

% --- pick state variable
if ismember('AwakeState', SLMGroupTableTrial.Properties.VariableNames)
    stateVar = 'AwakeState';
elseif ismember('Sensory', SLMGroupTableTrial.Properties.VariableNames)
    stateVar = 'Sensory';
else
    error('TblTrial2AveTbl: Missing AwakeState/Whisk column.');
end

% --- make sure Sensory is present (fix common capitalization)
if ~ismember('Sensory', SLMGroupTableTrial.Properties.VariableNames)
   error('TblTrial2AveTbl: Missing Sensory column.');
end

keysCommon  = {'Cell', stateVar, 'Sensory'};
keysNonSham = [keysCommon, {'Group'}];

% Split sham vs non-sham
TableNoSham = SLMGroupTableTrial(SLMGroupTableTrial.PowerZero==0, :);
TableSham   = SLMGroupTableTrial(SLMGroupTableTrial.PowerZero==1, :);

% Non-sham aggregation
if ~isempty(TableNoSham)
    Ave = groupsummary(TableNoSham, keysNonSham, "mean");
    % Ave.Properties.VariableNames{strcmp(Ave.Properties.VariableNames,'mean_Response')} = 'Response';
    Ave = groupsummaryBack2OldNames(TableNoSham, Ave, "mean");
% else
% 
% 
%     Ave = cell2table(cell(0,numel(keysNonSham)+2), 'VariableNames', [keysNonSham, {'GroupCount','Response'}]);
%     Ave.GroupCount = zeros(0,1);
end

% Sham aggregation (no 'Group' in keys)
if ~isempty(TableSham)
    Sham = groupsummary(TableSham, keysCommon, "mean");
    Sham.Properties.VariableNames{strcmp(Sham.Properties.VariableNames,'mean_Response')} = 'ShamOpto';
% else
%     Sham = cell2table(cell(0,numel(keysCommon)+2), 'VariableNames', [keysCommon, {'GroupCount','ShamOpto'}]);
end
if ismember('GroupCount', Sham.Properties.VariableNames)
    Sham.Properties.VariableNames{strcmp(Sham.Properties.VariableNames,'GroupCount')} = 'ShamCount';
else
    Sham.ShamCount = zeros(height(Sham),1);
end

% Left join so we never drop non-sham rows even if sham is missing
%Sham=Sham[:,{keysSham 'ShamOpto' 'ShamCount'}];

Sham = Sham(:, [keysCommon, {'ShamOpto', 'ShamCount'}]);
AveTable = innerjoin(Ave, Sham, 'Keys', keysCommon);

% Fill missing sham with 0 (Strategy A semantics when no sham exists)
if ~ismember('ShamOpto', AveTable.Properties.VariableNames)
    AveTable.ShamOpto = 0;
end
AveTable.ShamOpto = fillmissing(AveTable.ShamOpto, 'constant', 0);

% Normalized response (real group minus sham)
AveTable.ResponseNorm = AveTable.Response - AveTable.ShamOpto;
end

% Originalversion, Lu Zhang, Oct 21 2025

% function AveTable=TblTrial2AveTbl(SLMGroupTableTrial) 
% 
% 
% 
% TableNoShamTemp=SLMGroupTableTrial(SLMGroupTableTrial.PowerZero==0,:);
% TableShamTemp=SLMGroupTableTrial(SLMGroupTableTrial.PowerZero==1,:);
% GroupMethod='mean';
% AveTable = groupsummary(TableNoShamTemp, {'Cell', 'Group','Whisk','Session',}, GroupMethod);
% AveTable=groupsummaryBack2OldNames(TableNoShamTemp,AveTable,GroupMethod);
% ShamTable = groupsummary(TableShamTemp, {'Cell', 'Group', 'Whisk','Session'}, GroupMethod);
% ShamTable=groupsummaryBack2OldNames(TableShamTemp,ShamTable,GroupMethod);
% ShamTable.Properties.VariableNames{find(ismember(ShamTable.Properties.VariableNames,'Response')==1)}='ShamOpto';
% ShamTable=ShamTable(:,{'Cell','Whisk','Session','ShamOpto'});
% AveTable=innerjoin(AveTable,ShamTable,'Keys',{'Cell','Whisk','Session'});
% a=sum(AveTable2(AveTable2.Group==1,"ShamOpto")-AveTable2(AveTable2.Group==2,"ShamOpto"))
% AveTable.ResponseNorm=AveTable.Response-AveTable.ShamOpto;
% 
% 
% 
