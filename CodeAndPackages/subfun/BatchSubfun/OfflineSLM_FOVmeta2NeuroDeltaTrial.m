function [tbl, CovInfoTable] = OfflineSLM_FOVmeta2NeuroDeltaTrial(Input, ProcessPar, PSTHparam)
% STRATEGY  - TRIAL TABLE
% Build a per-trial, per-cell table of response deltas and covariates.
%
% Required fields are the same as in Strategy A, plus:
%   ProcessPar.OffTargetPixel
%   Input.ActCellListFOV{fov}, Input.ActCellFunGroupFOV{fov},
%   Input.MinTargetDistFOV{fov}, Input.rSpeed, Input.rStim
%
% Returns:
%   tbl           : tall table, one row per (cell x trial)
%   CovInfoTable  : event metadata replicated per (cell x trial) row
%
% NOTE: Window parity with Strategy A when PSTHparam.PostWinN == 1.

StructToVars(ProcessPar);
StructToVars(Input);

% Time windows (dynamic post window by PostWinN)
preIdx  = (PSTHparam.PreSLMCal - PSTHparam.PreTestFrame + 1) : PSTHparam.PreSLMCal;
shift   = (PSTHparam.PostWinN - 1) * PSTHparam.TestStepFrame;
postIdx = PSTHparam.PreSLMCal + (1:PSTHparam.TestStepFrame) + shift;

% Collectors
ResponseData = [];   TrialID = [];  FOVID = [];
CovTarget = [];      CovNonTarget = [];
CovCell = [];        CovCellFOV = [];
CovCellSpeedR = [];  CovCellStimR = [];
CovTargetSpeedR = []; CovTargetStimR = [];
CovPointSpeedR  = []; CovPointStimR  = [];
CovSpeed = [];
CovInfoTable = [];

CovCellMinDistXY=[];
CovCellMinDist=[];

% CovCellGroupDistXY=[];
% CovCellGroupDist=[];


nFOV = numel(AlignedNData);

for iFOV = 1:nFOV
    % FOV-scoped indices and tables
    iFOVEvent = find(AlignedInfoTable.iFOV == iFOV);
    T = AlignedInfoTableFOV{iFOV};
    NTrial = height(T);

    % FOV-scoped speed (align order!)
    SpeedFOV = AlignedSpeedMeta(:, iFOVEvent);
    preSp  = mean(SpeedFOV(preIdx,:), 1, 'omitnan');
    postSp = mean(SpeedFOV(postIdx,:), 1, 'omitnan');
    speedDelta = postSp - preSp; % 1 x trials

    % FOV-scoped neural responses
    X = AlignedNData{iFOV}; % [cells x time x trials]
    preAvg  = squeeze(mean(X(:, preIdx,  :), 2, 'omitnan')); % cells x trials
    postAvg = squeeze(mean(X(:, postIdx, :), 2, 'omitnan')); % cells x trials
    if isvector(preAvg), preAvg = preAvg(:); postAvg = postAvg(:); end
    RespDelta = postAvg - preAvg;                                % cells x trials

    NCellFOV = size(RespDelta,1);
    NeuroI   = find(NeuroPos3DMeta(:,4) == iFOV);
    rSpeedTemp = squeeze(rSpeed(NeuroI,1,1));
    rStimTemp  = squeeze(rStim(NeuroI,1,1));

    % Per-cell covariates replicated over trials
    CovCell        = [CovCell;        repmat(NeuroI, NTrial, 1)];
    CovCellFOV     = [CovCellFOV;     repmat((1:numel(NeuroI))', NTrial, 1)];
    CovCellSpeedR  = [CovCellSpeedR;  repmat(rSpeedTemp, NTrial, 1)];
    CovCellStimR   = [CovCellStimR;   repmat(rStimTemp,  NTrial, 1)];

    CovCellMinDistXY  = [CovCellMinDistXY;  repmat(MinTargetDistXYFOV{iFOV}, NTrial, 1)];
    CovCellMinDist   = [CovCellMinDist;   repmat(MinTargetDistFOV{iFOV},  NTrial, 1)];

    % CovCellGroupDistXY=[CovCellGroupDistXY;repmat(MinTargetDistXYFOV{iFOV}, NTrial, 1)];
    % CovCellGroupDist=[CovCellGroupDist;repmat(MinTargetDistXYFOV{iFOV}, NTrial, 1)];


    % Trial-level speed replicated over cells
    CovSpeedFOV = repmat(speedDelta(:)', NCellFOV, 1);
    CovSpeed    = [CovSpeed; CovSpeedFOV(:)];





    % Target/non-target flags per cell, replicated over trials
    TargetCells = ActCellListFOV{iFOV};
    CovTargetFOV = zeros(NCellFOV,1);
    if ~isempty(TargetCells)
        CovTargetFOV(TargetCells) = ActCellFunGroupFOV{iFOV};
    end
    CovTarget    = [CovTarget;    repmat(CovTargetFOV, NTrial, 1)];

    nonTargetInd = MinTargetDistFOV{iFOV} > PSTHparam.OffTargetUm;
    CovNonTargetFOV = zeros(NCellFOV,1);  CovNonTargetFOV(nonTargetInd) = 1;
    CovNonTarget = [CovNonTarget; repmat(CovNonTargetFOV, NTrial, 1)];

    % Group-level target R (by functional group among target cells)
    nFunGroups = max([ActCellFunGroupFOV{iFOV}(:); 0]);
    TargetSpeedR = nan(1, nFunGroups);
    TargetStimR  = nan(1, nFunGroups);
    for gi = 1:nFunGroups
        sel = ActCellListFOV{iFOV}(ActCellFunGroupFOV{iFOV}==gi);
        TargetSpeedR(gi) = mean(rSpeedTemp(sel), 'omitnan');
        TargetStimR(gi)  = mean(rStimTemp(sel),  'omitnan');
    end

    % Append row-wise per trial
    for k = 1:NTrial
        ResponseData = [ResponseData; RespDelta(:,k)];
        TrialID      = [TrialID;     repmat(k, NCellFOV, 1)];
        FOVID        = [FOVID;       repmat(iFOV, NCellFOV, 1)];

        CovInfoTable = [CovInfoTable; repmat(T(k,:), NCellFOV, 1)];

        if ~isnan(T.Group(k)) && T.Group(k)>=1 && T.Group(k)<=numel(TargetSpeedR)
            CovTargetSpeedR = [CovTargetSpeedR; repmat(TargetSpeedR(T.Group(k)), NCellFOV, 1)];
            CovTargetStimR  = [CovTargetStimR;  repmat(TargetStimR(T.Group(k)),  NCellFOV, 1)];
        else
            CovTargetSpeedR = [CovTargetSpeedR; nan(NCellFOV,1)];
            CovTargetStimR  = [CovTargetStimR;  nan(NCellFOV,1)];
        end

        if  T.PointTargetCell(k) > 0
            pt = T.PointTargetCell(k);
            CovPointSpeedR = [CovPointSpeedR; repmat(rSpeedTemp(pt), NCellFOV, 1)];
            CovPointStimR  = [CovPointStimR;  repmat(rStimTemp(pt),  NCellFOV, 1)];
        else
            CovPointSpeedR = [CovPointSpeedR; nan(NCellFOV,1)];
            CovPointStimR  = [CovPointStimR;  nan(NCellFOV,1)];
        end
    end
end

tbl = table( ...
    ResponseData(:), CovCellSpeedR(:), CovCellStimR(:), ...
    CovTargetSpeedR(:), CovTargetStimR(:), CovSpeed(:), ...
    CovCell(:), CovCellFOV(:), CovTarget(:), CovNonTarget(:), ...
    CovPointSpeedR(:), CovPointStimR(:), TrialID(:), FOVID(:), ...
    'VariableNames', {'Response','SpeedR','SensoryR', ...
                      'TargetSpeedR','TargetSensoryR','Speed', ...
                      'Cell','CellFOV','TargetCell','NonTargetCell', ...
                      'PointSpeedR','PointSensoryR','TrialID','Session'} );


%     'VariableNames', {'Response', 'SpeedR','SensoryR','TargetSpeedR','TargetSensoryR','Speed','Cell','CellFOV','TargetCell','NonTargetCell','PointSpeedR','PointSensoryR','TrialID','Session'});

end


%% Originalversion, Lu Zhang, Oct 21 2025


% function [tbl, CovInfoTable]=OfflineSLM_FOVmeta2NeuroDeltaTrial(Input,ProcessPar,PSTHparam) 
% 
% StructToVars(ProcessPar);
% StructToVars(Input);
% 
% GroupMetaName = [GroupLabel {'FakeSLM'}];
% GroupMetaColor = [GroupColor; PowerZeroColor];
% NCell = size(NeuroPos3DMeta,1);
% FunNum = length(GroupList) + 1;
% TimBinFrame = -PSTHparam.PreSLMCal:PSTHparam.PostSLMCal-1;
% 
% TestStepFrame=PSTHparam.TestStepFrame;
%         % clear statGroupRes
% 
% 
% ResponseDelta=[];
% DeltaSpeed=squeeze(mean(AlignedSpeedMeta(PSTHparam.PreSLMCal+[1:PSTHparam.TestStepFrame]+(PSTHparam.PostWinN-1)*TestStepFrame,:),1)-mean(AlignedSpeedMeta((PSTHparam.PreSLMCal-PSTHparam.PreTestFrame+1):PSTHparam.PreSLMCal,:),1));
% 
% CovInfoTable=[];
% CovTarget=[];    
% CovNonTarget=[];    
% 
% 
% ResponseData=[];
% TrialID=[];
% FOVID=[];
% CovTargetSpeedR=[];
% CovTargetStimR=[];
% CovPointSpeedR=[];
% CovPointStimR=[];
% CovCellSpeedR=[];
% CovCellStimR=[];
% 
% CovSpeed=[];
% CovCell=[];
% CovCellFOV=[];
% 
% 
% 
% for iFOV = 1:length(AlignedNData)
%     AlignedDataTemp=squeeze(mean(AlignedNData{iFOV}(:,PSTHparam.PreSLMCal+[1:PSTHparam.TestStepFrame]+(PSTHparam.PostWinN-1)*TestStepFrame,:),2)-mean(AlignedNData{iFOV}(:,(PSTHparam.PreSLMCal-PSTHparam.PreTestFrame+1):PSTHparam.PreSLMCal,:),2));
% 
% 
%     tempSpeed=DeltaSpeed(AlignedInfoTable.iFOV==iFOV);
% 
%     NCellFOV=size(AlignedDataTemp,1);
%     NeuroI=find(NeuroPos3DMeta(:,4)==iFOV);
%     NTrial=size(AlignedInfoTableFOV{iFOV},1);
% 
%     CellTargetTempFOV=GroupTargetCellMeta{iFOV};
% 
% 
%     TargetTemp=[];
%     TargetSpeedR=[];
%     TargetStimR=[];
%     rSpeedTemp=rSpeed(NeuroI,1,1);
%     rStimTemp=rStim(NeuroI,1,1);
% 
%     % for iGroup=1:length(CellTargetTempFOV)
%     %     TargetTemp=[TargetTemp;CellTargetTempFOV{iGroup}(:)];
%     %     TargetSpeedR(iGroup)=mean(rSpeedTemp(CellTargetTempFOV{iGroup}));
%     %     TargetStimR(iGroup)=mean(rStimTemp(CellTargetTempFOV{iGroup}));
%     % end
% 
%     TargetTemp=ActCellListFOV{iFOV};
%     for iGroup=1:length(CellTargetTempFOV)
%         % TargetTemp=[TargetTemp;CellTargetTempFOV{iGroup}(:)];
%         TargetSpeedR(iGroup)=mean(rSpeedTemp(TargetTemp(ActCellFunGroupFOV{iFOV}==iGroup)));
%         TargetStimR(iGroup)=mean(rStimTemp(TargetTemp(ActCellFunGroupFOV{iFOV}==iGroup)));
%     end
% 
% 
% 
%     CovTargetFOV=zeros(NCellFOV,1);
%     CovTargetFOV(TargetTemp,1)=ActCellFunGroupFOV{iFOV};    
%     CovTargetFOV=repmat(CovTargetFOV,NTrial,1);
%     CovTarget=[CovTarget;CovTargetFOV];
% 
% 
%     NonTargetInd=MinTargetDistFOV{iFOV}>PSTHparam.OffTargetPixel;
%     CovNonTargetFOV=zeros(NCellFOV,1);
%     CovNonTargetFOV(NonTargetInd,1)=1;    
%     CovNonTargetFOV=repmat(CovNonTargetFOV,NTrial,1);
%     CovNonTarget=[CovNonTarget;CovNonTargetFOV];
% 
% 
% 
% 
%     CovSpeedFOV=repmat(tempSpeed(:)',NCellFOV,1);
%     CovSpeedFOV=CovSpeedFOV(:);
%     CovSpeed=[CovSpeed;CovSpeedFOV];
% 
% 
%     CovCell=[CovCell;repmat(NeuroI,NTrial,1)];
%     CovCellFOV=[CovCellFOV;repmat([1:length(NeuroI)]',NTrial,1)];
% 
% 
%     CovCellSpeedR=[CovCellSpeedR;repmat(squeeze(rSpeed(NeuroI,1,1)),NTrial,1)];
%     CovCellStimR=[CovCellStimR;repmat(squeeze(rStim(NeuroI,1,1)),NTrial,1)];
% 
%     % CovTargetSpeedR=zeros(NCellFOV,NccccddTrial);
%     % CovTargetStimR=zeros(NCellFOV,NTrial);
%     % 
%     % CovTargetSpeedR=repmat(TargetSpeedR(AlignedInfoTableFOV{iFOV}.Group));
% 
% 
% 
% 
%     for iTrial=1:NTrial
%         ResponseData=[ResponseData;AlignedDataTemp(:,iTrial)];
%         TrialID=[TrialID;zeros(NCellFOV,1)+iTrial];
%         FOVID=[FOVID;zeros(NCellFOV,1)+iFOV];
% 
%         CovInfoTable=[CovInfoTable;repmat(AlignedInfoTableFOV{iFOV}(iTrial,:),NCellFOV,1)];
% 
%         if ~isnan(AlignedInfoTableFOV{iFOV}.Group(iTrial))
%            CovTargetSpeedR=[CovTargetSpeedR;zeros(NCellFOV,1)+TargetSpeedR(AlignedInfoTableFOV{iFOV}.Group(iTrial))];
%            CovTargetStimR=[CovTargetStimR;zeros(NCellFOV,1)+TargetStimR(AlignedInfoTableFOV{iFOV}.Group(iTrial))];
% 
%         else
%            CovTargetSpeedR=[CovTargetSpeedR;zeros(NCellFOV,1)+NaN];
%            CovTargetStimR=[CovTargetStimR;zeros(NCellFOV,1)+NaN];
% 
%         end
% 
%         if AlignedInfoTableFOV{iFOV}.PointTargetCell(iTrial)>0
%            CovPointSpeedR=[CovPointSpeedR;zeros(NCellFOV,1)+rSpeedTemp(AlignedInfoTableFOV{iFOV}.PointTargetCell(iTrial))];
%            CovPointStimR=[CovPointStimR;zeros(NCellFOV,1)+rStimTemp(AlignedInfoTableFOV{iFOV}.PointTargetCell(iTrial))];
% 
%         else
%            CovPointSpeedR=[CovPointSpeedR;zeros(NCellFOV,1)+NaN];
%            CovPointStimR=[CovPointStimR;zeros(NCellFOV,1)+NaN];
% 
%         end
% 
%     end
% 
% end
% 
% tbl = table(ResponseData(:), CovCellSpeedR(:), CovCellStimR(:), ...
%     CovTargetSpeedR(:), CovTargetStimR(:), CovSpeed(:),CovCell(:),CovCellFOV(:),CovTarget(:),CovNonTarget(:),CovPointSpeedR,CovPointStimR,TrialID,FOVID,...
%     'VariableNames', {'Response', 'SpeedR','SensoryR','TargetSpeedR','TargetSensoryR','Speed','Cell','CellFOV','TargetCell','NonTargetCell','PointSpeedR','PointSensoryR','TrialID','Session'});



