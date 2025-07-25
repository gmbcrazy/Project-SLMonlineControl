
function [tbl, CovInfoTable]=OfflineSLM_FOVmeta2NeuroDeltaTrial(Input,ProcessPar,PSTHparam)

StructToVars(ProcessPar);
StructToVars(Input);

GroupMetaName = [GroupLabel {'FakeSLM'}];
GroupMetaColor = [GroupColor; PowerZeroColor];
NCell = size(NeuroPos3DMeta,1);
FunNum = length(GroupList) + 1;
TimBinFrame = -PSTHparam.PreSLMCal:PSTHparam.PostSLMCal-1;

TestStepFrame=PSTHparam.TestStepFrame;
        % clear statGroupRes


ResponseDelta=[];
DeltaSpeed=squeeze(mean(AlignedSpeedMeta(PSTHparam.PreSLMCal+[1:PSTHparam.TestStepFrame]+(PSTHparam.PostWinN-1)*TestStepFrame,:),1)-mean(AlignedSpeedMeta(1:PSTHparam.PreSLMCal,:),1));

CovInfoTable=[];
CovTarget=[];    


ResponseData=[];
TrialID=[];
FOVID=[];
CovTargetSpeedR=[];
CovTargetStimR=[];
CovPointSpeedR=[];
CovPointStimR=[];
CovCellSpeedR=[];
CovCellStimR=[];

CovSpeed=[];
CovCell=[];
for iFOV = 1:length(AlignedNData)
    AlignedDataTemp=squeeze(mean(AlignedNData{iFOV}(:,PSTHparam.PreSLMCal+[1:PSTHparam.TestStepFrame]+(PSTHparam.PostWinN-1)*TestStepFrame,:),2)-mean(AlignedNData{iFOV}(:,1:PSTHparam.PreSLMCal,:),2));
    
    
    tempSpeed=DeltaSpeed(AlignedInfoTable.iFOV==iFOV);
    
    NCellFOV=size(AlignedDataTemp,1);
    NeuroI=find(NeuroPos3DMeta(:,4)==iFOV);
    NTrial=size(AlignedInfoTableFOV{iFOV},1);
    
    CellTargetTempFOV=GroupTargetCellMeta{iFOV};
    TargetTemp=[];
    TargetSpeedR=[];
    TargetStimR=[];
    rSpeedTemp=rSpeed(NeuroI,1,1);
    rStimTemp=rStim(NeuroI,1,1);

    for iGroup=1:length(CellTargetTempFOV)
        TargetTemp=[TargetTemp;CellTargetTempFOV{iGroup}(:)];
        TargetSpeedR(iGroup)=mean(rSpeedTemp(CellTargetTempFOV{iGroup}));
        TargetStimR(iGroup)=mean(rStimTemp(CellTargetTempFOV{iGroup}));
    end
    
    CovTargetFOV=zeros(NCellFOV,1);
    CovTargetFOV(TargetTemp,1)=1;
    
    CovTargetFOV=repmat(CovTargetFOV,NTrial,1);
    CovTarget=[CovTarget;CovTargetFOV];



    CovSpeedFOV=repmat(tempSpeed(:)',NCellFOV,1);
    CovSpeedFOV=CovSpeedFOV(:);
    CovSpeed=[CovSpeed;CovSpeedFOV];


    CovCell=[CovCell;repmat(NeuroI,NTrial,1)];
    CovCellSpeedR=[CovCellSpeedR;repmat(squeeze(rSpeed(NeuroI,1,1)),NTrial,1)];
    CovCellStimR=[CovCellStimR;repmat(squeeze(rStim(NeuroI,1,1)),NTrial,1)];

    % CovTargetSpeedR=zeros(NCellFOV,NccccddTrial);
    % CovTargetStimR=zeros(NCellFOV,NTrial);
    % 
    % CovTargetSpeedR=repmat(TargetSpeedR(AlignedInfoTableFOV{iFOV}.Group));




    for iTrial=1:NTrial
        ResponseData=[ResponseData;AlignedDataTemp(:,iTrial)];
        TrialID=[TrialID;zeros(NCellFOV,1)+iTrial];
        FOVID=[FOVID;zeros(NCellFOV,1)+iFOV];

        CovInfoTable=[CovInfoTable;repmat(AlignedInfoTableFOV{iFOV}(iTrial,:),NCellFOV,1)];

        if ~isnan(AlignedInfoTableFOV{iFOV}.Group(iTrial))
           CovTargetSpeedR=[CovTargetSpeedR;zeros(NCellFOV,1)+TargetSpeedR(AlignedInfoTableFOV{iFOV}.Group(iTrial))];
           CovTargetStimR=[CovTargetStimR;zeros(NCellFOV,1)+TargetStimR(AlignedInfoTableFOV{iFOV}.Group(iTrial))];
           
        else
           CovTargetSpeedR=[CovTargetSpeedR;zeros(NCellFOV,1)+NaN];
           CovTargetStimR=[CovTargetStimR;zeros(NCellFOV,1)+NaN];

        end

        if AlignedInfoTableFOV{iFOV}.PointTargetCell(iTrial)>0
           CovPointSpeedR=[CovPointSpeedR;zeros(NCellFOV,1)+rSpeedTemp(AlignedInfoTableFOV{iFOV}.PointTargetCell(iTrial))];
           CovPointStimR=[CovPointStimR;zeros(NCellFOV,1)+rStimTemp(AlignedInfoTableFOV{iFOV}.PointTargetCell(iTrial))];
           
        else
           CovPointSpeedR=[CovPointSpeedR;zeros(NCellFOV,1)+NaN];
           CovPointStimR=[CovPointStimR;zeros(NCellFOV,1)+NaN];

        end

    end

end

tbl = table(ResponseData(:), CovCellSpeedR(:), CovCellStimR(:), ...
    CovTargetSpeedR(:), CovTargetStimR(:), CovSpeed(:),CovCell(:),CovTarget(:),CovPointSpeedR,CovPointStimR,TrialID,FOVID,...
    'VariableNames', {'Response', 'SpeedR','StimR','TargetSpeedR','TargetStimR','Speed','Cell','TargetCell','PointSpeedR','PointStimR','TrialID','Session'});



