function [newtbl2,sumT1]=PowerTestTrialTBL2PowerTestTBL(Input,SLMPointTableTrial,activate_Pth)

r1Total=[];
r2Total=[];
% for iFOV=1:3
rgroup=[];
ngroup=[];
tblFOVclean=[];   %%refers to remove trials with power test when the power is not strong enough to activate one target cell.


SessionAll=unique(SLMPointTableTrial.Session);
iData=1;
length(Input(iData))

    for iFOV=1:length(Input(1).AlignedNData)
        % iFOV=SessionAll(jFOV);
        close all
        tblFOV=SLMPointTableTrial(SLMPointTableTrial.Session==iFOV,:);
        tempNeuroPos3D=Input(1).NeuroPos3DumMeta(Input(1).NeuroPos3DumMeta(:,4)==iFOV,1:3);
        % tempNeuroPos3D(:,1:3)=tempNeuroPos3D(:,1:3)*FOVscalePerPixel;

        tempTargetList=Input(iData).PowerTestCellListFOV{iFOV};
        tempPowerI=Input(iData).PowerTargetIFOV{iFOV};
        tempPower=unique(tblFOV.UncagingLaserPower);
        tempGroup=Input(iData).PowerTestCellListFunGroupFOV{iFOV};

        % 
        finalTargetList=Input(iData).ActCellListFOV{iFOV};   %%finalTargetList is actually the ActCellList, sucessufully detected and activated cells
        % 
        % finalTargetList=sort(finalTargetList);
        % i1=ismember(tempTargetList,finalTargetList);

        finalTargetList=sort(finalTargetList);
        i1=ismember(tempTargetList,finalTargetList);


        % tempPowerI=Input(iData).ActPowerIFOV{iFOV};
        % tempGroup=Input(iData).ActCellFunGroupFOV{iFOV};
        % tempTargetList=Input(iData).ActCellListFOV{iFOV};

    % tempTargetList=tempTargetList(i1);
    % tempPowerI=tempPowerI(i1);
    % tempGroup=tempGroup(i1);
    %Cells pass the tested might be more than the final included cell, this
    %part only inlucde the cell later used in SLM group stimuli

        %Only significant activated cell were considered
        validI=find(tempPowerI>0&i1>0);
        tempTargetList=tempTargetList(validI);
        tempPowerI=tempPowerI(validI);
        tempPower=tempPower(tempPowerI);
        tempGroup=tempGroup(validI);
        %Only significant activated cell were considered




    for iCell=1:length(validI)
    
        ActivateI=Input(iData).statCellRes{iFOV}(validI(iCell),tempPowerI(iCell)).p<activate_Pth;
        TrialID=unique(tblFOV.TrialID(tblFOV.PointTargetCell==tempTargetList(iCell)&tblFOV.UncagingLaserPower==tempPower(iCell)));
        tempData=[];
        for itrial=1:length(TrialID)
            tempData(:,itrial)=tblFOV.Response(tblFOV.TrialID==TrialID(itrial));
        end
    
        rgroup(end+1)=tempGroup(iCell);
    
        r1=corr(tempData);
        r1Total(end+1)=mean(r1(triu(true(size(r1)),1)));
        if sum(ActivateI)>5
           r2=corr(tempData(ActivateI,:));
           r2Total(end+1)=mean(r2(triu(true(size(r2)),1)));
        else
           r2Total(end+1)=NaN;
        end
    
        temptbl=tblFOV(ismember(tblFOV.TrialID,TrialID),:);
        temptbl.ActivateI=repmat(ActivateI(:),length(TrialID),1);
        temptbl.VecR=repmat(r1Total(end),size(temptbl,1),1);
        temptbl.ActVecR=repmat(r2Total(end),size(temptbl,1),1);
        

        % Coortbl=[repmat(tempNeuroPos3D,length(TrialID),1) repmat(tempNeuroPos3D(tempTargetList(iCell),:),size(temptbl,1),1)];

        Diff3Dtemp=tempNeuroPos3D-tempNeuroPos3D(tempTargetList(iCell),:);
        Dist3Dtemp=sqrt(sum(Diff3Dtemp.^2,2));
        Dist2Dtemp=sqrt(sum(Diff3Dtemp(:,1:2).^2,2));
        CoLayer=Diff3Dtemp(:,3)==0;

        spaMetrictemp= table(Dist3Dtemp(:),Dist2Dtemp(:),CoLayer(:),'VariableNames',{'Dist','DistXY','CoLayer'});
        
        Coortbl=[tempNeuroPos3D repmat(tempNeuroPos3D(tempTargetList(iCell),:),size(tempNeuroPos3D,1),1)];
        Coortbl=array2table(Coortbl,'VariableNames',{'X','Y','Z','PointX','PointY','PointZ'});

        temptbl=[temptbl repmat([spaMetrictemp Coortbl],length(TrialID),1)];
    
        tblFOVclean=[tblFOVclean;temptbl];
    
        
    end
    
    end




GroupMethod='mean';

% illegalSession=[5 6];
% tblFOVclean(ismember(tblFOVclean.Session,illegalSession),:)=[];



% Define positive and negative activation based on Response direction
% ActivateIPos: activated with positive response (>0)
% ActivateINeg: activated with negative response (<0)
tblFOVclean.ActivateIPos=tblFOVclean.ActivateI.*tblFOVclean.Response>0;
tblFOVclean.ActivateINeg=tblFOVclean.ActivateI.*tblFOVclean.Response<0;

% Multiply activation mask with SpeedR to get weighted activation values
% (e.g., Speed response only for active cells)
tblFOVclean.ActSpeedR=tblFOVclean.ActivateI.*tblFOVclean.SpeedR;
tblFOVclean.ActPosSpeedR=tblFOVclean.ActivateIPos.*tblFOVclean.SpeedR;
tblFOVclean.ActNegSpeedR=tblFOVclean.ActivateINeg.*tblFOVclean.SpeedR;

% Multiply activation mask with SensoryR to get weighted activation values
% (e.g., Sensory response only for active cells)
tblFOVclean.ActSensoryR=tblFOVclean.ActivateI.*tblFOVclean.SensoryR;
tblFOVclean.ActPosSensoryR=tblFOVclean.ActivateIPos.*tblFOVclean.SensoryR;
tblFOVclean.ActNegSensoryR=tblFOVclean.ActivateINeg.*tblFOVclean.SensoryR;

GroupMethod='mean';
% Keep only non-target cells for further analysis
tblFOVclean=tblFOVclean(tblFOVclean.NonTargetCell==1,:);



% Aggregate table by Session, PointTargetCell, and Cell
% Compute mean values for each group
tblFOVcleanavg = groupsummary(tblFOVclean, {'Session', 'PointTargetCell','Cell'}, GroupMethod);
newtbl2=groupsummaryBack2OldNames(tblFOVclean,tblFOVcleanavg,GroupMethod);

% Recalculate positive and negative activation flags after averaging
newtbl2.ActivateIPos=newtbl2.ActivateI.*newtbl2.Response>0;
newtbl2.ActivateINeg=newtbl2.ActivateI.*newtbl2.Response<0;



% Further average by Session and PointTargetCell (cell-level aggregation)
newtbl3 = groupsummary(newtbl2, {'Session', 'PointTargetCell'}, GroupMethod);
newtbl3=groupsummaryBack2OldNames(newtbl2,newtbl3,GroupMethod);

% Define grouping indices for custom splitapply functions
[G, sessionVals, TargetcellVals] = findgroups(newtbl2.Session, newtbl2.PointTargetCell);
% [G, sessionVals, TargetcellVals] = findgroups(newtbl2.Session, newtbl2.PointTargetCell);

% Define custom functions:
% myfun: computes weighted mean of y, weighted by x (sum(x.*y)/sum(x))
% myfun2: computes contrast between positive and negative activation
myfun = @(x,y) sum(x.*y) / sum(x);
myfun2 = @(x,y,z) sum((x.*z)-sum(y.*z)) / sum(x+y);

% Apply for each pair
AIabs_SpeedR  = splitapply(myfun, newtbl2.ActivateI,    newtbl2.SpeedR, G);
AIP_SpeedR = splitapply(myfun, newtbl2.ActivateIPos, newtbl2.SpeedR, G);
AIN_SpeedR = splitapply(myfun, newtbl2.ActivateINeg, newtbl2.SpeedR, G);
AI_SpeedR =  splitapply(myfun2, newtbl2.ActivateIPos,newtbl2.ActivateINeg, newtbl2.SpeedR, G);


AIabs_SensoryR   = splitapply(myfun, newtbl2.ActivateI,    newtbl2.SensoryR,  G);
AIP_SensoryR  = splitapply(myfun, newtbl2.ActivateIPos, newtbl2.SensoryR,  G);
AIN_SensoryR  = splitapply(myfun, newtbl2.ActivateINeg, newtbl2.SensoryR,  G);
AI_SensoryR =  splitapply(myfun2, newtbl2.ActivateIPos,newtbl2.ActivateINeg, newtbl2.SensoryR, G);
% Build result table using the grouping labels that match
newtbl4 = table( ...
    sessionVals, TargetcellVals, ...
    AIabs_SpeedR,AI_SpeedR, AIP_SpeedR, AIN_SpeedR, ...
    AIabs_SensoryR,AI_SensoryR, AIP_SensoryR, AIN_SensoryR, ...
    'VariableNames', {'Session','PointTargetCell', ...
                      'AabsSpeedR','ASpeedR','APSpeedR','ANSpeedR', ...
                      'AabsSensoryR','ASensoryR','APSensoryR','ANSensoryR'});

sumT1 = join(newtbl3, newtbl4, 'Keys', {'Session','PointTargetCell'});
sumT1(:,{'TargetSpeedR','TargetSensoryR'})=[];