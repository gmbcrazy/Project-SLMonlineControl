%% Run Initiall round of SLM test; MultiZseriesVersion

% XMLparam.RoundID = permute()
XMLparam.RoundID=randperm(XMLparam.TotalRounds,1);

XMLparam.PointList=PointsTest(iPP);
PSTHparam.TargetPos=Pos3Dneed(XMLparam.Point,:);
PSTHparam.CellStat=CaData.statCell{SLMIncludedIndFromIscell(XMLparam.Point)};


[SLMTrialInfoTemp SLMPowerFile(end+1)]=PV_LinkPowerTest_MultiZseries(XMLparam,PVparam);



for iPP=1:length(PointsTest)
    tic
XMLparam.Point=PointsTest(iPP);
PSTHparam.TargetPos=Pos3Dneed(XMLparam.Point,:);
PSTHparam.CellStat=CaData.statCell{SLMIncludedIndFromIscell(XMLparam.Point)};
[SLMTrialInfo(end+1,:) SLMTrialMap(:,:,:,end+1)]=PV_LinkExcuteXML(XMLparam,PVparam,confSet,PSTHparam);
% [SLMTrialInfo SLMTrialMap]=PV_LinkExcuteXML(XMLparam,PVparam,confSet,PSTHparam);

    toc

    if size(SLMTrialInfo,1)==size(SLMTrialMap,4)
       disp(['SLMTrialInfo matches SLMTrialMap'])
    elseif size(SLMTrialInfo,1)==size(SLMTrialMap,4)-1
           tempCheck=abs(SLMTrialMap(:,:,:,1));
           if sum(tempCheck(:))==0
              disp('1 zero redudnant data is found in SLMTrialMap, removed correctly');
               SLMTrialMap(:,:,:,1)=[];
           else
              disp('Warning: 1 non-zero redudnant data is found in SLMTrialMap, check it');
           end
    else
   
       disp(['Warning: SLMTrialInfo and SLMTrialMap doe not match']);

    end    
end