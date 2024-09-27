%% Run Initiall round of SLM test
for iPP=1:length(PointsTest)
XMLparam.Point=PointsTest(iPP);
PSTHparam.TargetPos=Pos3DNeed(XMLparam.Point,:);
PSTHparam.CellStat=CaData.statCell{SLMIncludedIndFromIscell(XMLparam.Point)};
[SLMTrialInfo(end+1,:) SLMTrialMap(:,:,:,end+1)]=PV_LinkExcuteXML(XMLparam,PVparam,confSet,PSTHparam);

[size(SLMTrialInfo,1) size(SLMTrialMap,4)]

end