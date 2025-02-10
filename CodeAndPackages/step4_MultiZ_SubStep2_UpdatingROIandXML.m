
% ROIall=[];
% TBLall=[];
% % 
% idRanges=[3;16];

[ROItemp, TBLtemp,PSTHtemp] = getSLMROI_BinMat(DataFolder, confSet, PSTHparam, Pos3Dneed, ROIparam, idRanges);
ROIall=cat(3,ROIall,ROItemp);
TBLall=[TBLall;TBLtemp];
SLMTrialInfo=TBLall;
SLMTrialInfo.Laser=PVpower2xmlPower(SLMTrialInfo.UncagingLaserPower);
FileIDrange=[];
[SLMRes,sampleN]=SLMResponse_ROIMultiZ(ROIall,SLMTrialInfo,ROIparam,SLMTestParam.TerminalTrialN,SumDataFolder,FileIDrange);
% close all
step4_MultiZ_SubStep1_PreTest  %% generate next points being test and update xml parameters.
