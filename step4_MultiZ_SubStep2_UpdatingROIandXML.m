
ROIall=[];
TBLall=[];

idRanges=[16;16];

[ROItemp, TBLtemp] = getSLMROI_BinMat(DataFolder, confSet, PSTHparam, Pos3Dneed, ROIparam, idRanges);
ROIall=cat(3,ROIall,ROItemp);
TBLall=[TBLall;TBLtemp];
SLMTrialInfo=TBLall;
SLMTrialInfo.Laser=PVpower2xmlPower(SLMTrialInfo.UncagingLaserPower);
FileIDrange=[];
[SLMRes,sampleN]=SLMResponse_ROIMultiZ(ROIall,SLMTrialInfo,ROIparam,minTrialN,SumDataFolder,FileIDrange);
% close all
step4_MultiZ_SubStep1_PreTest  %% generate next points being test and update xml parameters.

