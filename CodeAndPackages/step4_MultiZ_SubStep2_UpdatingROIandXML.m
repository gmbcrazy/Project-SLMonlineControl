
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
SLMRes(3,1)=0;
SLMRes(15,1)=0;
SLMRes(41,1)=0;
SLMRes(4,2)=0;
% SLMRes(13,2)=0;
SLMRes(41,2)=0;
% SLMRes(8,2)=0;

% SLMRes([3 7 32 35],2)=0;

% sampleN([10],:)=10;
step4_MultiZ_SubStep1_PreTest  %% generate next points being test and update xml parameters.
