
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
% SLMRes([20 27],1)=0;
SLMRes(19,1)=0;
sampleN(19,1)=4;
SLMRes([23 31 35],1)=0;

SLMRes([14 36],2)=0;

% SLMRes(28,2)=0;
% % SLMRes(13,2)=0;
% SLMRes(41,2)=0;
% SLMRes(8,2)=0;

% SLMRes([3 7 32 35],2)=0;

% sampleN([10],:)=10;
step4_MultiZ_SubStep1_PreTest  %% generate next points being test and update xml parameters.
