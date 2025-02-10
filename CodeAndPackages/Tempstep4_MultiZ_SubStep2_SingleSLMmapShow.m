
% ROIall=[];
% TBLall=[];
% % 
% idRanges=[3;16];

[ROItemp, TBLtemp,~,PSTHtempPlane] = getSLMROI_BinMat(DataFolder, confSet, PSTHparam, Pos3Dneed, ROIparam, idRanges);
ROIall=cat(3,ROIall,ROItemp);
TBLall=[TBLall;TBLtemp];
SLMTrialInfo=TBLall;
SLMTrialInfo.Laser=PVpower2xmlPower(SLMTrialInfo.UncagingLaserPower);
FileIDrange=[];
[SLMRes,sampleN]=SLMResponse_ROIMultiZ(ROIall,SLMTrialInfo,ROIparam,SLMTestParam.TerminalTrialN,SumDataFolder,FileIDrange);
% close all
step4_MultiZ_SubStep1_PreTest  %% generate next points being test and update xml parameters.

Point=19
Ind=TBLtemp.Point==Point&TBLtemp.UncagingLaserPower==440;

meanPlane=mean(PSTHtempPlane(:,:,Ind),3);
Zdepth=unique(Pos3Dneed(:,3));
[~,PlaneI]=ismember(Pos3Dneed(Point,3),Zdepth);
ImgClim=[-200 200];

figure;
MultPlaneIs2DShow1Plane(meanPlane, [], Pos3Dneed(Point,:), [], unique(Pos3Dneed(:,3)), PlaneI, [0 1 0], ImgClim)
set(gca,'clim',ROIparam.Clim,'xtick',[],'ytick',[]);
colormap(ROIparam.Colormap);

figure;
MultPlaneIs2DShow1Plane(meanPlane, [], Pos3Dneed(Point,:), [], unique(Pos3Dneed(:,3)), PlaneI, [0 1 0], ImgClim)
set(gca,'clim',ROIparam.Clim,'xtick',[],'ytick',[]);
colormap(ROIparam.Colormap);
