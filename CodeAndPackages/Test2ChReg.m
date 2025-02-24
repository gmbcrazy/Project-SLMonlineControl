
clear all
ProcessFolder='F:\LuSLMOnlineTest\SL0541-Emx1Ai96\09162024\';
Ch1Folder='F:\LuSLMOnlineTest\SL0541-Emx1Ai96\09162024\Ch1\TSeries-09162024-0937-081\';

% SessinInfo=[];
% SLMTable=SLMtargetAnalysis(ProcessFolder, Param,GlobalSavePath,SessinInfo, TargetZ);


ConfigFolder=ProcessFolder;
% ConfigFolder='C:\Users\zhangl33\Projects\Project-SLMonlineControl\config\';
ConfigFile='SLMsettingG7.yml';

[~,~,~,CaData,CaDataPlane,stat,yaml,confSet]=ROIToXYZ(ConfigFolder,ConfigFile);
Ch2Img=CaData.PlaneMeanImg;
Ch2Img=permute(Ch2Img,[2,1,3]);
nPlane=3;
Ind=[1:300];
Ch1Img=MeanFrameIndTifFolder(Ch1Folder,nPlane,Ind);

 Zdepth=confSet.ETL+confSet.scan_Z(1);
 colorCell=[0 1 0];
 ImgClim=[0 400]
 figure;
H1=MultiPlanes2DShow(Ch1Img, [], [], [], Zdepth, colorCell, ImgClim)


 figure;
H2=MultiPlanes2DShow(Ch2Img, [], [], [], Zdepth, colorCell, ImgClim)

