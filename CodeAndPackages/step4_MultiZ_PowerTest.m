clear all
% TestFile='TSeries-04222024-0926-040'
ConfigFolder='C:\Users\User\Project-SLMonlineControl\config\';


WorkingFolder='E:\LuSLMOnlineTest\L00121\10062025\'
% load('C:\Users\zhangl33\Projects\Project-SLMonlineControl\subfun\Color\colorMapPN3.mat');
load('C:\Users\User\Project-SLMonlineControl\subfun\Color\colorMapPN3.mat');
confSet = ReadYaml([WorkingFolder 'CurrentSLMsetting.yml']);

% ProcessFolder=[WorkingFolder 'SingleP\' 'Top49SpeedStimEdgeExc\'];

ProcessFolder = Get_ExpDataFolder(WorkingFolder, 'SpeedStimEdgeExc', {'AllIncluded','.gpl','.xml'})

% % ConfigFolder='C:\Users\User\Project-SLMonlineControl\config\';
% % 
% % SLMsettingFile='SLMsetting.yml';
% % confSet = ReadYaml([ConfigFolder '\' SLMsettingFile]);

step4_MultiZ_SubStep0_LoadData
%%
SLMTestParam.TerminalTrialN=4;    %<-------------------------------------------------------------------------------Edit, Trials # to define SLM responsive cells
SLMTestParam.ExcludeTrialN=2;     %<-------------------------------------------------------------------------------Edit, Trials # to define Non-SLM responsive cells
SLMTestParam.AllLaserPower=confSet.UncagingLaserPower;% Noted that, laser test levels is dependent on ROIparam.LaserPower, not SLMTestParam.AllLaserPower
PowerTestPVPar.nPlane=nPlane;            
PowerTestPVPar.ZRepetition=confSet.ZRepetition;      %%NumOfRepeition in each Zseries of Tseries in PV
PowerTestPVPar.Ziteration=confSet.Ziteration;                        %%NumOfZseries in Tseries in PV
PowerTestPVPar.InterMPRepetition=repmat(PowerTestPVPar.ZRepetition,1,PowerTestPVPar.Ziteration);
frameRepetition=PowerTestPVPar.ZRepetition*PowerTestPVPar.Ziteration;
PowerTestPVPar.maxFrame=nPlane*frameRepetition;

PowerTestTseries=dir([ConfigFolder 'PowerTest' num2str(confSet.Ziteration) 'Z' num2str(confSet.ZRepetition) '.env']);
LoadTSeriestoBruker(PowerTestTseries);


umPerPixel=mean([yaml.umPerlPixelX yaml.umPerlPixelY]);
%param for ROI neighbourhood to determine wether there is SLM response.
ROIparam.TotalSLMPos3D=Pos3Dneed;    %%such that ROIparam.PointAll=1:size(Pos3Dneed,1)
ROIparam.PointAll=PointAll;
% ROIparam.PlaneZ=PlaneZ;
ROIparam.CellSize=20;                %%normal neuron diameter by um;        
ROIparam.threshold_percentage=0.25;   %%thereshold to define responsive fields SLM responsive heatmap: percentage*Peak rate
ROIparam.thNum=10;                   %%Minimal single responsive field by pixels
ROIparam.max_distance=ceil(ROIparam.CellSize/2/umPerPixel);  %% 1/3 diameter of a cell by pixel as maximal response region-SLM center distance
ROIparam.min_region_size=5;
ROIparam.PeakTh=180;
ROIparam.min_merged_region_size=50;  %%Minimal total size of responsive fields by pixels
ROIparam.contourMethod='boundaries';      %%Method to detect ROI boader of responsive fields
ROIparam.NeighbourHfWidthPixel=20;   %%PixFromMedCenter: Number of pixels from the median center of each ROI to get the ROI neighborhood.
ROIparam.umPerPixel=mean([yaml.umPerlPixelX yaml.umPerlPixelY]);  
ROIparam.Colormap=colorMapPN1;                  
ROIparam.LaserPower=confSet.UncagingLaserPower;   
ROIparam.Clim=[-400;400];
ROIparam.PointsTest=ROIparam.PointAll;

PSTHparam.PreSLMCal=10;        %<----------------------------------------------------------------------------------Edit,Frame # before SLM to calculate baseline map
PSTHparam.PostSLMCal=3;        %<----------------------------------------------------------------------------------Edit,Frame # after SLM to calculate responsive map
PSTHparam.YLim=[-50 600];       % %<-----------------------------------------For Suite2p based ROI signal only method
PSTHparam.pTh=0.05;             % %<-----------------------------------------For Suite2p based ROI signal only method
PSTHparam.TestMethod='ranksum'; % %<-----------------------------------------For Suite2p based ROI signal only method
PSTHparam.MPFrameJump=2; % %<-----------------------------------------For Suite2p based ROI signal only method


XMLparam.SwitchXMLPostMPFrame=6;                     %%<-----------------------------------------------------------MarkPoint switching occurs after 10 Repetitions of nplanes of Zseries.
XMLparam.ProcessFolder=ProcessFolder;
XMLparam.AllPointList=1:size(AllTestPoints3D,1);

tempAllList=repmat(XMLparam.AllPointList(:),1,SLMTestParam.TerminalTrialN); %Each Point needs to be test for at most SLMTestParam.ExcludeTrialN times
tempAllVector=tempAllList(:)
XMLparam.PointList=tempAllVector(1:PowerTestPVPar.Ziteration-1) %%<-----------nP, NumOfTestedPoints, nP + 1 = NumOfZseries 
XMLparam.TotalRounds= confSet.Repetition;
XMLparam.Laser=confSet.UncagingLaserPower(1);

%%
SLMRes=zeros(length(XMLparam.AllPointList),length(ROIparam.LaserPower));
sampleN=SLMRes;

% Zdepth=confSet.ETL+confSet.scan_Z;
% PointPlane=ismember(Pos3Dneed(:,3),PlaneZ);
% SLMRes(PointPlane==3,1)=0;                   %No need to test low power for deepest plane, usually do not work
% sampleN(PointPlane==3,1)=100;                %No need to test low power for deepest plane, usually do not work



OutTBLAll=[];               %%SLM trial information across all testing files;
PSTHall=[];
ROIall=[];
iCount=1;
SLMTrialInfo=[];                  %Inital response information, automatically updated after each single trial test
SLMTrialMap=[];                   %Inital response map, automatically updated after each single trial test
clear SLMTable;
SLMTable(:,1)=round(1:size(Pos3Dneed,1));
SLMTable(:,2)=NaN;

%% 
numGPUs=0;      %%Do not use GPU, assume in general the aquisition PC has no GPU. 
FileType=2;   %Choose a specific bin file as reference for motion correction
% ProcessFolder='E:\LuSLMOnlineTest\SL0777-Ai203\12122024\'
% RefFile=[DataFolder 'TSeries-12132024-1247-023.bin'];
% [RegOps, RegImg] = LoadRegRefFile(RefFile, FileType, numGPUs, [512,512,3,30]);
% 
FileType=0;   %Choose a pre-recorded multi-tif files for motion correction
RefFile=[];
RefFile=[WorkingFolder 'RegRef1Ch1Ch2\'];
[RegOps, RegImg] = LoadRegRefFile(RefFile, FileType,numGPUs);
RegImg=RegImg(:,:,4:6);

RefFile=[WorkingFolder 'RegRef2Ch2\'];
[RegOps, RegImg] = LoadRegRefFile(RefFile, FileType,numGPUs);

% MultiMatrix3DHeatmap(RegImg)
%%
% FileType=1;   %Choose suite2p folder, using ops.meanImg for motion correction
% RefFile=[WorkingFolder 'suite2p\'];
% [RegOps, RegImg] = LoadRegRefFile(RefFile, FileType,numGPUs);

XMLparam.DoRegistration=1;
XMLparam.RegRefOps=RegOps;
XMLparam.RegRefImg=RegImg;  
XMLparam.RegRefSource=RefFile;  

disp('Referrence Img for Motion correction updated')


ROIall=[];
TBLall=[];
step4_MultiZ_SubStep1_PreTest  %% generate next points being test and update xml parameters.


%%
%% Mannual control when necessary
% XMLparam.PointList=[7 9 7 16 17 7 21 22 17 30];
% XMLparam.Laser(1:10)=[1.35 1.35 1.35 1.35 1.55 1.35 1.6 1.55 1.55 1.6];
% XMLparam.RoundID=randperm(XMLparam.TotalRounds,1);
XMLparam.PointList=[ 25    10    38    26    11     9    25    32    24    26];
% % % % % % % % % % XMLparam.PointList(10)=[29]
% XMLparam.Laser(1:10)=repmat(1.4,1,10);
% % 
pause(10)

[XMLTableTemp,FileGenerateInfoTemp]=PV_LinkPowerTest_MultiZseries(XMLparam,PowerTestPVPar);


% idRanges=[6;12];
% idRanges=[4;27];
idRanges=[FileGenerateInfoTemp.FileID;FileGenerateInfoTemp.FileID];   %Automatic update the new File ID to calculate ROIs
step4_MultiZ_SubStep2_UpdatingROIandXML
close all



%% Save results and generate final MarkPoints and Functional groups
SLMTableOrigin=SLMTable;
% PostSLMTable;
FalsePositiveID=input('Mannual correction: index in SLMTable with false positive error: ');
SLMTable(FalsePositiveID,2)=nan;

% SLMTable([37],2)=1.5;
% SLMTable([26 27 28 34],2)=1.5;


SMLTablePowerPV=xmlPower2PVpower(SLMTable(:,2));
refPVpower=max(SMLTablePowerPV);
CellPerGroup=10;
[Group, FinalPos3D, FinalCellstat, FinalFunScore, confSetFinal] = SLMWeightsAssignToFunGroups(FunScore, CellPerGroup, Pos3Dneed, Cellstat, SLMIncludedIndFromIscell, SLMTable, NonTargets, refPVpower, confSet);
save([ProcessFolder 'SLMFunGroup.mat'],'Group','FinalPos3D','FinalCellstat','FinalFunScore','confSetFinal','SLMTableOrigin','SLMTable','ROIparam','SLMRes','sampleN','SLMTestParam','SLMIncludedIndFromIscell','FunScore','yaml','Cellstat');
XYZtoMarkPointFunGroup_MultiZ(ProcessFolder,FinalPos3D,Group,yaml,confSetFinal);


%%

GroupLabel={'L','S','N'};
% nGroup=length(SLMPosInfo.Group);
GroupColor=[255 51 153;91 20 212;121 247 111]/255;
        PlotParam.RowPlot=1;
        PlotParam.RowColNum=1;
        PlotParam.RowColID=1;
        PlotParam.EdgeParam=[0.06 0.1 0.06 0.06 0.06 0.06];
        PlotParam.CellCenterWith=1.5;
        PlotParam.CellBoundaryWidth=0.5;
        PlotParam.PlotCenter=1;

[~,~,~,CaData,~,~,~,~]=ROIToXYZ(confSetFinal.save_path0,'CurrentSLMsetting.yml');
[~, ~, ~, FincalCellBoundary] = Suite2pCellIDMapFromStat(FinalCellstat, [confSetFinal.SLM_Pixels_Y confSetFinal.SLM_Pixels_X]);



colorCell=repmat([1 1 0],length(FinalCellstat),1);
for iGroup=1:length(Group)
    I1=find(FinalFunScore(:,1)==iGroup);
    colorCell(I1,:)=repmat(GroupColor(iGroup,:),length(I1),1);
end

MeanImg=AmpNormalizeDim(double(CaData.PlaneMeanImg),3,[0.5 99.5]);
ImgClim=[0 1];
H=MultiPlanes2DShow(permute(MeanImg,[2 1 3]), FincalCellBoundary, FinalPos3D, [], confSetFinal.ETL+confSetFinal.scan_Z(1), colorCell, ImgClim,PlotParam);
bar=colorbar(H(3));
bar.Location='eastoutside';
bar.Position=[0.95,0.3,0.01,0.3];
bar.Label.String='F.';
bar.Ticks=[0 1];
papersizePX=[0 0 30 9];
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'PaperPosition',papersizePX,'PaperSize',papersizePX(3:4));

print(gcf, [ProcessFolder 'FinalSLMtarget.svg'], '-dsvg', '-painters');
print(gcf, [ProcessFolder 'FinalSLMtarget.tif'], '-dtiffn', '-painters');          
close all


figure;
plot(0,0);set(gca,'xlim',[0 10],'ylim',[0 10],'xcolor','w','ycolor','w')
hold on;
text(0,6,['Save the xml into SavedMarkPoints in Bruker!'],'color','r','fontsize',13);



%% Copy all power test .xml and .gpl files to a new folder 'PowerTest'
PowerTestFolder=[ProcessFolder 'PowerTest\'];
mkdir(PowerTestFolder);

PowerTestList=SearchFile(ProcessFolder,'GPoint*',0);
for i=1:length(PowerTestList)
    copyfile([ProcessFolder PowerTestList(i).name], PowerTestFolder);
end

for i=1:length(PowerTestList)
    delete([ProcessFolder '\*GPoint*']);
end


