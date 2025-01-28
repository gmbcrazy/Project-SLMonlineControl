clear all
% TestFile='TSeries-04222024-0926-040'
WorkingFolder='E:\LuSLMOnlineTest\SL0838-Ai203\01242025\'
% load('C:\Users\zhangl33\Projects\Project-SLMonlineControl\subfun\Color\colorMapPN3.mat');
% load('C:\Users\User\Project-SLMonlineControl\subfun\Color\colorMapPN3.mat');
confSet = ReadYaml([WorkingFolder 'CurrentSLMsetting.yml']);

ProcessFolder=[WorkingFolder 'SingleP\' 'Top19SpeedStimEdgeExc\'];
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


umPerPixel=mean([yaml.umPerlPixelX yaml.umPerlPixelY]);
%param for ROI neighbourhood to determine wether there is SLM response.
ROIparam.TotalSLMPos3D=Pos3Dneed;    %%such that ROIparam.PointAll=1:size(Pos3Dneed,1)
ROIparam.PointAll=PointAll;
% ROIparam.PlaneZ=PlaneZ;
ROIparam.CellSize=20;                %%normal neuron diameter by um;        
ROIparam.threshold_percentage=0.3;   %%thereshold to define responsive fields SLM responsive heatmap: percentage*Peak rate
ROIparam.thNum=10;                   %%Minimal single responsive field by pixels
ROIparam.max_distance=ceil(ROIparam.CellSize/2/umPerPixel);  %% 1/3 diameter of a cell by pixel as maximal response region-SLM center distance
ROIparam.min_region_size=5;
ROIparam.PeakTh=200;
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


XMLparam.SwitchXMLPostMPFrame=8;                     %%<-----------------------------------------------------------MarkPoint switching occurs after 10 Repetitions of nplanes of Zseries.
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
OutTBLAll=[];               %%SLM trial information across all testing files;
PSTHall=[];
ROIall=[];
iCount=1;
minTrialN=SLMTestParam.TerminalTrialN;
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
RefFile=[WorkingFolder 'RefReg2\'];
% 
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
XMLparam.PointList=[1 4 7 19 20 22 1 4 7 19];
% XMLparam.Laser(1:10)=repmat(1.6,1,10);


[XMLTableTemp,FileGenerateInfoTemp]=PV_LinkPowerTest_MultiZseries(XMLparam,PowerTestPVPar);

% idRanges=[6;12];
% idRanges=[6 18;16 18];
idRanges=[FileGenerateInfoTemp.FileID;FileGenerateInfoTemp.FileID];   %Automatic update the new File ID to calculate ROIs

idRanges=[40;40];   %Automatic update the new File ID to calculate ROIs

step4_MultiZ_SubStep2_UpdatingROIandXML
close all



%% Save results and generate final MarkPoints and Functional groups
SLMTableOrigin=SLMTable;
% PostSLMTable;
FalsePositiveID=input('Mannual correction: index in SLMTable with false positive error: ');
SLMTable(FalsePositiveID,:)=nan;



SMLTablePowerPV=xmlPower2PVpower(SLMTable(:,2));
refPVpower=max(SMLTablePowerPV);
CellPerGroup=10;
[Group, FinalPos3D, FinalCellstat, FinalFunScore, confSetFinal] = SLMWeightsAssignToFunGroups(FunScore, CellPerGroup, Pos3Dneed, Cellstat, SLMIncludedIndFromIscell, SLMTable, NonTargets, refPVpower, confSet);
save([ProcessFolder 'SLMFunGroup.mat'],'Group','FinalPos3D','FinalCellstat','FinalFunScore','confSetFinal','SLMTableOrigin','SLMTable','ROIparam','SLMRes','sampleN','SLMTestParam','SLMIncludedIndFromIscell','FunScore','yaml','Cellstat');
XYZtoMarkPointFunGroup_MultiZ(ProcessFolder,FinalPos3D,Group,yaml,confSetFinal);







