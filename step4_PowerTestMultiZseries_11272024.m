clear all
% TestFile='TSeries-04222024-0926-040'

load('C:\Users\User\Project-SLMonlineControl\subfun\Color\colorMapPN3.mat');

ConfigFolder='C:\Users\User\Project-SLMonlineControl\config\';

SLMsettingFile='SLMsetting.yml';
confSet = ReadYaml([ConfigFolder '\' SLMsettingFile]);

nPlane=length(confSet.ETL)
% DataFolder='F:\LuSLMOnlineTest\04222024\Data\'
ProcessFolder='E:\LuSLMOnlineTest\SL0777-Ai203\12092024\SingleP\Top8SpeedStimEdgeExc\';
DataFolder=[ProcessFolder 'Data\'];
mkdir(DataFolder);
DataLogFolder=[ProcessFolder 'DataLog\'];
SumDataFolder=[ProcessFolder 'DataSum\'];
mkdir(SumDataFolder);

load([ProcessFolder 'SLMIncludedIndFromIscell.mat'],'Pos3Dneed','yaml');
AllTestPoints3D=Pos3Dneed;
PointAll=1:size(AllTestPoints3D,1);

%%
SLMTestParam.TerminalTrialN=4;    %<-------------------------------------------------------------------------------Edit, Trials # to define SLM responsive cells
SLMTestParam.ExcludeTrialN=2;     %<-------------------------------------------------------------------------------Edit, Trials # to define Non-SLM responsive cells
SLMTestParam.AllLaserPower=confSet.UncagingLaserPower;% Noted that, laser test levels is dependent on ROIparam.LaserPower, not SLMTestParam.AllLaserPower

PowerTestPVPar.nPlane=nPlane;            
PowerTestPVPar.ZRepetition=31;                       %%<----------------------------------------------------- -----NumOfRepeition in each Zseries of Tseries in PV
PowerTestPVPar.Ziteration=11;                        %%<-----------------------------------------------------------NumOfZseries in Tseries in PV
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
ROIparam.max_distance=ceil(ROIparam.CellSize/3/umPerPixel);  %% 1/3 diameter of a cell by pixel as maximal response region-SLM center distance
ROIparam.min_region_size=5;
ROIparam.PeakTh=200;
ROIparam.min_merged_region_size=40;  %%Minimal total size of responsive fields by pixels
ROIparam.contourMethod='boundaries';      %%Method to detect ROI boader of responsive fields
ROIparam.NeighbourHfWidthPixel=20;   %%PixFromMedCenter: Number of pixels from the median center of each ROI to get the ROI neighborhood.
ROIparam.umPerPixel=mean([yaml.umPerlPixelX yaml.umPerlPixelY]);  
ROIparam.Colormap=colorMapPN1;                  
ROIparam.LaserPower=confSet.UncagingLaserPower(1:3);   
ROIparam.Clim=[-400;400];

PSTHparam.PreSLMCal=15;        %<----------------------------------------------------------------------------------Edit,Frame # before SLM to calculate baseline map
PSTHparam.PostSLMCal=3;        %<----------------------------------------------------------------------------------Edit,Frame # after SLM to calculate responsive map
PSTHparam.YLim=[-50 600];       % %<-----------------------------------------For Suite2p based ROI signal only method
PSTHparam.pTh=0.05;             % %<-----------------------------------------For Suite2p based ROI signal only method
PSTHparam.TestMethod='ranksum'; % %<-----------------------------------------For Suite2p based ROI signal only method


XMLparam.SwitchXMLPostMPFrame=6;                     %%<-----------------------------------------------------------MarkPoint switching occurs after 10 Repetitions of nplanes of Zseries.
XMLparam.ProcessFolder=ProcessFolder;
XMLparam.AllPointList=1:size(AllTestPoints3D,1);

tempAllList=repmat(XMLparam.AllPointList(:),1,SLMTestParam.TerminalTrialN); %Each Point needs to be test for at most SLMTestParam.ExcludeTrialN times
tempAllVector=tempAllList(:)
XMLparam.PointList=tempAllVector(1:PowerTestPVPar.Ziteration-1) %%<-----------nP, NumOfTestedPoints, nP + 1 = NumOfZseries 
XMLparam.TotalRounds= confSet.Repetition;


%%
SLMRes=zeros(length(XMLparam.AllPointList),length(ROIparam.LaserPower));
sampleN=SLMRes;
OutTBLAll=[];               %%SLM trial information across all testing files;
PSTHall=[];
ROIall=[];
iCount=1;
minTrialN=4;
SLMTrialInfo=[];                  %Inital response information, automatically updated after each single trial test
SLMTrialMap=[];                   %Inital response map, automatically updated after each single trial test
clear SLMTable;
SLMTable(:,1)=round(1:size(Pos3Dneed,1));
SLMTable(:,2)=NaN;

%% 
step4_MultiZ_SubStep1_PreTest  %% generate next points being test and update xml parameters.



[XMLTable{iCount},FileGenerateInfo(iCount)]=PV_LinkPowerTest_MultiZseries(XMLparam,PowerTestPVPar);


step4_MultiZ_SubStep2_PostTest  %% generate next points being test and update xml parameters.






%%
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
XYZtoMarkPointFunGroup(ProcessFolder,FinalPos3D,Group,yaml,confSetFinal);







