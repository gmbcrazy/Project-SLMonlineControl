%% Load Data
clear all
% ProcessFolder='F:\LuSLMOnlineTest\04222024\SingleP\30PixelFromEdgeExc\';
load('C:\Users\User\Project-SLMonlineControl\subfun\Color\colorMapPN3.mat');
ConfigFolder='C:\Users\User\Project-SLMonlineControl\config\';
ConfigFile='SLMsetting.yml';%<----------------------------------------------------------------------------------Edit, configuration file
[~,~,~,CaData,CaDataPlane,stat,yaml,confSet]=ROIToXYZ(ConfigFolder);
umPerPixel=mean([yaml.umPerlPixelX yaml.umPerlPixelY]);
ProcessFolder='E:\LuSLMOnlineTest\SL0777-Ai203\10302014\SingleP\Top13SpeedStimEdgeExc\';%<----------------------Edit, Data folder

step4_SubStep1_LoadData;



nPlane=length(confSet.ETL);
frameRepetion=PreMarkPointRepetition+PostMarkPointRepetition; %%Total repepitions of Z series in T series setting.


PVparam.maxFrame=nPlane*frameRepetion;
PVparam.BreakPointFrame=PreMarkPointRepetition*nPlane;
PVparam.maxFrame=nPlane*frameRepetion;

%param for calculate the PSTH heatmap for online analysis
PSTHparam.PreInd=PreMarkPointRepetition-PreSLMCal:PreMarkPointRepetition;
PSTHparam.PostInd=PreMarkPointRepetition+1:PreMarkPointRepetition+PostSLMCal;
PSTHparam.Plot=1;
PSTHparam.SmoothSD=1;
PSTHparam.ColorMap=colorMapPN1;
PSTHparam.Clim=[-400 400];

%param for xml files
XMLparam.ProcessFolder=ProcessFolder;
XMLparam.TotalRounds=confSet.NumTrial;
PointAll=1:size(Pos3Dneed,1);     %All possible MarkPoints for testing
PointsTest=PointAll;              %Initial test Points, this would be updated automatically later


%% Intiate 1st round test
PointAll=1:size(Pos3Dneed,1);     %All possible MarkPoints for testing
PointsTest=PointAll;              %Initial test Points, this would be updated automatically later

SLMTrialInfo=[];                  %Inital response information, automatically updated after each single trial test
SLMTrialMap=[];                   %Inital response map, automatically updated after each single trial test
clear SLMTable;
SLMTable(:,1)=round(1:size(Pos3Dneed,1));
SLMTable(:,2)=NaN;

FileIDrange=[];                   %<-------------------------------------------------------------------------------Edit, BinFile ID range to calculate SLMresponse
minTrialN=1;
% PointTest=PointAll;

%% Keep alternating following 2 lines to do all necessary SLM tests and online analysis
XMLparam.SwitchXMLPostMPFrame=10;                   %%<-----------MarkPoint switching occurs after 10 Repetitions of nplanes of Zseries.
XMLparam.ProcessFolder=ProcessFolder;

XMLparam.RoundID=1;
XMLparam.PointList=[1 2 1 3 2];                     %%<-----------nP, NumOfTestedPoints, nP + 1 = NumOfZseries 
XMLparam.Laser=[repmat(1.5,1,5)] ;                  %%<-----------laser values, could be 1 value of a vector with length of nP



PowerTestPVPar.nPlane=nPlane;            
PowerTestPVPar.ZRepetition=31;                      %%<-----------NumOfRepeition in each Zseries of Tseries in PV
PowerTestPVPar.Ziteration=6;                        %%<-----------NumOfZseries in Tseries in PV
PowerTestPVPar.InterMPRepetition=repmat(PowerTestPVPar.ZRepetition,1,PowerTestPVPar.Ziteration);
frameRepetition=PowerTestPVPar.ZRepetition*PowerTestPVPar.Ziteration;
PowerTestPVPar.maxFrame=nPlane*frameRepetition;
% PowerTestPVPar.BreakPointFrame=PowerTestPVPar.InterMPRepetition(1:end-1)*nPlane;
% PowerTestPVPar.InterMPFrame=[40 60 30 20]*nPlane;
% PowerTestPVPar.TrialMPSwitch=length(PowerTestPVPar.InterMPRepetition)-1;


if sum(abs([length(XMLparam.Laser) length(XMLparam.PointList)]+1-PowerTestPVPar.Ziteration))~=0
   disp('Check whether # Point, Laser levels and Zseries match')
end
iRound=1;

[XMLTable{iRound},FileGenerateInfo(iRound)]=PV_LinkPowerTest_MultiZseries(XMLparam,PowerTestPVPar)
iRound=iRound+1;






%% 
% RandomDelayInterval=[0 1]; %%Random delay is induced after each trial of stimulation.
% PointRepetition=1;  %%Trial Number per each xml MarkPoint stimulation.
PreMarkPointRepetition=34;    %<----------------------------------------------------------------------------------Edit,Frame # before SLM in PV
PostMarkPointRepetition=6;   %<----------------------------------------------------------------------------------Edit,Frame # after SLM in PV
PreSLMCal=15;                 %<----------------------------------------------------------------------------------Edit,Frame # before SLM to calculate baseline map
PostSLMCal=3;                 %<----------------------------------------------------------------------------------Edit,Frame # before SLM to calculate responsive map
% ROIparam.LaserPower=confSet.UncagingLaserPower;  %Laser power to test, using all possible power levels
% ROIparam.LaserPower=confSet.UncagingLaserPower([2 3]);  %<--------------------------------------------------------Edit,It is not necessary to test all possible power levels



ROIparam.LaserPower=confSet.UncagingLaserPower;
ROIparam.min_merged_region_size=30;
ROIparam.threshold_percentage=0.25;
ROIparam.PointsTest=PointsTest;
ROIparam.Clim=[-400 400];
ROIparam.max_distance=8;


SLMTestParam.TerminalTrialN=4;    %<-------------------------------------------------------------------------------Edit, Trials # to define SLM responsive cells
SLMTestParam.ExcludeTrialN=1;     %<-------------------------------------------------------------------------------Edit, Trials # to define Non-SLM responsive cells
SLMTestParam.AllLaserPower=confSet.UncagingLaserPower;


% %param for ROI neighbourhood to determine wether there is SLM response.
% ROIparam.TotalSLMPos3D=Pos3Dneed;    %%such that ROIparam.PointAll=1:size(Pos3Dneed,1)
% ROIparam.PointAll=PointAll;
% ROIparam.PlaneZ=PlaneZ;
% ROIparam.CellSize=20;                %%normal neuron diameter by um;        
% ROIparam.threshold_percentage=0.3;   %%thereshold to define responsive fields SLM responsive heatmap: percentage*Peak rate
% ROIparam.thNum=15;                   %%Minimal single responsive field by pixels
% ROIparam.max_distance=ceil(ROIparam.CellSize*2/3/umPerPixel);  %% 2/3 diameter of a cell by pixel as maximal response region-SLM center distance
% ROIparam.min_region_size=10;
% ROIparam.PeakTh=200;
% ROIparam.min_merged_region_size=20;  %%Minimal total size of responsive fields by pixels
% ROIparam.contourMethod='perim';      %%Method to detect ROI boader of responsive fields
% ROIparam.NeighbourHfWidthPixel=20;   %%PixFromMedCenter: Number of pixels from the median center of each ROI to get the ROI neighborhood.
% ROIparam.umPerPixel=mean([yaml.umPerlPixelX yaml.umPerlPixelY]);
% ROIparam.Colormap=colorMapPN1;                  
% % ROIparam.LaserPower=confSet.UncagingLaserPower;   





% FileIDrange=[1;400];             %<------------------------------------------------------------------------------Edit, BinFile ID range to calculate SLMresponse
[SLMRes,sampleN]=SLMResponseROIMap(SLMTrialMap,SLMTrialInfo,ROIparam,minTrialN,SumDataFolder,FileIDrange);


ROIparam.LaserPower=confSet.UncagingLaserPower([1 2 3]);  %<-------------------------------------------------------Edit,It is not necessary to test all possible power levels
[UpdateXml, SLMTable, PointsTest, XMLparam, InfoListByLaser]=step3Fun_NextSLMtest(SLMRes,sampleN,ROIparam,XMLparam,SLMTestParam,SLMTable);

% PointsTest=[1 2 3 4];
% XMLparam.Laser=1.5;
% PSTHparam.Plot=1;

%% Save results and generate final MarkPoints and Functional groups
SLMTableOrigin=SLMTable;
% PostSLMTable;
FalsePositiveID=input('Mannual correction: index in SLMTable with false positive error: ');
SLMTable(FalsePositiveID,:)=nan;

save([SumDataFolder 'SLMResponseTable.mat'],'SLMTableOrigin','SLMTable','ROIparam','SLMRes','sampleN','SLMTestParam','SLMIncludedIndFromIscell','FunScore','yaml','Cellstat');

% load([SumDataFolder 'SLMResponseTable.mat'],'SLMTable','ROIparam','SLMRes','sampleN','SLMTestParam','SLMIncludedIndFromIscell','FunScore');

SMLTablePowerPV=xmlPower2PVpower(SLMTable(:,2));
refPVpower=max(SMLTablePowerPV);

CellPerGroup=10;
[Group, FinalPos3D, FinalCellstat, FinalFunScore, confSetFinal] = SLMWeightsAssignToFunGroups(FunScore, CellPerGroup, Pos3Dneed, Cellstat, SLMIncludedIndFromIscell, SLMTable, NonTargets, refPVpower, confSet);
save([ProcessFolder 'SLMFunGroup.mat'],'Group','FinalPos3D','FinalCellstat','FinalFunScore','confSetFinal','SLMTableOrigin','SLMTable','ROIparam','SLMRes','sampleN','SLMTestParam','SLMIncludedIndFromIscell','FunScore','yaml','Cellstat');
XYZtoMarkPointFunGroup(ProcessFolder,FinalPos3D,Group,yaml,confSetFinal);

