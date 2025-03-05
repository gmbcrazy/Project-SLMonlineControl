

%% Load Data
clear all
% ProcessFolder='F:\LuSLMOnlineTest\04222024\SingleP\30PixelFromEdgeExc\';
load('C:\Users\User\Project-SLMonlineControl\subfun\Color\colorMapPN3.mat');
ConfigFolder='C:\Users\User\Project-SLMonlineControl\config\';
WorkingFolder='E:\LuSLMOnlineTest\NoAnimalTest\01232025\';%<---------------------------Edit, Data folder
PreDefTseriesFolder=[ConfigFolder 'PreGenerateTseriesMultiZ\'];
PreDefTmat='SpontBeh5T_Z11Frame550.mat';                  %<---------------------------Edit, TseriesPreDefined
PreDefFolder=[PreDefTseriesFolder PreDefTmat(1:end-4) '\'];


ConfigFile='SLMsetting.yml';%<---------------------------------------------------------Edit, configuration file
[~,~,~,CaData,CaDataPlane,stat,yaml,confSet]=ROIToXYZ(ConfigFolder);

confSet = ReadYaml([ConfigFolder '\' ConfigFile]);

umPerPixel=mean([yaml.umPerlPixelX yaml.umPerlPixelY]);
ProcessFolder=[WorkingFolder 'SingleP\Top12SpeedStimEdgeExc\'];%<----------------------Edit, Data folder
SumDataFolder=[ProcessFolder '\DataSum\'];
mkdir(SumDataFolder)
DataLogFolder=[ProcessFolder 'DataLog\'];
mkdir(DataLogFolder)
DataFolder=[ProcessFolder 'Data\'];
% mkdir(DataLogFolder)

load([ProcessFolder 'SLMFunGroup.mat'],'Group','FinalPos3D','FinalCellstat','FinalFunScore','confSetFinal','SLMTableOrigin','SLMTable','ROIparam','SLMRes','sampleN','SLMTestParam','SLMIncludedIndFromIscell','FunScore','yaml','Cellstat');
load([PreDefTseriesFolder PreDefTmat],'TSeriesBrukerTBL');              

TSeriesENVFile=dir([PreDefFolder '*.env']);




% PreMarkPointRepetition=25;    %<----------------------------------------------------------------------------------Edit,Frame # before SLM in PV
% PostMarkPointRepetition=10;   %<----------------------------------------------------------------------------------Edit,Frame # after SLM in PV
PreSLMCal=15;                   %<----------------------------------------------------------------------------------Edit,Frame # before SLM to calculate baseline map
PostSLMCal=3;                   %<----------------------------------------------------------------------------------Edit,Frame # before SLM to calculate responsive map
nPlane=length(confSet.ETL);




XMLparam.ProcessFolder=ProcessFolder;
XMLparam.TotalRounds=confSet.NumTrial;
% PointAll=1:size(Pos3Dneed,1);
XMLparam.SwitchXMLPostMPFrame=10;
XMLparam.ProcessFolder=ProcessFolder;
%param for calculate the PSTH heatmap for online analysis
% PSTHparam.PreInd=PreMarkPointRepetition-PreSLMCal:PreMarkPointRepetition;
% PSTHparam.PostInd=PreMarkPointRepetition+1:PreMarkPointRepetition+PostSLMCal;
% PSTHparam.Plot=1;
PSTHparam.SmoothSD=1;
PSTHparam.ColorMap=colorMapPN1;
PSTHparam.Clim=[-400 400];
PSTHparam.PreSLMCal=15;
PSTHparam.PostSLMCal=3;
PSTHparam.MPFrameJump=2;





%% 
numGPUs=0;      %%Do not use GPU, assume in general the aquisition PC has no GPU. 
FileType=2;   %Choose a specific bin file as reference for motion correction
% ProcessFolder='E:\LuSLMOnlineTest\SL0777-Ai203\12122024\'
% RefFile=[DataFolder 'TSeries-12132024-1247-023.bin'];
% [RegOps, RegImg] = LoadRegRefFile(RefFile, FileType, numGPUs, [512,512,3,30]);
% 
% FileType=0;   %Choose a pre-recorded multi-tif files for motion correction
% RefFile=[];
% RefFile='E:\LuSLMOnlineTest\SL0777-Ai203\12122024\TSeries-12122024-0938-000\';
% 
% [RegOps, RegImg] = LoadRegRefFile(RefFile, FileType);
% MultiMatrix3DHeatmap(RegImg)
numGPUs=0;      %%Do not use GPU, assume in general the aquisition PC has no GPU. 
FileType=2;   %Choose a specific bin file as reference for motion correction
% ProcessFolder='E:\LuSLMOnlineTest\SL0777-Ai203\12122024\'
% RefFile=[DataFolder 'TSeries-12132024-1247-023.bin'];
% [RegOps, RegImg] = LoadRegRefFile(RefFile, FileType, numGPUs, [512,512,3,30]);
% 
FileType=0;   %Choose a pre-recorded multi-tif files for motion correction
RefFile=[];
RefFile=[WorkingFolder 'RegRef1\'];
% 
[RegOps, RegImg] = LoadRegRefFile(RefFile, FileType,numGPUs);

% FileType=1;   %Choose suite2p folder, using ops.meanImg for motion correction
% RefFile=[WorkingFolder 'suite2p\'];
% [RegOps, RegImg] = LoadRegRefFile(RefFile, FileType,numGPUs);

XMLparam.DoRegistration=1;
XMLparam.RegRefOps=RegOps;
XMLparam.RegRefImg=RegImg;  
XMLparam.RegRefSource=RefFile;  

disp('Referrence Img for Motion correction updated')


%%
XMLTable=[];
PSTHmap=[];
CountExp=1;
TotalGroupIDs=[1 2 3];   %% All possible Functional Group IDs.

XMLparam.LoadGPL=1;

% pause(10)
nPlane=3;
for TseriesID=3:3
PVparam=BrukerTBLtoPVparm(TSeriesBrukerTBL{TseriesID},nPlane);   %%Update Tseries
LoadTSeriestoBruker(TSeriesENVFile(TseriesID))                   %%Update PVparam with Current Tseries
% [~,~]=PV_LinkExcuteXMLFunGroup(XMLparam,PVparam);
pause(4.5);
'Ready'
[~,~]=PV_LinkExcuteDefTseries_XMLFunGroup(XMLparam,PVparam)
pause(6)
end

PVparam.TseriesTBL

