

%% Load Data
clear all
% ProcessFolder='F:\LuSLMOnlineTest\04222024\SingleP\30PixelFromEdgeExc\';
load('C:\Users\User\Project-SLMonlineControl\subfun\Color\colorMapPN3.mat');
ConfigFolder='C:\Users\User\Project-SLMonlineControl\config\';
ConfigFile='SLMsetting.yml';%<----------------------------------------------------------------------------------Edit, configuration file
[~,~,~,CaData,CaDataPlane,stat,yaml,confSet]=ROIToXYZ(ConfigFolder);
umPerPixel=mean([yaml.umPerlPixelX yaml.umPerlPixelY]);
ProcessFolder='F:\LuSLMOnlineTest\SL0541-Emx1Ai96\10032024\SingleP\Top13SpeedStimEdgeExc\';%<----------------------Edit, Data folder


% PreMarkPointRepetition=25;    %<----------------------------------------------------------------------------------Edit,Frame # before SLM in PV
% PostMarkPointRepetition=10;   %<----------------------------------------------------------------------------------Edit,Frame # after SLM in PV
PreSLMCal=15;                 %<----------------------------------------------------------------------------------Edit,Frame # before SLM to calculate baseline map
PostSLMCal=3;                 %<----------------------------------------------------------------------------------Edit,Frame # before SLM to calculate responsive map

nPlane=length(confSet.ETL);

PVparam.InterMPRepetion=[40 60 30 20];
frameRepetion=sum(PVparam.InterMPRepetion); %%Total repepitions of Z series in T series setting;


PVparam.maxFrame=nPlane*frameRepetion;
PVparam.BreakPointFrame=PVparam.InterMPRepetion(1:end-1)*nPlane;
% PVparam.InterMPFrame=[40 60 30 20]*nPlane;
PVparam.TrialMPSwitch=length(PVparam.InterMPRepetion)-1;
PVparam.nPlane=nPlane;




%param for calculate the PSTH heatmap for online analysis
% PSTHparam.PreInd=PreMarkPointRepetition-PreSLMCal:PreMarkPointRepetition;
% PSTHparam.PostInd=PreMarkPointRepetition+1:PreMarkPointRepetition+PostSLMCal;
% PSTHparam.Plot=1;
% PSTHparam.SmoothSD=1;
% PSTHparam.ColorMap=colorMapPN1;
% PSTHparam.Clim=[-400 400];

%param for xml files
% XMLparam.Laser=1.5;               %<-------------------------------------------------------------------------------Edit, starting laser power to test    
% XMLparam.RoundID=1;               %starting round

XMLparam.ProcessFolder=ProcessFolder;
XMLparam.TotalRounds=confSet.NumTrial;
% PointAll=1:size(Pos3Dneed,1);


%param for ROI neighbourhood to determine wether there is SLM response.
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
% ROIparam.LaserPower=confSet.UncagingLaserPower;   

XMLparam.ShamPossibility=0.4;
XMLparam.SwitchXMLPostMPFrame=10;
XMLparam.ProcessFolder=ProcessFolder;



% PreMarkPointRepetition=[40 60 50];
% PostMarkPointRepetition=20;




%%
PV_LinkExcuteXMLFunGroup(XMLparam,PVparam,confSet)

