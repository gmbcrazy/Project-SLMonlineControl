%% Load Data
clear all
% ProcessFolder='F:\LuSLMOnlineTest\04222024\SingleP\30PixelFromEdgeExc\';
% load('C:\Users\zhangl33\Projects\GenMatCode\Plotfun\Color\colorMapPN3.mat');
ConfigFolder='C:\Users\User\Project-SLMonlineControl\config\';
ConfigFile='SLMsettingG7.yml';%<----------------------------------------------------------------------------------Edit, configuration file
[~,~,~,CaData,CaDataPlane,stat,yaml,confSet]=ROIToXYZ(ConfigFolder);
umPerPixel=mean([yaml.umPerlPixelX yaml.umPerlPixely]);
load('C:\Users\User\Project-SLMonlineControl\subfun\Color\colorMapPN3.mat');
ProcessFolder='F:\LuSLMOnlineTest\SL0702-G7-Ai203\09192024\SingleP\Top15SpeedStimEdgeExc\';%<---------------------Edit, Data folder

step4_SubStep1_LoadData;



%% 
% RandomDelayInterval=[0 1]; %%Random delay is induced after each trial of stimulation.
% PointRepetition=1;  %%Trial Number per each xml MarkPoint stimulation.
PreMarkPointRepetition=40;    %<----------------------------------------------------------------------------------Edit,Frame # before SLM in PV
PostMarkPointRepetition=20;   %<----------------------------------------------------------------------------------Edit,Frame # after SLM in PV
PreSLMCal=20;                 %<----------------------------------------------------------------------------------Edit,Frame # before SLM to calculate baseline map
PostSLMCal=3;                 %<----------------------------------------------------------------------------------Edit,Frame # before SLM to calculate responsive map
% ROIparam.LaserPower=confSet.UncagingLaserPower;  %Laser power to test, using all possible power levels
% ROIparam.LaserPower=confSet.UncagingLaserPower([2 3]);  %<--------------------------------------------------------Edit,It is not necessary to test all possible power levels

step4_SubStep2_PreparingParam;




%% Intiate 1st round test
XMLparam.Laser=1.5;               %<-------------------------------------------------------------------------------Edit, starting laser power to test    
XMLparam.RoundID=1;               %starting round
PointAll=1:size(Pos3Dneed,1);     %All possible MarkPoints for testing
PointsTest=PointAll;              %Initial test Points, this would be updated automatically later
SLMTrialInfo=[];                  %Inital response information, automatically updated after each single trial test
SLMTrialMap=[];                   %Inital response map, automatically updated after each single trial test
clear SLMTable;
SLMTable(:,1)=round(1:size(Pos3Dneed,1));
SLMTable(:,2)=NaN;




%% Keep alternating following 2 lines to do all necessary SLM tests and online analysis
step4_SubStep3_InitiateTest

ROIparam.LaserPower=confSet.UncagingLaserPower([2 3]);  %<--------------------------------------------------------Edit,It is not necessary to test all possible power levels
[UpdateXml, SLMTable, PointsTest, XMLparam, InfoListByLaser]=step3Fun_NextSLMtest(SLMRes,sampleN,ROIparam,XMLparam,SLMTestParam,SLMTable);

